#!/usr/bin/env python3
"""Shared fail-closed parsers for checksum-pinned pathway source releases."""
from __future__ import annotations

import gzip
import hashlib
import math
import re
import zipfile
from collections import defaultdict
from pathlib import Path


GENE_ATTRIBUTE = re.compile(r'(?:^|;\s*)([A-Za-z_]+) "([^"]*)"')
GO_NAMESPACES = {"biological_process", "molecular_function", "cellular_component"}


def gtf_attributes(text: str) -> dict[str, str]:
    return {key: value for key, value in GENE_ATTRIBUTE.findall(text)}


def load_unique_gencode_symbols(path: Path) -> tuple[dict[str, str], dict[str, int]]:
    """Return exact unique symbol->versionless Ensembl mappings from gene rows."""
    symbols: dict[str, set[str]] = defaultdict(set)
    gene_rows = 0
    with path.open(encoding="utf-8") as handle:
        for line_number, line in enumerate(handle, start=1):
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 9:
                raise ValueError(f"malformed GENCODE GTF row at line {line_number}")
            if fields[2] != "gene":
                continue
            attributes = gtf_attributes(fields[8])
            gene_id = attributes.get("gene_id", "").split(".", 1)[0]
            symbol = attributes.get("gene_name", "")
            if not re.fullmatch(r"ENSG[0-9]{11}", gene_id) or not symbol:
                raise ValueError(f"GENCODE gene row lacks a valid ID or symbol at line {line_number}")
            symbols[symbol].add(gene_id)
            gene_rows += 1
    unique = {symbol: next(iter(ids)) for symbol, ids in symbols.items() if len(ids) == 1}
    return unique, {
        "gene_rows": gene_rows,
        "distinct_symbols": len(symbols),
        "unique_symbols": len(unique),
        "ambiguous_symbols": sum(len(ids) > 1 for ids in symbols.values()),
    }


def load_reactome_sets(
    path: Path, symbol_map: dict[str, str], minimum_size: int, maximum_size: int,
) -> tuple[dict[str, tuple[str, set[str]]], dict[str, int | str]]:
    with zipfile.ZipFile(path) as archive:
        files = [info for info in archive.infolist() if not info.is_dir()]
        if [info.filename for info in files] != ["ReactomePathways.gmt"]:
            raise ValueError("Reactome archive does not contain exactly ReactomePathways.gmt")
        payload = archive.read(files[0])
    sets: dict[str, tuple[str, set[str]]] = {}
    source_genes: set[str] = set()
    mapped_genes: set[str] = set()
    eligible_gene_total = 0
    for line_number, line in enumerate(payload.decode("utf-8").splitlines(), start=1):
        fields = line.split("\t")
        if len(fields) < 3 or not fields[0] or not re.fullmatch(r"R-HSA-[0-9]+", fields[1]):
            raise ValueError(f"invalid Reactome GMT row at line {line_number}")
        name, identity = fields[:2]
        if identity in sets:
            raise ValueError(f"duplicate Reactome pathway ID: {identity}")
        source_symbols = set(fields[2:])
        source_genes.update(source_symbols)
        genes = {symbol_map[symbol] for symbol in source_symbols if symbol in symbol_map}
        mapped_genes.update(genes)
        if minimum_size <= len(genes) <= maximum_size:
            sets[identity] = (name, genes)
            eligible_gene_total += len(genes)
    return sets, {
        "archive_member": "ReactomePathways.gmt",
        "archive_member_bytes": len(payload),
        "archive_member_sha256": hashlib.sha256(payload).hexdigest(),
        "source_pathways": len(payload.splitlines()),
        "source_symbols": len(source_genes),
        "mapped_genes": len(mapped_genes),
        "eligible_sets": len(sets),
        "eligible_gene_memberships": eligible_gene_total,
    }


def parse_go_ontology(path: Path) -> tuple[dict[str, dict[str, object]], dict[str, str]]:
    text = path.read_text(encoding="utf-8")
    metadata: dict[str, str] = {}
    for line in text.split("\n[Term]\n", 1)[0].splitlines():
        if line.startswith("data-version: "):
            metadata["data_version"] = line.split(": ", 1)[1]
    terms: dict[str, dict[str, object]] = {}
    for block in text.split("\n[Term]\n")[1:]:
        identity = name = namespace = ""
        obsolete = False
        parents: set[str] = set()
        for line in block.splitlines():
            if line.startswith("["):
                break
            if line.startswith("id: "):
                identity = line[4:]
            elif line.startswith("name: "):
                name = line[6:]
            elif line.startswith("namespace: "):
                namespace = line[11:]
            elif line == "is_obsolete: true":
                obsolete = True
            elif line.startswith("is_a: "):
                parents.add(line.split()[1])
            elif line.startswith("relationship: part_of "):
                parents.add(line.split()[2])
        if identity and not obsolete and namespace in GO_NAMESPACES:
            terms[identity] = {"name": name, "namespace": namespace, "parents": parents}
    for term in terms.values():
        term["parents"] = {parent for parent in term["parents"] if parent in terms}
    return terms, metadata


def go_ancestors(identity: str, terms: dict[str, dict[str, object]], cache: dict[str, set[str]], active: set[str]) -> set[str]:
    if identity in cache:
        return cache[identity]
    if identity in active:
        raise ValueError(f"GO basic graph contains a cycle involving {identity}")
    active.add(identity)
    ancestors = {identity}
    for parent in terms[identity]["parents"]:
        ancestors.update(go_ancestors(str(parent), terms, cache, active))
    active.remove(identity)
    cache[identity] = ancestors
    return ancestors


def load_go_sets(
    ontology_path: Path, annotation_path: Path, symbol_map: dict[str, str],
    minimum_size: int, maximum_size: int,
) -> tuple[dict[str, tuple[str, set[str]]], dict[str, int | str]]:
    terms, metadata = parse_go_ontology(ontology_path)
    direct: dict[str, set[str]] = defaultdict(set)
    rows = excluded_not = excluded_unmapped = 0
    comments: list[str] = []
    with gzip.open(annotation_path, "rt", encoding="utf-8") as handle:
        for line_number, line in enumerate(handle, start=1):
            if line.startswith("!"):
                comments.append(line.rstrip("\n"))
                continue
            fields = line.rstrip("\n").split("\t")
            rows += 1
            if (
                len(fields) != 17 or fields[12].split("|", 1)[0] != "taxon:9606"
                or fields[4] not in terms
            ):
                raise ValueError(f"invalid human GO GAF row at line {line_number}")
            if "NOT" in fields[3].split("|"):
                excluded_not += 1
                continue
            symbol = fields[2]
            if symbol not in symbol_map:
                excluded_unmapped += 1
                continue
            direct[fields[4]].add(symbol_map[symbol])
    propagated: dict[str, set[str]] = defaultdict(set)
    ancestor_cache: dict[str, set[str]] = {}
    for identity, genes in direct.items():
        for ancestor in go_ancestors(identity, terms, ancestor_cache, set()):
            propagated[ancestor].update(genes)
    sets = {
        identity: (str(terms[identity]["name"]), genes)
        for identity, genes in propagated.items()
        if minimum_size <= len(genes) <= maximum_size
    }
    generated = next((line.split(": ", 1)[1] for line in comments if line.startswith("!date-generated: ")), "")
    go_version = next((line.split(": ", 1)[1] for line in comments if line.startswith("!go-version: ")), "")
    return sets, {
        "ontology_data_version": metadata.get("data_version", ""),
        "active_ontology_terms": len(terms),
        "gaf_rows": rows,
        "gaf_date_generated": generated,
        "gaf_go_version": go_version,
        "excluded_not_rows": excluded_not,
        "excluded_unmapped_rows": excluded_unmapped,
        "direct_annotated_terms": len(direct),
        "propagated_terms": len(propagated),
        "mapped_genes": len(set().union(*propagated.values())) if propagated else 0,
        "eligible_sets": len(sets),
        "eligible_gene_memberships": sum(len(genes) for _, genes in sets.values()),
    }


def hypergeometric_right_tail(overlap: int, universe: int, set_size: int, selected: int) -> float:
    """P[X >= overlap] for X~Hypergeometric(universe, set_size, selected)."""
    if not (0 <= overlap <= min(set_size, selected) <= universe and set_size <= universe):
        raise ValueError("invalid hypergeometric parameters")
    if overlap == 0:
        return 1.0
    upper = min(set_size, selected)
    logs = [
        math.lgamma(set_size + 1) - math.lgamma(value + 1) - math.lgamma(set_size - value + 1)
        + math.lgamma(universe - set_size + 1)
        - math.lgamma(selected - value + 1)
        - math.lgamma(universe - set_size - selected + value + 1)
        - math.lgamma(universe + 1) + math.lgamma(selected + 1)
        + math.lgamma(universe - selected + 1)
        for value in range(overlap, upper + 1)
        if 0 <= selected - value <= universe - set_size
    ]
    maximum = max(logs)
    return min(1.0, math.exp(maximum) * sum(math.exp(value - maximum) for value in logs))
