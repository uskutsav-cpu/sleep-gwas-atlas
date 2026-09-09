"""Transparent evidence aggregation and competitive locus-level cell annotation."""
from __future__ import annotations
import collections
import math
from pathlib import Path
from .io import read_tsv, require, write_tsv
from .stats import bh, stratified_enrichment


def integrate(loci_file, molecular_file, out):
    """Do not convert proximity/expression into a causal-gene claim.

    Inputs are audited, normalized evidence exports (see docs/INPUTS.md), not
    arbitrary gene lists. A negative/ambiguous locus is retained, never dropped.
    """
    loci=list(read_tsv(loci_file,["pair_id","locus_id","status","trait_h4","pip_max"]))
    molecular=list(read_tsv(molecular_file,["pair_id","locus_id","gene_id","evidence_type",
                                         "molecular_h4","tissue","cell_type","source"]))
    keys={(r["pair_id"],r["locus_id"]) for r in loci}
    require(len(keys)==len(loci),"Duplicate locus keys")
    group=collections.defaultdict(list)
    allowed={"eqtl_coloc","sqtl_coloc","coding","chromatin_link","expression","nearest_gene"}
    for r in molecular:
        key=r["pair_id"],r["locus_id"]
        require(key in keys,"Molecular evidence not linked to an analyzed locus")
        require(r["evidence_type"] in allowed and r["source"],"Unknown/unsourced molecular evidence")
        group[key].append(r)
    result=[]
    def probability(v):
        if v in {"","NA"}: return None
        p=float(v); require(math.isfinite(p) and 0<=p<=1,"Invalid evidence probability"); return p
    for locus in loci:
        key=locus["pair_id"],locus["locus_id"]
        h4=probability(locus["trait_h4"]); pip=probability(locus["pip_max"])
        rows=group[key]
        genes=sorted({r["gene_id"] for r in rows}) or [""]
        for gene in genes:
            gr=[r for r in rows if r["gene_id"]==gene]
            types={r["evidence_type"] for r in gr}
            qtl=[]
            for r in gr:
                h=probability(r["molecular_h4"])
                if r["evidence_type"] in {"eqtl_coloc","sqtl_coloc"} and h is not None:
                    qtl.append(h)
            # Display categories, not inferential guarantees or an opaque weighted score.
            if locus["status"]!="PASS": tier="UNRESOLVED_OR_BLOCKED"
            elif h4 is None or h4<.8: tier="REGIONAL_OR_DISTINCT_SIGNAL"
            elif qtl and max(qtl)>=.8: tier="SHARED_SIGNAL_WITH_MOLECULAR_SUPPORT"
            elif "coding" in types: tier="SHARED_SIGNAL_WITH_CODING_ANNOTATION"
            else: tier="SHARED_SIGNAL_GENE_UNRESOLVED"
            result.append(dict(pair_id=key[0],locus_id=key[1],gene_id=gene,tier=tier,
                    trait_h4="NA" if h4 is None else h4,pip_max="NA" if pip is None else pip,
                    molecular_h4_max=max(qtl) if qtl else "NA",evidence_types=";".join(sorted(types)),
                    tissues=";".join(sorted({r["tissue"] for r in gr if r["tissue"]})),
                    cell_types=";".join(sorted({r["cell_type"] for r in gr if r["cell_type"]})),
                    sources=";".join(sorted({r["source"] for r in gr}))))
    fields=["pair_id","locus_id","gene_id","tier","trait_h4","pip_max","molecular_h4_max",
            "evidence_types","tissues","cell_types","sources"]
    write_tsv(out,fields,result)
    return result


def cell_enrichment(table, out, *, permutations=10000, seed=20260908):
    """Test precomputed cell annotation scores using independent matched loci.

    Not a raw-scRNA analysis and not a replacement for S-LDSC/MAGMA covariates.
    The unit table must be prepared with reviewed LD blocks and matched control
    loci; the primary tests across ALL pair x cell combinations form one family.
    """
    rows=list(read_tsv(table,["pair_id","cell_type","locus_id","score","selected","stratum"]))
    groups=collections.defaultdict(list)
    seen=set()
    for r in rows:
        key=r["pair_id"],r["cell_type"],r["locus_id"]
        require(key not in seen,"Locus counted more than once in a cell test")
        require(r["selected"] in {"0","1"},"selected must be 0 or 1")
        seen.add(key); groups[key[:2]].append(r)
    results=[]
    for (pair,cell),rs in sorted(groups.items()):
        ans=stratified_enrichment([float(r["score"]) for r in rs],
                [r["selected"]=="1" for r in rs],[r["stratum"] for r in rs],
                permutations=permutations,seed=seed)
        results.append(dict(pair_id=pair,cell_type=cell,**ans))
    require(results,"No enrichment tests")
    q=bh([r["p"] for r in results])
    for r,v in zip(results,q):r["family_fdr"]=v
    write_tsv(out,list(results[0]),results)
    return results


def cross_disorder(evidence_file, out):
    """Count distinct loci, not SNPs or annotation rows; unknown != absent."""
    rows=list(read_tsv(evidence_file,["pair_id","locus_id","gene_id","tier"]))
    groups=collections.defaultdict(lambda:collections.defaultdict(set))
    for r in rows:
        groups[r["pair_id"]][r["tier"]].add(r["locus_id"])
    results=[]
    for pair,tiers in sorted(groups.items()):
        for tier,loci in sorted(tiers.items()):
            results.append(dict(pair_id=pair,tier=tier,n_independent_blocks=len(loci)))
    write_tsv(out,["pair_id","tier","n_independent_blocks"],results)
    return results
