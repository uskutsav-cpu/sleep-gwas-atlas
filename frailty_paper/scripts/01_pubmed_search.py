#!/usr/bin/env python3

from __future__ import annotations

import argparse
import csv
import datetime as dt
import json
import os
import time
from pathlib import Path

import requests
from lxml import etree


BASE = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"

MAX_PUBMED = 9990
BATCH = 200


def req(session, endpoint, params, api_key=None):

    params = dict(params)

    if api_key:
        params["api_key"] = api_key

    last_error = None

    for attempt in range(6):

        try:
            r = session.get(
                f"{BASE}/{endpoint}",
                params=params,
                timeout=90,
            )
        except requests.exceptions.RequestException as exc:
            last_error = exc
            if attempt == 5:
                raise
            time.sleep(2 ** attempt)
            continue

        if r.status_code == 200:

            time.sleep(
                0.11 if api_key else 0.36
            )

            return r

        if r.status_code in {
            429,
            500,
            502,
            503,
            504,
        }:

            time.sleep(2 ** attempt)

            continue

        r.raise_for_status()

    raise RuntimeError(
        f"NCBI request failed after retries: {endpoint}"
    ) from last_error


def text(el):

    if el is None:
        return ""

    return "".join(
        el.itertext()
    ).strip()


def parse_article(article, query_id):

    if article.tag == "PubmedArticle":
        record = article.find("MedlineCitation")
        content = record.find("Article") if record is not None else None
        pmid = text(record.find("PMID")) if record is not None else ""
        title = text(content.find("ArticleTitle")) if content is not None else ""
        journal = text(content.find("Journal/Title")) if content is not None else ""
        pubdate = content.find("Journal/JournalIssue/PubDate") if content is not None else None
        abstract_nodes = content.findall("Abstract/AbstractText") if content is not None else []
    elif article.tag == "PubmedBookArticle":
        record = article.find("BookDocument")
        content = record.find("Book") if record is not None else None
        pmid = text(record.find("PMID")) if record is not None else ""
        title = (
            text(record.find("ArticleTitle")) or text(content.find("BookTitle"))
            if content is not None else text(record.find("ArticleTitle"))
        ) if record is not None else ""
        journal = (
            text(content.find("CollectionTitle")) or text(content.find("BookTitle"))
            if content is not None else ""
        )
        pubdate = content.find("PubDate") if content is not None else None
        abstract_nodes = record.findall("Abstract/AbstractText") if record is not None else []
    else:
        raise ValueError(f"Unsupported PubMed record type: {article.tag}")

    year = (
        text(pubdate.find("Year")) or text(pubdate.find("MedlineDate"))[:4]
        if pubdate is not None else ""
    )
    abstract = " ".join(text(node) for node in abstract_nodes)

    doi = ""

    pubmed_data = article.find(
        "PubmedData"
    )

    if pubmed_data is not None:

        for aid in pubmed_data.findall(
            "ArticleIdList/ArticleId"
        ):

            if aid.get("IdType") == "doi":

                doi = text(aid)

                break

    return {

        "pmid": pmid,
        "doi": doi,
        "year": year,
        "journal": journal,
        "title": title,
        "abstract": abstract,
        "query_id": query_id,

    }


def dated_term(
    base_query,
    start,
    end,
):

    return (

        f'({base_query}) AND '

        f'("{start:%Y/%m/%d}"'
        f'[Date - Publication] : '

        f'"{end:%Y/%m/%d}"'
        f'[Date - Publication])'

    )


def esearch(
    session,
    term,
    email,
    api_key,
    history=False,
    retmax=0,
):

    params = {

        "db": "pubmed",
        "term": term,
        "retmode": "json",
        "retmax": retmax,

    }

    if history:

        params["usehistory"] = "y"

    if email:

        params["email"] = email

    for attempt in range(6):
        payload = req(
            session,
            "esearch.fcgi",
            params,
            api_key,
        ).json()
        result = payload.get("esearchresult", {})
        if isinstance(result, dict) and "count" in result:
            return result

        error = result.get("ERROR", "missing esearchresult.count") if isinstance(result, dict) else str(result)
        transient = "temporarily unavailable" in error.lower() or "search backend failed" in error.lower()
        if transient and attempt < 5:
            time.sleep(2 ** attempt)
            continue
        raise RuntimeError(f"PubMed ESearch failed: {error}")

    raise RuntimeError("PubMed ESearch failed after retries")


def parse_batch(xml_bytes, query_id, expected_count):

    root = etree.fromstring(xml_bytes)
    articles = [
        article
        for article in root
        if article.tag in {"PubmedArticle", "PubmedBookArticle"}
    ]

    if len(articles) != expected_count:
        raise ValueError(
            f"Expected {expected_count} PubMed articles; found {len(articles)}."
        )

    return [parse_article(article, query_id) for article in articles]


def split_segments(
    session,
    base_query,
    start,
    end,
    email,
    api_key,
):

    term = dated_term(
        base_query,
        start,
        end,
    )

    result = esearch(
        session,
        term,
        email,
        api_key,
        history=False,
    )

    count = int(
        result["count"]
    )

    if count <= MAX_PUBMED:

        yield (
            start,
            end,
            count,
        )

        return

    if start >= end:

        raise RuntimeError(
            f"More than {MAX_PUBMED} "
            f"records occur on {start}. "
            "Query needs refinement."
        )

    span = (
        end - start
    ).days

    mid = (
        start
        +
        dt.timedelta(
            days=span // 2
        )
    )

    yield from split_segments(

        session,
        base_query,
        start,
        mid,
        email,
        api_key,

    )

    yield from split_segments(

        session,
        base_query,
        mid + dt.timedelta(days=1),
        end,
        email,
        api_key,

    )


def fetch_segment(
    session,
    qid,
    base_query,
    start,
    end,
    qdir,
    email,
    api_key,
    all_rows,
):

    term = dated_term(
        base_query,
        start,
        end,
    )

    result = esearch(
        session,
        term,
        email,
        api_key,
        history=True,
    )

    count = int(
        result["count"]
    )

    if count == 0:

        return 0

    if count > MAX_PUBMED:

        raise RuntimeError(
            "Internal error: "
            f"unsplit segment has "
            f"{count} records."
        )

    id_result = esearch(
        session,
        term,
        email,
        api_key,
        history=False,
        retmax=count,
    )
    segment_ids = [str(value) for value in id_result.get("idlist", [])]
    if len(segment_ids) != count or len(set(segment_ids)) != count:
        raise RuntimeError(
            f"ESearch ID snapshot mismatch for {qid} {start}..{end}: "
            f"count={count}, unique_ids={len(set(segment_ids))}"
        )

    webenv = result["webenv"]
    qkey = result["querykey"]

    segdir = (
        qdir
        /
        (
            "segment_"
            f"{start.isoformat()}_"
            f"{end.isoformat()}"
        )
    )

    segdir.mkdir(
        parents=True,
        exist_ok=True,
    )

    (
        segdir
        /
        "query.txt"
    ).write_text(
        term + "\n",
        encoding="utf-8",
    )

    for retstart in range(
        0,
        count,
        BATCH,
    ):

        xml_path = (
            segdir
            /
            f"batch_{retstart:06d}.xml"
        )
        expected_count = min(BATCH, count - retstart)
        expected_ids = segment_ids[retstart:retstart + expected_count]

        if xml_path.is_file():
            try:
                cached_rows = parse_batch(
                    xml_path.read_bytes(),
                    qid,
                    expected_count,
                )
            except (OSError, ValueError, etree.XMLSyntaxError):
                cached_rows = None
            cached_ids = [row["pmid"] for row in cached_rows] if cached_rows is not None else []
            if (cached_rows is not None and len(set(cached_ids)) == expected_count
                    and set(cached_ids) == set(expected_ids)):
                all_rows.extend(cached_rows)
                print(f"[{qid}] reuse {xml_path}", flush=True)
                continue

        p = {

            "db": "pubmed",
            "query_key": qkey,
            "WebEnv": webenv,
            "retstart": retstart,
            "retmax": BATCH,
            "retmode": "xml",

        }

        if email:

            p["email"] = email

        rr = req(
            session,
            "efetch.fcgi",
            p,
            api_key,
        )

        batch_rows = parse_batch(
            rr.content,
            qid,
            expected_count,
        )
        fetched_ids = [row["pmid"] for row in batch_rows]
        if len(set(fetched_ids)) != expected_count or set(fetched_ids) != set(expected_ids):
            raise RuntimeError(
                f"EFetch IDs do not match ESearch for {qid} {start}..{end} "
                f"at retstart={retstart}"
            )
        temporary_path = xml_path.with_suffix(".xml.tmp")
        temporary_path.write_bytes(rr.content)
        temporary_path.replace(xml_path)
        all_rows.extend(batch_rows)

    ids_path = segdir / "esearch_ids.txt"
    ids_tmp = ids_path.with_suffix(".txt.tmp")
    ids_tmp.write_text("\n".join(segment_ids) + "\n", encoding="utf-8")
    ids_tmp.replace(ids_path)

    return count


def main():

    ap = argparse.ArgumentParser()

    ap.add_argument(
        "--config",
        required=True,
    )

    ap.add_argument(
        "--outdir",
        required=True,
    )

    ap.add_argument(
        "--email",
        default=os.environ.get(
            "NCBI_EMAIL",
            "",
        ),
    )

    ap.add_argument(
        "--api-key",
        default=os.environ.get(
            "NCBI_API_KEY",
            "",
        ),
    )

    args = ap.parse_args()

    cfg = json.loads(
        Path(
            args.config
        ).read_text()
    )

    out = Path(
        args.outdir
    )

    out.mkdir(
        parents=True,
        exist_ok=True,
    )

    session = (
        requests.Session()
    )

    session.headers.update({

        "User-Agent":
        "sleep-frailty-systematic-review/1.1"

    })

    all_rows = []
    log_rows = []

    as_of = dt.datetime.strptime(
        cfg.get(
            "as_of_date",
            "2026/09/22",
        ),
        "%Y/%m/%d",
    ).date()

    start_all = dt.date(
        1900,
        1,
        1,
    )

    for qid, base_query in (
        cfg["queries"].items()
    ):

        full_term = dated_term(
            base_query,
            start_all,
            as_of,
        )

        total = int(
            esearch(
                session,
                full_term,
                args.email,
                args.api_key,
                False,
            )["count"]
        )

        print(
            f"[{qid}] "
            f"PubMed count="
            f"{total:,}",
            flush=True,
        )

        qdir = out / qid

        qdir.mkdir(
            exist_ok=True
        )

        (
            qdir / "query.txt"
        ).write_text(
            full_term + "\n",
            encoding="utf-8",
        )

        segments = list(
            split_segments(
                session,
                base_query,
                start_all,
                as_of,
                args.email,
                args.api_key,
            )
        )

        print(
            f"[{qid}] retrieving "
            f"in {len(segments)} "
            "date segment(s)",
            flush=True,
        )

        fetched = 0

        for i, (
            seg_start,
            seg_end,
            seg_count,
        ) in enumerate(
            segments,
            1,
        ):

            print(
                f"[{qid}] segment "
                f"{i}/{len(segments)} "
                f"{seg_start}.."
                f"{seg_end}: "
                f"{seg_count:,}",
                flush=True,
            )

            fetched += fetch_segment(

                session,
                qid,
                base_query,
                seg_start,
                seg_end,
                qdir,
                args.email,
                args.api_key,
                all_rows,

            )

        log_rows.append({

            "query_id": qid,
            "count": total,
            "fetched": fetched,
            "segments": len(
                segments
            ),
            "date_run":
                dt.date.today()
                .isoformat(),
            "as_of_date":
                as_of.isoformat(),
            "query":
                full_term,

        })

    by_pmid = {}

    for row in all_rows:

        key = (
            row["pmid"]
            or
            (
                row["doi"].lower()
                if row["doi"]
                else
                row["title"].lower()
            )
        )

        if key not in by_pmid:

            by_pmid[key] = (
                row.copy()
            )

            by_pmid[key][
                "query_ids"
            ] = {
                row["query_id"]
            }

        else:

            by_pmid[key][
                "query_ids"
            ].add(
                row["query_id"]
            )

    master = []

    for row in (
        by_pmid.values()
    ):

        row["query_ids"] = (
            ";".join(
                sorted(
                    row[
                        "query_ids"
                    ]
                )
            )
        )

        row.pop(
            "query_id",
            None,
        )

        master.append(
            row
        )

    master.sort(
        key=lambda x: (
            x.get("year", ""),
            x.get("pmid", ""),
        ),
        reverse=True,
    )

    fields = [

        "pmid",
        "doi",
        "year",
        "journal",
        "title",
        "abstract",
        "query_ids",

    ]

    with (
        out /
        "pubmed_master.csv"
    ).open(
        "w",
        newline="",
        encoding="utf-8",
    ) as f:

        w = csv.DictWriter(
            f,
            fieldnames=fields,
        )

        w.writeheader()
        w.writerows(master)

    log_fields = [

        "query_id",
        "count",
        "fetched",
        "segments",
        "date_run",
        "as_of_date",
        "query",

    ]

    with (
        out /
        "search_log.tsv"
    ).open(
        "w",
        newline="",
        encoding="utf-8",
    ) as f:

        w = csv.DictWriter(
            f,
            fieldnames=log_fields,
            delimiter="\t",
        )

        w.writeheader()
        w.writerows(
            log_rows
        )

    print(
        "PubMed complete: "
        f"{len(master):,} "
        "unique records across "
        f"{len(log_rows)} searches",
        flush=True,
    )


if __name__ == "__main__":
    main()
