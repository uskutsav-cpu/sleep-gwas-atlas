#!/usr/bin/env python3
"""Count item identifiers in the public Pasquet et al. FI supplement.

Pass the downloaded supplementary DOCX path as the sole argument. This checks
the published table's structure only; it does not construct participant FIs.
"""
import collections
import hashlib
import sys
import xml.etree.ElementTree as ET
import zipfile
from pathlib import Path

source = Path(sys.argv[1])
data = source.read_bytes()
ns = "{http://schemas.openxmlformats.org/wordprocessingml/2006/main}"
with zipfile.ZipFile(source) as archive:
    root = ET.fromstring(archive.read("word/document.xml"))
tables = list(root.iter(ns + "tbl"))
table = tables[2]
groups = collections.defaultdict(list)
group = None
for row in table.findall(ns + "tr"):
    cells = [" ".join(" ".join(t.text or "" for t in c.iter(ns + "t")).split())
             for c in row.findall(ns + "tc")]
    if len(cells) == 1 and cells[0]:
        group = cells[0]
    elif len(cells) == 3 and cells[0] and not cells[1]:
        group = cells[0]
    elif len(cells) == 4 and cells[1] and cells[1] != "All 93 variables presented above":
        if group is None:
            raise ValueError("Item encountered before category")
        groups[group].append(cells[1])
items = [x for values in groups.values() for x in values]
assert len(items) == len(set(items)) == 93, (len(items), len(set(items)))
required = {"DEP_RSTLS_COM", "DEP_FFRT_COM", "DEP_GTGO_COM", "GS_EXAM_AVG_COM", "WLK_TIME_COM"}
assert required <= set(items)
print("SHA256:", hashlib.sha256(data).hexdigest())
for name, values in groups.items():
    print(f"{name}: {len(values)}")
print(f"Total distinct listed item labels: {len(items)}")
print(f"Maximum retained items, sleep removed: {len(items)-1}; sleep and exhaustion removed: {len(items)-3}")
