from __future__ import annotations

import importlib.util
import tempfile
import unittest
import zipfile
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
SCRIPT = ROOT / "frailty_paper/scripts/crosswalk_hfrs_codes_to_finngen_catalog.py"
SPEC = importlib.util.spec_from_file_location("crosswalk_hfrs_codes", SCRIPT)
CROSSWALK = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
SPEC.loader.exec_module(CROSSWALK)


class HfrsFinnGenCrosswalkTests(unittest.TestCase):
    def test_column_reference_conversion(self):
        self.assertEqual(CROSSWALK._column_index("A1"), 0)
        self.assertEqual(CROSSWALK._column_index("Z4"), 25)
        self.assertEqual(CROSSWALK._column_index("AA1"), 26)
        self.assertEqual(CROSSWALK._column_index("XFD1"), 16383)

    def test_xml_reader_ignores_incorrect_dimension(self):
        workbook = b'''<?xml version="1.0" encoding="UTF-8"?>
        <workbook xmlns="http://schemas.openxmlformats.org/spreadsheetml/2006/main"
          xmlns:r="http://schemas.openxmlformats.org/officeDocument/2006/relationships">
          <sheets><sheet name="Sheet 1" sheetId="1" r:id="rId1"/></sheets>
        </workbook>'''
        rels = b'''<?xml version="1.0" encoding="UTF-8"?>
        <Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships">
          <Relationship Id="rId1" Target="worksheets/sheet1.xml" Type="worksheet"/>
        </Relationships>'''
        worksheet = b'''<?xml version="1.0" encoding="UTF-8"?>
        <worksheet xmlns="http://schemas.openxmlformats.org/spreadsheetml/2006/main">
          <dimension ref="A1:A1"/><sheetData>
            <row r="1"><c r="A1" t="inlineStr"><is><t>TAGS</t></is></c></row>
            <row r="2"><c r="A2" t="inlineStr"><is><t>#F5</t></is></c>
              <c r="C2" t="inlineStr"><is><t>F00</t></is></c></row>
          </sheetData>
        </worksheet>'''

        with tempfile.TemporaryDirectory() as temp_dir:
            path = Path(temp_dir) / "malformed-dimension.xlsx"
            with zipfile.ZipFile(path, "w") as archive:
                archive.writestr("xl/workbook.xml", workbook)
                archive.writestr("xl/_rels/workbook.xml.rels", rels)
                archive.writestr("xl/worksheets/sheet1.xml", worksheet)
            rows, digest, size = CROSSWALK.read_xlsx_sheet(path, "Sheet 1")

        self.assertEqual(rows, [["TAGS"], ["#F5", "", "F00"]])
        self.assertEqual(len(digest), 64)
        self.assertGreater(size, 0)

    def test_code_token_split_preserves_pattern_as_distinct_token(self):
        tokens = {
            item.strip()
            for item in CROSSWALK.TOKEN_SPLIT.split("F00|G30[01]&I69")
            if item.strip()
        }
        self.assertEqual(tokens, {"F00", "G30[01]", "I69"})
        self.assertNotIn("G30", tokens)


if __name__ == "__main__":
    unittest.main()
