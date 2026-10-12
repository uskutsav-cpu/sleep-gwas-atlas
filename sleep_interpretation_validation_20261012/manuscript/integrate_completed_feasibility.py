#!/usr/bin/env python3
"""Edit a separate manuscript copy with actual feasibility status only."""
import argparse
from pathlib import Path
from docx import Document

parser = argparse.ArgumentParser()
parser.add_argument('--input', type=Path, required=True)
parser.add_argument('--output', type=Path, required=True)
args = parser.parse_args()
if args.input.resolve() == args.output.resolve():
    raise ValueError('Preserve the preceding manuscript; write a separate version.')
doc = Document(args.input)
old = ('A separate computational benchmark is being prepared to evaluate defined '
       'interpretation and retrieval tasks; a preserved human feasibility pilot '
       'remains pending authorized execution. Neither preparation establishes improved performance.')
new = ('A registered ten-task Codex feasibility evaluation found no matched-evidence '
       'accuracy gain. Its separate validation addendum reports completed results '
       'and limitations; the human pilot remains pending authorized execution.')
matches = 0
for paragraph in doc.paragraphs:
    if old in paragraph.text:
        paragraph.text = paragraph.text.replace(old, new)
        matches += 1
if matches != 1:
    raise ValueError('Unexpected manuscript status paragraph; do not silently substitute.')
args.output.parent.mkdir(parents=True, exist_ok=True)
doc.save(args.output)
print(args.output)
