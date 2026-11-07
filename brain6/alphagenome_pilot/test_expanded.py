#!/usr/bin/env python3
"""Audit the live expanded outputs without calling the remote service."""
import csv,gzip,hashlib,json,math,struct,unittest
from pathlib import Path
P=Path(__file__).resolve().parent;E=P/'results/avi_mirror/variant_evidence_avi.tsv';D=P/'results/jbrowse_expanded'
def sha(b):return hashlib.sha256(b).hexdigest()
class TestExpanded(unittest.TestCase):
 def test_inputs_and_completed_windows(self):
  with E.open() as f:rows=list(csv.DictReader(f,delimiter='\t'))
  self.assertEqual(len(rows),425);self.assertEqual(len({r['atlas_variant'] for r in rows}),425)
  self.assertEqual(json.loads((D/'input_binding.json').read_text())['evidence_sha256'],sha(E.read_bytes()))
  found=0
  for r in rows:
   folder=D/'variants'/sha(r['atlas_variant'].encode())[:20];p=folder/'provenance.json'
   if not p.exists():continue
   found+=1;m=json.loads(p.read_text());v=m['variant'];self.assertEqual(r['atlas_variant'],f"{v['chromosome']}:{v['position']}:{v['reference_bases']}>{v['alternate_bases']}")
   b=(folder/'track_window.f32.gz').read_bytes();self.assertEqual(m['track_window_sha256'],sha(b));self.assertEqual(len(gzip.decompress(b)),64032)
   a=struct.unpack('<16008f',gzip.decompress(b));self.assertEqual(len(m['tracks']),4)
   for i,t in enumerate(m['tracks']):
    ref=sum(a[i*4002+750:i*4002+1251]);alt=sum(a[i*4002+2751:i*4002+3252]);self.assertAlmostEqual(ref,t['ref_sum_501bp'],places=7);self.assertAlmostEqual(alt,t['alt_sum_501bp'],places=7)
    self.assertAlmostEqual(math.log2((alt+1)/(ref+1)),t['local_log2_ratio_plus_one'],places=10)
   self.assertNotIn('https://',json.dumps(m))
  self.assertGreater(found,0);print('Completed variant windows verified:',found)
if __name__=='__main__':unittest.main(verbosity=2)
