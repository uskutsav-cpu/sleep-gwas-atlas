"""Source-free regression tests; no claim of GWAS biological validation."""
import gzip,importlib.util,math,pathlib,struct,tempfile,types,unittest
from unittest import mock
import numpy as np
HERE=pathlib.Path(__file__).resolve().parent
import sys
sys.path.insert(0,str(HERE))
from acquire_regional_ld import fetch,index_ranges
from reconstruct_signed_ld import complete_bgzf,scalar_pearson
from prepare_finemap_inputs import orientation
from validate_native_coloc import weights

def bgzf(payload):
 gz=gzip.compress(payload,mtime=0)
 head=b'\x1f\x8b\x08\x04'+gz[4:10]+struct.pack('<H',6)+b'BC'+struct.pack('<H',2)
 n=len(head)+2+len(gz[10:]);return head+struct.pack('<H',n-1)+gz[10:]

class RegressionTests(unittest.TestCase):
 def test_complete_bgzf_and_truncated_member(self):
  a,b=bgzf(b'first\n'),bgzf(b'second\n')
  self.assertEqual(list(complete_bgzf(a+b)),[b'first\n',b'second\n'])
  self.assertEqual(list(complete_bgzf(a+b[:-2])),[b'first\n'])
 def test_bgzf_CRC_failure_is_not_ignored(self):
  a=bytearray(bgzf(b'first\n'));a[-8]^=1
  with self.assertRaises(gzip.BadGzipFile):list(complete_bgzf(bytes(a)))
 def test_alt_orientation_exact_and_swapped(self):
  self.assertEqual(orientation('C','T','T','C'),('EXACT',1))
  self.assertEqual(orientation('T','C','T','C'),('SWAP',-1))
 def test_alt_orientation_complement_and_swapped(self):
  self.assertEqual(orientation('G','A','T','C'),('COMPLEMENT',1))
  self.assertEqual(orientation('A','G','T','C'),('COMPLEMENT_SWAP',-1))
 def test_palindromic_alleles_remain_excluded(self):
  self.assertIsNone(orientation('G','C','C','G'))
  self.assertIsNone(orientation('A','T','T','A'))
 def test_exact_variant_alleles_required(self):
  self.assertIsNone(orientation('A','C','T','C'))
 def test_signed_dosage_flip_changes_sign(self):
  a=np.array([0,0,1,1,2,2]);b=np.array([0,1,1,1,1,2]);r=scalar_pearson(a,b)
  self.assertAlmostEqual(scalar_pearson(2-a,b),-r,places=14)
  self.assertAlmostEqual(r,float(np.corrcoef(a,b)[0,1]),places=14)
 def test_coloc_two_variant_closed_form(self):
  p1,p2,p12=1e-4,2e-4,1e-5
  direct=np.array([1,5*p1,9*p2,22*p1*p2,23*p12]);direct/=direct.sum()
  np.testing.assert_allclose(weights(np.log([2.,3.]),np.log([4.,5.]),p1,p2,p12),direct,rtol=1e-14,atol=1e-14)
 def test_whole_file_200_response_cannot_be_accepted_as_range(self):
  with tempfile.TemporaryDirectory() as tmp:
   p=pathlib.Path(tmp)/'body.bgzf'
   def run(*args,**kw):
    p.write_bytes(b'0123456789');p.with_suffix('.bgzf.headers').write_text('HTTP/2 200\nContent-Length:10\n');return types.SimpleNamespace(returncode=0,stderr='')
   with mock.patch('acquire_regional_ld.subprocess.run',side_effect=run):
    with self.assertRaises(RuntimeError):fetch('https://example.invalid/file',p,10,[0,9])
 def test_matching_206_range_accepts_exact_bytes(self):
  with tempfile.TemporaryDirectory() as tmp:
   p=pathlib.Path(tmp)/'body.bgzf'
   def run(*args,**kw):
    p.write_bytes(b'0123456789');p.with_suffix('.bgzf.headers').write_text('HTTP/2 206\nContent-Range: bytes 0-9/100\n');return types.SimpleNamespace(returncode=0,stderr='')
   with mock.patch('acquire_regional_ld.subprocess.run',side_effect=run):r=fetch('https://example.invalid/file',p,10,[0,9])
   self.assertEqual(r['bytes'],10)
if __name__=='__main__':unittest.main()
