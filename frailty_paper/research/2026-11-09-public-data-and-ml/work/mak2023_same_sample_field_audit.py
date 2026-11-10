import pathlib,itertools,collections,json
root=pathlib.Path('/Volumes/Extreme SSD/Codex Archive/mak2023_same_cohort');stats=collections.Counter();examples=[]
with (root/'fi.fastGWA').open('rb') as f,(root/'fp.fastGWA').open('rb') as g:
 assert f.readline()==g.readline()==b'CHR\tSNP\tPOS\tA1\tA2\tN\tAF1\tBETA\tSE\tP\n'
 for a,b in itertools.zip_longest(f,g):
  stats['row_count']+=1
  if a is None or b is None:stats['missing_row']+=1;continue
  x=a.split(b'\t');y=b.split(b'\t');assert len(x)==len(y)==10
  stats['variant_allele_fields_equal' if x[:5]==y[:5] else 'variant_allele_fields_differ']+=1
  stats['N_text_equal' if x[5]==y[5] else 'N_text_differ']+=1
  stats['AF1_text_equal' if x[6]==y[6] else 'AF1_text_differ']+=1
  if len(examples)<5 and (x[5]!=y[5] or x[6]!=y[6]):examples.append({'variant':x[1].decode(),'FI_N':x[5].decode(),'FP_N':y[5].decode(),'FI_AF1':x[6].decode(),'FP_AF1':y[6].decode()})
res={'counts':dict(stats),'examples_of_N_or_AF1_disagreements':examples,'source_checksums_verified_in_previous_run':{'fi':'5cd10928bde27e48b82d38f097696473','fp':'3ec3de20038f3022eaa9163347b0e194'},'interpretation':'Row-wise sample count and frequency equality support a shared source genotype panel but do not identify individuals or prove identical phenotype samples or outcome coding; imputation and rounding can yield equality despite differing data.'}
pathlib.Path('outputs/mak2023_paired_source_sample_field_audit.json').write_text(json.dumps(res,indent=2)+'\n');print(json.dumps(res,indent=2))
