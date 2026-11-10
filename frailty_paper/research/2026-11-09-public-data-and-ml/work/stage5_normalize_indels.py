import csv,json,subprocess,hashlib
rows=[x for x in csv.DictReader(open('outputs/stage5_variant_input_manifest.csv')) if len(x['reference_allele'])!=len(x['alternate_allele'])];out=[];receipts=[]
for offset in range(0,len(rows),20):
 block=rows[offset:offset+20];regions=[f"{x['chr'][3:]}:{int(x['position_GRCh38'])-100}..{int(x['position_GRCh38'])+len(x['reference_allele'])+99}:1" for x in block];body=subprocess.run(['curl','--fail','--silent','--show-error','--max-time','60','-X','POST','https://rest.ensembl.org/sequence/region/human','-H','Content-Type: application/json','-H','Accept: application/json','--data-binary','@-'],input=json.dumps({'regions':regions}).encode(),capture_output=True,check=True).stdout;d={x['query']:x['seq'].upper() for x in json.loads(body)};assert set(d)==set(regions);receipts.append({'n':len(regions),'sha256':hashlib.sha256(body).hexdigest(),'regions':regions})
 for x,query in zip(block,regions):
  seq=d[query];oldpos=int(x['position_GRCh38']);start=oldpos-100;pos=oldpos;ref=x['reference_allele'];alt=x['alternate_allele'];assert seq[oldpos-start:oldpos-start+len(ref)]==ref
  def edited(p,r,a):
   i=p-start;assert seq[i:i+len(r)]==r;return seq[:i]+a+seq[i+len(r):]
  original=edited(pos,ref,alt);shifts=0
  while ref[-1]==alt[-1]:
   if min(len(ref),len(alt))==1:
    if pos<=start:break
    prev=seq[pos-start-1];ref=prev+ref;alt=prev+alt;pos-=1;shifts+=1
   ref=ref[:-1];alt=alt[:-1]
  while len(ref)>1 and len(alt)>1 and ref[0]==alt[0]:ref=ref[1:];alt=alt[1:];pos+=1
  assert edited(pos,ref,alt)==original
  status='normalized_within_window' if pos>start else 'requires_longer_reference_window'
  out.append({'source_variant_id':x['input_label'],'chr':x['chr'],'source_position_GRCh38':oldpos,'source_ref':x['reference_allele'],'source_alt':x['alternate_allele'],'normalized_position_GRCh38':pos,'normalized_ref':ref,'normalized_alt':alt,'position_shift_bp':pos-oldpos,'normalization_status':status,'reference_haplotype_equivalence':'verified within +/-100 bp fetched sequence','aliases':x['aliases']})
with open('outputs/stage5_indel_normalization.csv','w',newline='') as f:w=csv.DictWriter(f,fieldnames=out[0]);w.writeheader();w.writerows(out)
json.dump({'method':'left normalize anchored alleles by repeat-aware suffix shifting, then shared-prefix trimming; verify equivalent edited haplotype and reference sequence in +/-100 bp Ensembl GRCh38 context','source':'https://rest.ensembl.org/sequence/region/human','batches':receipts,'input_indels':len(rows),'shifted':sum(int(x['position_shift_bp'])!=0 for x in out),'requires_longer_window':sum(x['normalization_status']!='normalized_within_window' for x in out)},open('outputs/stage5_indel_normalization_receipt.json','w'),indent=2)
print('indels',len(out),'shifted',sum(int(x['position_shift_bp'])!=0 for x in out),'unresolved',sum(x['normalization_status']!='normalized_within_window' for x in out));print('FOXP2',next(x for x in out if x['source_variant_id']=='chr7_114534306_T_TA'))
