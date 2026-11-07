"""Run the frozen source-audit logic for the dementia-excluded release."""
from pathlib import Path

source = Path(__file__).resolve().parent.parent / 'stream.py'
code = source.read_text()
replacements = {
 '/Volumes/Brain6Work/physical-function-source-sensitivity-20261028/questionnaire.sumstats.gz': '/Volumes/Brain6Work/ukb_hfrs_public_20261028.sumstats.gz',
 '/Volumes/Brain6Work/ukb_hfrs_public_20261028.sumstats.gz': '/Volumes/Extreme SSD/brain6-ukb-hfrs-without-dementia-20261028.sumstats.gz',
 'https://ndownloader.figshare.com/files/61648825': 'https://ndownloader.figshare.com/files/61648822',
 '4eafec0dc6a5e07612d22873ae45bb79': '68d413550cbc2b3d7b4857ad36b7b4f7',
 'c66d2fe5eb78edb5613a8b6d06259200e798e2ab065f520c518a1fca1d1c7d50': '477f245c86aeea045d1226fc827a409b6d875615fdd800597c2611b96d34e78a',
 '280558723': '280547193',
 "'UKB_hfrs_gwas.txt.gz'": "'UKB_hfrswodementia_gwas.txt.gz'",
}
# Replace simultaneously to avoid substituting the newly introduced common path.
for i,(old,new) in enumerate(replacements.items()):
 assert old in code, old
 code=code.replace(old,f'@@REPLACE{i}@@')
for i,(old,new) in enumerate(replacements.items()):code=code.replace(f'@@REPLACE{i}@@',new)
code=code.replace("with gzip.open(target,'xt') as f:", "print(json.dumps({'source_bytes':size,'source_md5':md5.hexdigest(),'source_sha256':digest.hexdigest(),'gzip_complete':z.eof,'counts':dict(counts),'n_counts':dict(n_counts)},indent=2),flush=True)\nwith gzip.open(target,'xt') as f:")
exec(compile(code,str(source),'exec'))
