#!/usr/bin/env python3
"""Append verified reference metadata and generate journal abstract variant."""
import json,re,csv
from pathlib import Path
from audit_evidence import AREA
M=AREA/'MANUSCRIPT';p=M/'brain6_manuscript.tex';s=p.read_text();s=s.split('\\begin{thebibliography}{99}')[0]+'\\begin{thebibliography}{99}\n'
refs=json.loads((AREA/'review/verified_reference_metadata.json').read_text())
keys={'10.1038/ng.3406':'ldsc','10.1038/s41588-022-01017-y':'lava','10.1371/journal.pgen.1009218':'placo','10.1186/s13059-021-02478-w':'supergnova','10.1371/journal.pgen.1009440':'colocsusie','10.1371/journal.pgen.1004383':'colocabf','10.1371/journal.pmed.0040296':'strobe','10.1007/s10654-008-9302-y':'strega','10.1093/sleep/zsae209':'jia','10.1093/sleep/zsaf317':'xue','10.1038/s41398-026-04166-4':'zu','10.1007/s00335-026-10232-5':'lin','10.1101/2025.10.18.25338281':'schipper'}
def esc(x):return str(x).replace('&','\\&').replace('%','\\%').replace('_','\\_').replace('#','\\#')
found=set();records=[]
for r in refs:
 doi=r.get('doi','');key=keys.get(doi)
 if not key:continue
 title=r['title'];authors=r.get('authorString') or 'Authors listed in linked primary record';year=(r.get('journalInfo') or {}).get('yearOfPublication') or r.get('pubYear') or '2025';journal=((r.get('journalInfo') or {}).get('journal') or {}).get('title') or ('medRxiv preprint v1 (not peer reviewed)' if key=='schipper' else 'See primary record')
 s+=f"\\bibitem{{{key}}} {esc(authors)} {esc(title)} {esc(journal)}. {year}. \\href{{https://doi.org/{doi}}}{{doi: {esc(doi)}}}.\n"
 found.add(key);records.append(dict(key=key,doi=doi,title=title,metadata_id=r['id'],metadata_source=r['source'],verification='PRIMARY_INDEXED_METADATA;not_all_fulltexts_reviewed'))
assert found==set(keys.values()),set(keys.values())-found
s+='\\end{thebibliography}\n'
s+=r'''
\section*{Figure legends}
\textbf{Figure 1. Global genetic-correlation display.} All 72 inherited estimates across six disorders and 12 sleep traits. Stars indicate significance from the original 396-test FDR family; the displayed subset is not a new correction family. Colour represents the recorded estimate, not a causal effect. Complete standard errors, intercepts, source definitions and q values accompany the figure.

\textbf{Figure 2. Canonical LAVA non-estimation.} Counts for all seven traits; the dashed reference is the whole-family maximum of 873, not a per-trait allowance. Long sleep alone exceeds that ceiling. The family remains failed; NOT\_RUN is not a biological null.

\textbf{Figure 3. Secondary local covariance.} Three corrected-significant blocks with normal-approximation 95\% intervals from recorded covariance variance. These same-source estimates are not independent replication. The chr11 ADHD derived correlation exceeds 1; covariance values are displayed without clipping. Selection of significant blocks is descriptive; the complete 8,465-slot numerical table retains every status.

\textbf{Figure 4. Single-signal colocalization prior sensitivity.} All six estimated insomnia--ADHD intervals at default $p_{12}=10^{-5}$ and lower $10^{-6}$. The 0.80 line marks descriptive ABF-model support; it is not a replication or causal-validation threshold. Nineteen other pair candidates were unestimated and are retained in source tables. Chr5 is already a published shared region.

\section*{Supplementary materials}
The companion supplement supplies historical source methods, all numerical source tables, frozen provenance, claim/source ledgers, reporting checklist and new residual diagnostics. Figures are provided separately as vector PDF/SVG with numerical source files. The investigator-review PDF follows the separate-figure submission layout; no copyrighted third-party figure is reproduced.
\end{document}
'''
p.write_text(s)
with (M/'REFERENCES_VERIFIED.tsv').open('w') as f:
 w=csv.DictWriter(f,fieldnames=list(records[0]),delimiter='\t');w.writeheader();w.writerows(records)
abstract='''Genetic overlap between sleep traits and psychiatric disorders is established, but cross-trait screening does not validate a shared causal locus. We audited a six-disorder sleep atlas using frozen statistical families and explicit source, linkage-disequilibrium, replication and molecular gates. The 72-comparison global display retained 35 significant entries from the original 396-test false-discovery-rate family. Twenty-five pair-specific candidates occupied 20 geographic regions. Canonical LAVA had 13,745 tested and 3,720 unestimated cells among 17,465 evaluations, exceeding the frozen maximum of 873 and remaining unpromoted. Three secondary covariance blocks passed the complete 8,465-slot family threshold, without validating a protected candidate for the same pair. One significant block had a derived correlation above 1. No assessed LD matrix passed the multi-signal gate. Of six single-signal colocalization estimates, three had default-prior H4 support and one retained it under the lower prior; that chr5 region had already been reported. Four component-GWAS–QTL tests and an 18-region positional sensitivity established no effector or false-discovery-rate-significant enrichment. Independent two-trait locus replication was not established. Numerical checks reproduced available table arithmetic, but full native-input replay remains incomplete. The atlas supports a descriptive map and bounded secondary-validation analysis, without a new shared causal variant or mechanism. Failed eligibility and same-source agreement must remain distinct from biological nulls and independent replication.'''
assert 150<=len(abstract.split())<=250
(M/'MP_ABSTRACT.txt').write_text(abstract+'\n')
print('Verified references',len(records),'MP abstract words',len(abstract.split()))
