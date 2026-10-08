import fs from 'node:fs/promises';
import path from 'node:path';
import { Workbook, SpreadsheetFile } from '@oai/artifact-tool';
const input=JSON.parse(await fs.readFile(process.argv[2],'utf8'));
const packageDir=path.resolve(process.argv[3]);
const wb=Workbook.create();
const summary=wb.worksheets.add('Overview');summary.showGridLines=false;
function col(n){let s='';for(n++;n>0;n=Math.floor((n-1)/26))s=String.fromCharCode(65+(n-1)%26)+s;return s;}
const receipt=[];
for(const data of input){
  const sh=wb.worksheets.add(data.sheet);sh.showGridLines=false;
  const last=col(data.headers.length-1), nr=data.rows.length+1;
  sh.getRange(`A1:${last}${nr}`).values=[data.headers,...data.rows];
  sh.getRange(`A1:${last}${nr}`).format.font={name:'Arial',size:10};
  sh.getRange(`A1:${last}1`).format={fill:'#24445D',font:{name:'Arial',size:10,bold:true,color:'#FFFFFF'},rowHeight:42,wrapText:true};
  sh.getRange(`A1:${last}${nr}`).format.columnWidth=22;
  sh.getRange(`A2:${last}${nr}`).format.rowHeight=18;
  sh.freezePanes.freezeRows(1);sh.freezePanes.freezeColumns(2);
  for(let c=0;c<data.headers.length;c++){
    const h=data.headers[c],letter=col(c),values=data.rows.map(r=>r[c]);
    if(values.some(v=>typeof v==='number'))sh.getRange(`${letter}2:${letter}${nr}`).setNumberFormat(/(^p$|_p$|fdr|bonferroni|p_from|p_rounding)/i.test(h)?'0.000E+00':(/count|size|cases|controls|denominator|_n$|_rows|overlap.*snp/i.test(h)?'#,##0':'0.000000'));
    if(/name|definition|reason|qualification|status|interpretation|source|path|evidence|doi|queries|remaining/.test(h))sh.getRange(`${letter}1:${letter}${nr}`).format.columnWidth=48;
  }
  receipt.push({sheet:data.sheet,source:data.source,rows:data.rows.length,columns:data.headers.length});
}
summary.getRange('A2').values=[['Sleep GWAS evidence workbook']];
summary.getRange('A2').format.font={name:'Arial',size:14,bold:true};
summary.getRange('A4:C11').values=[['Measure','Count','Qualification'],
 ['Discovery comparisons',null,'Fixed exploratory family'],['BH FDR < 0.05',null,'Archived P; independently recalculated'],
 ['Replication candidates',null,'All outcomes retained'],['Outcome-side positives',null,'Sleep GWAS reused'],
 ['Fully independent two-trait positives',0,'No independent sleep input'],['Nominal heterogeneity flags',null,'Zero covariance assumed'],
 ['Native source reruns',0,'GWAS and LDSC references unavailable']];
const global=input.find(x=>x.sheet==='Global correlations'),rep=input.find(x=>x.sheet==='Replication family');
const qcol=col(global.headers.indexOf('extension_fdr')),rcol=col(rep.headers.indexOf('audit_replication_class'));
summary.getRange('B5').formulas=[[`=COUNTA('Global correlations'!A2:A${global.rows.length+1})`]];
summary.getRange('B6').formulas=[[`=COUNTIFS('Global correlations'!${qcol}2:${qcol}${global.rows.length+1},"<0.05")`]];
summary.getRange('B7').formulas=[[`=COUNTA('Replication family'!A2:A${rep.rows.length+1})`]];
summary.getRange('B8').formulas=[[`=COUNTIFS('Replication family'!${rcol}2:${rcol}${rep.rows.length+1},"EXTERNAL_OUTCOME_SIDE_REPLICATION")`]];
summary.getRange('B10').formulas=[[`=COUNTA('Heterogeneity'!A2:A${input.find(x=>x.sheet==='Heterogeneity').rows.length+1})`]];
summary.getRange('A4:C11').format.font={name:'Arial',size:10};
summary.getRange('A4:C4').format={fill:'#24445D',font:{name:'Arial',size:10,bold:true,color:'#FFFFFF'}};
summary.getRange('A4:A11').format.columnWidth=38;summary.getRange('B4:B11').format.columnWidth=10;summary.getRange('C4:C11').format.columnWidth=52;
summary.getRange('A4:C11').format.rowHeight=24;
summary.getRange('B5:B11').setNumberFormat('#,##0');
summary.getRange('A13:C15').values=[['Scope','Archived precision; source-free audits',''],['Novelty','See prior-art classifications and exact tables',''],['Data','Canonical TSV files preserve all source fields','']];
summary.getRange('A13:C15').format.font={name:'Arial',size:10};
summary.getRange('B13:C15').merge(true);summary.getRange('A13:C15').format.rowHeight=24;
wb.recalculate();
const inspect=await wb.inspect({kind:'table',range:'Overview!A4:C11',include:'values,formulas',tableMaxRows:8,tableMaxCols:3,maxChars:4000});
await fs.writeFile(path.join(packageDir,'logs/workbook_inspect.jsonl'),inspect.ndjson);
const errors=await wb.inspect({kind:'match',searchTerm:'#REF!|#DIV/0!|#VALUE!|#NAME\\?|#NUM!|#NULL!|#SPILL!|#CALC!',options:{useRegex:true,maxResults:100},maxChars:2000});
await fs.writeFile(path.join(packageDir,'logs/workbook_error_scan.jsonl'),errors.ndjson);
const counts=summary.getRange('B5:B11').values.flat();
if(JSON.stringify(counts)!==JSON.stringify([1200,603,217,23,0,7,0]))throw new Error('Workbook count mismatch '+JSON.stringify(counts));
const image=await wb.render({sheetName:'Overview',autoCrop:'all',scale:1.6,format:'png'});
await fs.writeFile(path.join(packageDir,'logs/workbook_preview.png'),new Uint8Array(await image.arrayBuffer()));
const exported=await SpreadsheetFile.exportXlsx(wb);await exported.save(path.join(packageDir,'tables/supplementary_data.xlsx'));
await fs.writeFile(path.join(packageDir,'logs/workbook_receipt.json'),JSON.stringify({sheets:receipt,overview_counts:counts,recalculated:true,missing_markers_preserved:true,exported:true,static_TSV_export:true},null,2));
console.log(JSON.stringify({worksheets:receipt.length+1,rows:receipt.reduce((n,x)=>n+x.rows,0),counts}));
