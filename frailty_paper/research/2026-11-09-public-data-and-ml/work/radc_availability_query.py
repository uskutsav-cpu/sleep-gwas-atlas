#!/usr/bin/env python3
"""Fetch aggregate RADC MAP overlap reports; no participant data or authentication."""
import csv, subprocess, tempfile, os
from html.parser import HTMLParser

URL='https://www.radc.rush.edu/reports/dataAvailabilityCriteria.htm'
class FirstRow(HTMLParser):
 def __init__(self): super().__init__(); self.body=False; self.tr=False; self.td=False; self.cells=[]; self.done=False; self.checked=False; self.value=''
 def handle_starttag(self, tag, attrs):
  a=dict(attrs)
  if tag=='tbody' and not self.done: self.body=True
  if self.body and tag=='tr' and not self.done: self.tr=True
  if self.tr and tag=='td': self.td=True; self.checked=False; self.value=''
  if self.td and tag=='span' and 'fa-check' in a.get('class',''): self.checked=True
 def handle_data(self,data):
  if self.td:self.value+=data
 def handle_endtag(self,tag):
  if tag=='td' and self.td:self.cells.append(('check' if self.checked else self.value.strip()));self.td=False
  if tag=='tr' and self.tr:self.done=True;self.tr=False
  if tag=='tbody':self.body=False
cookie=tempfile.NamedTemporaryFile(delete=False); cookie.close()
subprocess.run(['curl','-fsSL','-c',cookie.name,URL],stdout=subprocess.DEVNULL,check=True)
base={'studyGroup.map.value':'1','omicsGenotypingGroup.conjunction':'OR','omicsGenotypingGroup.gwasAffy.value':'1','omicsGenotypingGroup.gwasIllumina.value':'1','omicsGenotypingGroup.wgsQc.value':'1'}
queries={'questionnaire':{'substudyGroup.sleepQuestionnaire.value':'1'},'actigraphy':{'substudyGroup.actigraphy.value':'1'},'both':{'substudyGroup.conjunction':'AND','substudyGroup.sleepQuestionnaire.value':'1','substudyGroup.actigraphy.value':'1'}}
with open('outputs/stage6_radc_public_overlap_counts.csv','w',newline='') as out:
 w=csv.writer(out); w.writerow(['study','substudy','genotyping_any','participants'])
 for k,extra in queries.items():
  cmd=['curl','-fsSL','-b',cookie.name,'-X','POST',URL]
  for name,value in (base|extra).items(): cmd += ['--data-urlencode',name+'='+value]
  p=FirstRow(); p.feed(subprocess.check_output(cmd,text=True))
  if p.cells[:3]!=['check']*3 or not p.cells[3].isdigit(): raise RuntimeError('Unexpected report layout: '+repr(p.cells))
  n=int(p.cells[3]);w.writerow(['MAP',k,'Affymetrix OR Illumina OR WGS',n]);print(k,n)

 for platform in ('gwasAffy','gwasIllumina','wgsQc'):
  fields={'studyGroup.map.value':'1','substudyGroup.sleepQuestionnaire.value':'1',f'omicsGenotypingGroup.{platform}.value':'1'}
  cmd=['curl','-fsSL','-b',cookie.name,'-X','POST',URL]
  for name,value in fields.items(): cmd += ['--data-urlencode',name+'='+value]
  p=FirstRow();p.feed(subprocess.check_output(cmd,text=True))
  if p.cells[:3]!=['check']*3 or not p.cells[3].isdigit(): raise RuntimeError('Unexpected report layout: '+repr(p.cells))
  n=int(p.cells[3]);w.writerow(['MAP','questionnaire',platform,n]);print(platform,n)
os.unlink(cookie.name)
