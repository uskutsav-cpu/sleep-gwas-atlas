"""Read the publicly linked original publisher supplement, bounded download."""
import hashlib,io,json,re,ssl,urllib.request,zipfile,xml.etree.ElementTree as E
from pathlib import Path
import certifi

ROOT=Path(__file__).resolve().parents[2]
OUT=ROOT/'outputs/sleep_frailty_ml_native_stable/benchmark_exposure'
ARTICLE='https://www.nature.com/articles/s41586-024-08430-9'
NAME='41586_2024_8430_MOESM1_ESM.docx'


def main():
    ctx=ssl.create_default_context(cafile=certifi.where())
    with urllib.request.urlopen(ARTICLE,context=ctx,timeout=30) as response:html=response.read(2000000)
    urls=[u for u in re.findall(r'href="([^"]+)"',html.decode()) if NAME in u and u.startswith('https://static-content.springer.com/')]
    assert len(set(urls))==1;url=urls[0]
    with urllib.request.urlopen(url,context=ctx,timeout=30) as response:
        data=response.read(1000001);final=response.url
    assert len(data)<1000000 and data.startswith(b'PK')
    with zipfile.ZipFile(io.BytesIO(data)) as archive:xml=archive.read('word/document.xml')
    root=E.fromstring(xml);ns={'w':'http://schemas.openxmlformats.org/wordprocessingml/2006/main'}
    paragraphs=[''.join(p.itertext()) for p in []]
    paragraphs=[''.join(p.itertext()).strip() for p in root.findall('.//w:p',ns)]
    # Word XML text is only inside w:t. Exclude field instruction or style text.
    paragraphs=[''.join(t.text or '' for t in p.findall('.//w:t',ns)).strip() for p in root.findall('.//w:p',ns)]
    target=OUT/NAME;target.write_bytes(data);plain=OUT/'reader_supplementary_methods.txt';plain.write_text('\n'.join(paragraphs)+'\n')
    receipt={'passed':True,'article_url':ARTICLE,'observed_publisher_attachment_url':url,'final_download_url':final,
       'article_HTML_sha256':hashlib.sha256(html).hexdigest(),'bytes':len(data),'sha256':hashlib.sha256(data).hexdigest(),
       'extracted_text_sha256':hashlib.sha256(plain.read_bytes()).hexdigest(),'paragraphs':len(paragraphs),
       'XML_media_basename_matches_Europe_PMC_fullTextXML':NAME,'reserved_study_outcomes_read':False}
    (OUT/'reader_supplementary_methods_retrieval.json').write_text(json.dumps(receipt,indent=2)+'\n')
    print(json.dumps(receipt),flush=True)


if __name__=='__main__':main()
