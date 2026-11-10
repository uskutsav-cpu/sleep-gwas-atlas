import math,json,pathlib
matrices={'modified_SOF':[[39338,14354,854],[43545,29467,5314],[15551,26346,11923]],'VA_FI':[[41051,11790,2646],[30107,33626,16102],[6846,19212,29308]]}
result={'source':'Qazi et al. JAHA 2024 Table 1; rows MVP-FI robust/prefrail/frail; columns comparator robust/prefrail/frail','source_url':'https://doi.org/10.1161/JAHA.123.033111','analyses':{}}
for label,m in matrices.items():
 n=sum(map(sum,m));a=m[2][2];b=sum(m[2])-a;c=sum(row[2] for row in m)-a;d=n-a-b-c;phi=(a*d-b*c)/math.sqrt((a+b)*(c+d)*(a+c)*(b+d));observed=(a+d)/n;expected=((a+b)*(a+c)+(c+d)*(b+d))/(n*n);k=(observed-expected)/(1-expected)
 result['analyses'][label]={'source_3x3':m,'n_complete':n,'both_frail':a,'MVP_FI_only_frail':b,'comparator_only_frail':c,'neither_frail':d,'MVP_FI_frail':a+b,'comparator_frail':a+c,'phi':phi,'agreement':observed,'kappa':k,'jaccard':a/(a+b+c)}
print(json.dumps(result,indent=2));pathlib.Path('outputs/mvp_frailty_agreement.json').write_text(json.dumps(result,indent=2)+'\n')
