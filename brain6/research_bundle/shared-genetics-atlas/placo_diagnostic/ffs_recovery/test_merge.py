import io
from merge_chunks import combine
out=io.BytesIO();assert combine([io.BytesIO(b'header\na\nb\n'),io.BytesIO(b'header\nc\n')],out,[2,1])==[2,1];assert out.getvalue()==b'header\na\nb\nc\n'
for inputs,expected in [([b'header\na\n',b'wrong\nc\n'],[1,1]),([b'header\na\n',b'header\nc\n'],[1,2]),([b'header\na',b'header\nc\n'],[1,1])]:
 try:combine([io.BytesIO(x) for x in inputs],io.BytesIO(),expected)
 except ValueError:pass
 else:raise AssertionError('must reject mismatched header, count, or row termination')
print('PASS: one header, ordered rows, and failure cases')
