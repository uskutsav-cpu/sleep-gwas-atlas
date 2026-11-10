import runpy,zipfile,pathlib
m=runpy.run_path('work/gstride_remotezip_inventory.py');z=zipfile.ZipFile(m['Range']())
for n,out in [('GSTRIDE_database/Database_register.csv','work/gstride_register.csv'),('GSTRIDE_database/Gait_analysis_w_Python/README.txt','work/gstride_gait_readme.txt')]:
 b=z.read(n);pathlib.Path(out).write_bytes(b);print(out,len(b))
