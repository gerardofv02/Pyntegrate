import os
import gffutils

path = os.path.dirname(os.path.abspath(__file__))
print(path)

db = gffutils.create_db(data=path+"/datos/gencode.vM34.annotation.gtf",dbfn=path+"/datos/gencode.vM25.annotation.gtf.db")