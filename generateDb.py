import os
import gffutils

path = os.path.dirname(os.path.abspath(__file__))
print(path)

db = gffutils.create_db(data=path+"/data/gencode.vM25.annotation_lowdata.gtf",dbfn=path+"/data/gencode.vM25.annotation_lowdata.gtf.db")