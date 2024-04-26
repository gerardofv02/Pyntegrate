import gffutils
import os

path = os.path.dirname(os.path.abspath(__file__))
print(path)
db = gffutils.create_db(data=path+"/Homo_sapiens.GRCh38.109.gtf",dbfn=path+"/data/Homo_sapiens.GRCh38.109onlychr1.gtf.db")