import pandas as pd

df_ip = pd.read_csv("SRR1204544_sort_nondup.bw.bed.txt", delimiter="\t")
df_input = pd.read_csv("SRR1204546_sort_nondup.bw.bed.txt", delimiter="\t")

df_ip["Annotation"] = df_ip["Annotation"].str.replace(r'\(.*', '', regex=True)
df_input["Annotation"] = df_input["Annotation"].str.replace(r'\(.*', '', regex=True)

print("Cantidad de vecees de elementos IP: \n", df_ip["Annotation"].value_counts())
print("\nCantidad de vecees de elementos INPUT:\n ", df_input["Annotation"].value_counts())