import pandas as pd
import Pyntegrate
import matplotlib.pyplot as plt

df = pd.read_csv("SRR1204546_sort_nondup.bw.bed.txt",delimiter="\t")

fig = Pyntegrate.plotutils.annotate_peaks(df)

plt.show()


# df["Annotation"] = df["Annotation"].str.replace(r'\(.*', '', regex=True)

# unique_values = df["Annotation"].unique()

# with open("unique_annotation_values.txt", "w") as output_file:
#     for value in unique_values:
#         output_file.write(f'{value}\n')

