import Pyntegrate

## Variables: self,column_name_padj,log2FoldChange_upper, log2FoldChange_lower, nlog10_upper,text_nlog10,text_log2FoldChange


df = Pyntegrate.results_table.DEseq2ResultsPrueba("./desq2_Sanchez.csv",  import_kwargs=dict(index_col=1))
print(df)



fig = df.volcano_plot(column_name_padj="padj", log2FoldChange_upper = 1.5,log2FoldChange_lower=-1.5,nlog10_upper=2,text_nlog10=40,text_log2FoldChange=6 )
##Por defecto sin nombres pero poner opciones para usarlos 

fig.show()