
.libPaths(c("~/R/x86_64-pc-linux-gnu-library/4.4", .libPaths()))

library(MOSim)
library(rtracklayer)
library(GenomicRanges)
omics_list<-c("DNase-seq")

data("sampleData")

custom_chipseq <- sampleData$SimDNaseseq$data
if (!is.data.frame(custom_chipseq)) {
    stop("custom_chipseq is not a dataframe.")
}

head(custom_chipseq)
print(str(custom_chipseq))


pos_names <- rownames(custom_chipseq)
counts <- custom_chipseq$Counts

split_data <- strsplit(pos_names, '_', fixed = TRUE)

# Crear un dataframe a partir de los datos divididos
custom_chipseq_split <- do.call(rbind, split_data)
custom_chipseq_split <- data.frame(custom_chipseq_split)

# Renombrar las columnas correctamente
colnames(custom_chipseq_split) <- c("seqnames", "start", "end")

# Convertir start y end a numéricos
custom_chipseq_split$start <- as.numeric(as.character(custom_chipseq_split$start))
custom_chipseq_split$end <- as.numeric(as.character(custom_chipseq_split$end))

# Añadir la columna de counts
custom_chipseq_split$score <- counts

# Añadir el prefijo "chr" a la columna de secuencias
custom_chipseq_split$seqnames <- paste0("chr", custom_chipseq_split$seqnames)
# Crear un vector con las longitudes de las secuencias cromosómicas
seqlengths <- c(
  chr1 = 195471971, chr2 = 182113224, chr3 = 160039680,
  chr4 = 156508116, chr5 = 151834684, chr6 = 149736546,
  chr7 = 145441459, chr8 = 129401213, chr9 = 124595110,
  chr10 = 130694993, chr11 = 122082543, chr12 = 120129022,
  chr13 = 120421639, chr14 = 124902244, chr15 = 104043685,
  chr16 = 98207768, chr17 = 94987271, chr18 = 90702639,
  chr19 = 61431566, chrX = 171031299, chrY = 91744698,
  chrM = 16150
)

gr <- GRanges(seqnames = custom_chipseq_split$seqnames,
              ranges = IRanges(start = custom_chipseq_split$start, end = custom_chipseq_split$end),
              score = custom_chipseq_split$score,
              seqlengths = seqlengths)



file_name <- "test_3_DNase_4.bw"
export.bw(gr, file_name)
