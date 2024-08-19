.libPaths(c("~/R/x86_64-pc-linux-gnu-library/4.4", .libPaths()))

library(MOSim)
library(rtracklayer)
library(GenomicRanges)

omics_list <- c("ChIP-seq")

omics_options <- list(omicSim("ChIP-seq", totalFeatures = 1000))

sim <- mosim(omics = omics_list, omicsOptions = omics_options)

res <- omicResults(sim)

str(res)
summary(res)

lapply(res, dim)

seqlengths <- c(
  chr1 = 195471971, chr2 = 182113224, chr3 = 160039680,
  chr4 = 156508116, chr5 = 151834684, chr6 = 149736546,
  chr7 = 145441459, chr8 = 129401213, chr9 = 124595110,
  chr10 = 130694993, chr11 = 122082543, chr12 = 120129022,
  chr13 = 120421639, chr14 = 124902244, chr15 = 104043685,
  chr16 = 98207768, chr17 = 94987271, chr18 = 90702639,
  chr19 = 61431566, chrX = 171031299, chrY = 91744698,
  chrM = 16150,  
  chrMT = 16150 
)
seqinfo <- Seqinfo(names(seqlengths), seqlengths)

for (omic in omics_list) {
  if (omic %in% names(res)) {
    file_name <- paste0("test_2_", omic, ".bw")
    cat("Procesando", omic, "\n")
    df <- res[[omic]]
    if (!is.data.frame(df)) {
        stop("Los resultados no están en un formato de data frame.")
    }
    df$Score <- rowMeans(df, na.rm = TRUE)
    pos_names <- rownames(df)
    if (length(pos_names) != length(df$Score)) {
      stop("La longitud de pos_names y Score no coincide.")
    }

    split_data <- strsplit(pos_names, '_', fixed = TRUE)
    custom_chipseq_split <- do.call(rbind, split_data)
    custom_chipseq_split <- data.frame(custom_chipseq_split)
    colnames(custom_chipseq_split) <- c("seqnames", "start", "end")
    custom_chipseq_split$start <- as.numeric(as.character(custom_chipseq_split$start))
    custom_chipseq_split$end <- as.numeric(as.character(custom_chipseq_split$end))
    custom_chipseq_split$seqnames <- paste0("chr", custom_chipseq_split$seqnames)
    missing_seqnames <- setdiff(custom_chipseq_split$seqnames, seqnames(seqinfo))
      
    if (length(missing_seqnames) > 0) {
      cat("Las siguientes secuencias no están en seqinfo:\n")
      print(missing_seqnames)
      stop("Asegúrate de que todos los nombres de secuencias están en seqinfo.")
    }

    custom_chipseq_split$score <- as.numeric(as.character(df$Score))
    custom_chipseq_split$score[is.na(custom_chipseq_split$score)] <- 0
    gr <- GRanges(seqnames = custom_chipseq_split$seqnames,
                  ranges = IRanges(start = custom_chipseq_split$start, end = custom_chipseq_split$end),
                  score = custom_chipseq_split$score,
                  seqinfo = seqinfo)

    gr_collapsed <- reduce(gr, with.revmap=TRUE)
    score_agg <- sapply(gr_collapsed$revmap, function(x) sum(gr$score[x]))
    gr_collapsed$score <- score_agg
    gr_collapsed$revmap <- NULL
    export.bw(gr_collapsed, file_name)

  } else {
    cat("Las columnas necesarias no están presentes en los datos para", omic, "\n")
  }
}
