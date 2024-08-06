# Establecer las rutas de las bibliotecas
.libPaths(c("~/R/x86_64-pc-linux-gnu-library/4.4", .libPaths()))

# Cargar las bibliotecas necesarias
library(MOSim)
library(rtracklayer)
library(GenomicRanges)

# Definir la lista de omics
omics_list <- c("RNA-seq")

# Definir las opciones para la simulación de DNase-seq
omics_options <- list(omicSim("RNA-seq", totalFeatures = 1000,
                              regulatorEffect = list('activator' = 0.68, 'repressor' = 0.3, 'NE' = 0.02)))

# Realizar la simulación
sim <- mosim(omics = omics_list, omicsOptions = omics_options)

# Obtener los resultados de la simulación
res <- omicResults(sim)

# Inspeccionar la estructura de los resultados
str(res)
summary(res)

# Mostrar las dimensiones de los resultados
lapply(res, dim)

# Definir la longitud de las secuencias (cromosomas)
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

# Crear el objeto Seqinfo
seqinfo <- Seqinfo(names(seqlengths), seqlengths)

# Procesar y guardar los resultados en formato BigWig
for (omic in omics_list) {
  if (omic %in% names(res)) {
    file_name <- paste0("test_2_", omic, ".bw")
    
    # Inspeccionar la estructura de los datos específicos del omic
    cat("Procesando", omic, "\n")
    print(head(res[[omic]]))
    
    # Extraer y transformar los datos
    df <- res[[omic]]
    
    # Verificar que df es un data frame
    if (!is.data.frame(df)) {
        stop("custom_chipseq no es un dataframe.")
    }

    pos_names <- rownames(df)
    counts <- df$Counts

    # Verificar las longitudes de pos_names y counts
    cat("Longitud de pos_names:", length(pos_names), "\n")
    cat("Longitud de counts:", length(counts), "\n")
    
    # Imprimir ejemplos de pos_names y counts
    cat("Ejemplos de pos_names:", head(pos_names), "\n")
    cat("Ejemplos de counts:", head(counts), "\n")

    # Verificar que pos_names y counts tienen la misma longitud
    if (length(pos_names) != length(counts)) {
      stop("La longitud de pos_names y counts no coincide.")
    }

    split_data <- strsplit(pos_names, '_', fixed = TRUE)

    # Crear un dataframe a partir de los datos divididos
    custom_chipseq_split <- do.call(rbind, split_data)
    custom_chipseq_split <- data.frame(custom_chipseq_split)

    # Renombrar las columnas correctamente
    colnames(custom_chipseq_split) <- c("seqnames", "start", "end")

    # Convertir start y end a numéricos
    custom_chipseq_split$start <- as.numeric(as.character(custom_chipseq_split$start))
    custom_chipseq_split$end <- as.numeric(as.character(custom_chipseq_split$end))

    # Verificar la longitud de custom_chipseq_split antes de añadir counts
    if (nrow(custom_chipseq_split) != length(counts)) {
      stop("La longitud de custom_chipseq_split y counts no coincide.")
    }

    # Añadir la columna de counts y reemplazar NA con 0
    custom_chipseq_split$score <- as.numeric(as.character(counts))
    custom_chipseq_split$score[is.na(custom_chipseq_split$score)] <- 0

    # Añadir el prefijo "chr" a la columna de secuencias
    custom_chipseq_split$seqnames <- paste0("chr", custom_chipseq_split$seqnames)

    # Crear el objeto GRanges
    gr <- GRanges(seqnames = custom_chipseq_split$seqnames,
                  ranges = IRanges(start = custom_chipseq_split$start, end = custom_chipseq_split$end),
                  score = custom_chipseq_split$score,
                  seqlengths = seqlengths)

    # Exportar a formato BigWig
    file_name <- "test_3_RNA_2.bw"
    export.bw(gr, file_name)

  } else {
    cat("Las columnas necesarias no están presentes en los datos para", omic, "\n")
  }
}
