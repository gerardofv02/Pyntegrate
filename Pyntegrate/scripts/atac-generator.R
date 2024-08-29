.libPaths(c("~/R/x86_64-pc-linux-gnu-library/4.4", .libPaths()))
library(MOSim)

# Definir la lista de omics
omics_list <- c("RNA-seq")
omics_options <- list(omicSim("RNA-seq", totalFeatures = 10, regulatorEffect=list('activator'=1, 'repressor'=0, 'NE'=0 ), depth=10))

# Realizar la simulación
sim <- mosim(omics = omics_list, omicsOptions = omics_options)

# Obtener los resultados de la simulación
res <- tryCatch(
  {
    omicResults(sim, "RNA-seq")
  },
  error = function(e) {
    cat("Error en omicResults: ", e$message, "\n")
    NULL
  }
)

# Verificar si res es NULL debido a un error en omicResults
if (is.null(res)) {
  stop("No se pudo obtener los resultados de omicResults. Verifica los mensajes de error anteriores.")
}

# Verificar la estructura de los resultados
cat("Estructura de los resultados res:\n")
str(res)

# Preservar la columna de IDs de genes
gene_ids <- rownames(res)

# Redondear los valores de conteo
res_rounded <- as.data.frame(lapply(res, round))

# Volver a asignar las IDs de genes como nombres de fila
rownames(res_rounded) <- gene_ids

# Añadir la columna de IDs de genes al data.frame redondeado
res_rounded$GeneID <- gene_ids

# Reordenar para tener GeneID como la primera columna
res_rounded <- res_rounded[, c(ncol(res_rounded), 1:(ncol(res_rounded) - 1))]

# Verificar los primeros elementos después del redondeo
cat("Primeros elementos después del redondeo:\n")
print(head(res_rounded))

# Guardar los resultados redondeados en un archivo
file_name <- "test_2_options_RNA-seq.csv"
write.table(res_rounded, file = file_name, sep = "\t", row.names = FALSE, col.names = TRUE)
cat("Resultados guardados en:", file_name, "\n")
