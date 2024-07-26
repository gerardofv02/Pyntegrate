omic_list <- sc_omicData(c("scRNA-seq", "scATAC-seq"))
cell_types <- list('Treg' = c(1:10),'cDC' = c(11:20),'CD4_TEM' = c(21:30),
      'Memory_B' = c(31:40))
sim <- scMOSim(omic_list, cell_types, numberReps = 2, 
               numberGroups = 2, diffGenes = list(c(0.2, 0.3)), feature_no = 8000, 
               clusters = 3, mean = c(2*10^6, 1*10^6,2*10^6, 1*10^6), 
               sd = c(5*10^5, 2*10^5, 5*10^5, 2*10^5), 
               regulatorEffect = list(c(0.1, 0.2), c(0.1, 0.2)))

print(sim)