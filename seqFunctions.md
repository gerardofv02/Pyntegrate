# Signal Functions

Una vez se tienen las señales calculadas de nuestros archivos, se procede a explicar la variedad de funciones que existen dentro de la librería de Pyntegrate y como usarlas. La mayoría de estas funciones se encuentran dentro del archivo de SeqSignalAnalysis. Además, cada una de estas funciones tiene su propia información. Para poder verla se ha de poner el siguiente script de python:
```python
help(Pyntegrate.FilesFunction.Function)
```
Lo que devuelve la función para calcular las señales de los tipos, es un array de intervales, cuyo valor es un objeto con los siguientes campos:

- Values: Valores de las señales
- gene_name: Nombre del gen del intervalo
- gene_id: id del gen del intervalo
- strand: Indica la hebra de donde se ha sacado el intervalo
- chr: Chromosoma del intervalo
- start: Comienzo del intervalo
- stop: Fin del intervalo

Para comenzar, se va a explicar como hacer una limpieza de los datos de las señales, para así poder hacer análisis conjuntos entre distintos tipos y posteriormente se explicarán las funciones que sirven para generar las gráficas y hacer el análisis.

## Funciones de limpieza de datos

### CHIP-seq/RNA-seq o DNase-seq/RNA-seq

Esta función ha sido creada para que se devuelvan los genes coincidentes entre estos tipos. Al tener CHIP-seq y DNase-seq un procesado de datos similar, se podría usar la misma función para ambas y así encontrar los genes coincidentes.

Esta función es llamada **chip_or_atac_genes_not_used_with_rna** y tiene como parámetros de entrada las siguientes variables:

- data: clase DEseq2ResultsPrueba o parecida : Donde se almacena los datos de RNA-seq
- signal: Array de objetos: Donde se almacena el array de objetos mencionado anteriormente donde está indicada la información de la señal
- id: String: Indica por qué valor filtrar (id, gene_name,...)
- rna_column_name: String: Nombre de la columna de RNA para aplicar el filtrado
- delete_no_peaks: Boolean: Si es True, elimina aquellos genes de la señal donde el valor de la media de los valores de la señal sea 0

### Valores de la señal

En algunas funciones, es neceario unicamente añadir los valores de las señales de CHIP-seq y ATAC-seq. Para ello, se ha creado esta función donde únicamente devuelve un array con solo los valores de las señales, sin que sea un objeto.

Esta función es llamada **value_array_simple** y tiene el siguiente parámetro de entrada:

- array: Array de objetos: Donde se almacena el array de objetos mencionado anteriormente donde está indicada la información de la señal

###  CHIP-seq/DNase-seq

Esta función ha sido creada para que se devuelvan los genes coincidentes entre CHIP-seq y DNase-seq.

Esta función es llamada **chip_genes_not_used_with_atac** y tiene como parámetros de entrada las siguientes variables:

- atac_signal: Array de objetos: Donde se almacena el array de objetos mencionado anteriormente donde está indicada la información de la señal DNase
- chip_signal: Array de objetos: Donde se almacena el array de objetos mencionado anteriormente donde está indicada la información de la señal CHIP
- by: String: Indica por qué valor filtrar (id, gene_name,...)

### Calcular picos CHIP-seq

Cuando se tiene un archivo de IP y otro de INPUT en CHIP-seq, para hacer el cálculo de los picos, se tiene una función específica en Pyntegrate.  Esta función es llamada **calculate_peaks_with_gene_name** y tiene las siguientes variables de entrada:

- arrays_ip: Array de objetos: Donde se almacena el array de objetos mencionado anteriormente donde está indicada la información de la señal IP
- arrays_input: Array de objetos: Donde se almacena el array de objetos mencionado anteriormente donde está indicada la información de la señal INPUT
