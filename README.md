# Pyntegrate

Pyntegrate es una librería para análisis de datos genómicos, la cual se ha creado con la base de otra librería llamada [metaseq](https://github.com/daler/metaseq) creada por Ryan Dale. La actualización de dicha librería viene dada por la necesidad de análisis entre más tipos de datos. Los tipos de datos aceptados por Pyntegrate son:

1.  ChIP-seq
2.  RNA-seq
3.  DNase-seq/ATAC-seq

Teniendo estos tipos de datos aceptados por la librería, también se han creado distintas funciones para poder hacer gráficos con ellos y entre ellos para poder analizarlos entre sí.

Esta librería también hace uso del software llamado [HOMER](http://homer.ucsd.edu/homer/). Es especialmente usado para identificar los enriquecimiento de los picos de CHIP-seq y para diferenciar los distintos tipos de regiones que tiene el ADN según el intervalo.


## Instrucciones de instalación de Pyntegrate

Para instalar Pyntegrate es necesario tener github instalado.
Si se está haciendo uso de Linux se deberían de hacer los siguientes comandos:
```console
sudo apt update
sudo apt install git
```

Una vez se tiene git instalado, se debería de ejecutar el siguiente comando:

```console
pip install git+https://github.com/gerardofv02/Pyntegrate.git
```

## Instalación de HOMER

Para hacer uso de las funciones usadas de HOMER (ya que no son obligatorias de usar para usar la librería), se debería de seguir la instalación del software siguiente las instrucciones de [este link](http://homer.ucsd.edu/homer/introduction/install.html)

## 