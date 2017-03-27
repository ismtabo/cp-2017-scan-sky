Paralelización del código secuencial
====


Explicación general
----

### Objetivo

Reconocer los distintos cuerpos celestes de una imagen tomada desde un satélite.

### Prodecimiento

Problema similares: '_Connected Component Labelling_', '_Blob detection_'.

Ambos problemas se presentan en el análisis de imágenes o en la visión artificial, cuyo objetivo 
es reconocer las distintas regiones bordeadas dentro de una imagen.

Se denomina '_blob_' a cada una de las regiones independientes del análisis resultante. Esta 
misma nomenclatura la utilizaremos en la práctica al querer diferenciar entre los '_blob_s' 
(bloques o regiones dentro de la imagen) y los bloques de la matriz (posteriormente se explica el 
análisis en bloques de la matriz del problema).

### Solución secuencial

La solución que se propone tiene las siguientes secciones:

- Carga desde fichero de una imagen preprocesada y representada en forma de matriz compuesta por 
bloques de colores etiquetados de `1` a `16`, `0` indica que no hay color(fondo de la imagen). Con halo exterior para evitar iteraciones colisiones sobre los bordes de la imagen.
- Inicialización de matriz de etiquetas en matriz auxiliar.
- Propagación de etiquetas sobre matriz auxilar
- Conteo etiquetas raíz (conteo de bloques o _blobs_).

#### Carga de fichero

El fichero se carga desde un fichero de prueba. La imagen representada en el fichero está 
preprocesada, teniendo una matriz de bloques de colores discretizados a 16 colores.

La inicialización del algoritmo general consiste en la carga a una matriz local ampliando la 
matriz con un borde, cuyo objetivo consiste en evitar colisiones o condicionales para revisar si en 
una iteración se encuentra uno en el borde de la imagen, produciendo violaciones de segmento
(`segment fault`).

#### Inicialización de matriz de etiquetas

Como matriz auxiliar se utiliza una matriz de etiquetas.

Entenderemos '_etiqueta_' como el valor de la casilla dentro de la matriz, este valor se calcula 
mediante: `i*nCols+j`. Esta función es equivalente al número de la casilla si la matriz fuese 
representada como un único vector.

Funciones principales
----

### `computation`

Función encargada de la copia del identificador blob al que perteneces (esquina superior izquierda).

Esta revisará las casillas vecinas en el siguiente orden: `superior, inferior, izquierda, derecha`.
Comparando el color de ambas casillas en la matriz original, si estas son iguales se recoge el menor 
de los índices que se encuentran en la matriz copia.

Esto quiere decir que el índice que se va propagando en la matriz será el de menor índice, 
en ese caso el que esté en la esquina superior izquierda del '_blob_'.



