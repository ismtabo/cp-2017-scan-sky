Ideas
====

CSR - Compress Row Storage
----

Para la reducción de iteraciones sobre los elementos de la matriz se puede hacer uso de
la característica de matriz dispera de este tipo de imágenes, que contendran muchas
casillas('pixeles') vacíos.

Una CSR o Compres Row Storage es una estructura de datos que comprime una matriz dispre-
sa según la siguiente estructura:

- Vector de filas: `rows`
- Vector de columnas: `cols`
- Vector de valores: `data`

#### Vector de filas:

Este vector se compone de una lista de índices que indica la posición inicial de la lis-
ta de índices de columnas con elementos en el `cols`.

#### Vector de columnas:

Este vector se compone de una lista de índices de las columnas con elementos según cada
fila.

#### Vector de valores:

Este vector recoge la lista de valores no nulos que contiene la matriz.

### ¿Cómo se relacionan los elementos anteriores?

Para acceder a los elementos de la matriz lo primero que debe comprobarse es el valor que contenga
la posición correspondiente a la fila dada en `rows`, al igual que el valor siguiente a su posi-
ción.

Estos valores nos darán el subarray de `cols` que contenga los índices de las columnas de la fila
que no contienen valores nulos. 

De la misma manera, también indica el subarray de `data` donde se encuentran los valores no nulos
de la fila.

#### Ejemplo

Ante la matriz:

```
10 00 00 00 -2
03 09 00 00 00
00 07 08 07 00
03 00 08 07 05
00 08 00 09 13
```

El correspondiente CSR es:

```
rows:  0  2 4 7 11 14
cols:  0  4 0 1  1  2 3 0 2 3 4 1 3  4
data: 10 -2 3 9  7  8 7 3 8 7 5 8 9 13
```

### Operaciones especiales

#### Acceso a casilla

Pseudocódigo:

```
Dado: 
M - matriz dispersa
i - indice de fila
j - indice de columna
Devolver:
d - dato correspondiente a M[i][j]

M[i][j]: d
Inicio
    rinit <- rows[i] 
    rend <- rows[i+1]
    
    para i<-rinit hasta rend hacer
        si cols[i] == n:
            di <- i
            devolver data[rinit + di]
        si cols[i] > n:
            devolver 0
    fin_para
    
    devolver 0
Fin
```

#### Iteración sobre elementos no nulos de matriz

Pseudocódigo:

```
Dado:
M - matriz dispersa

Devolver:
void

Inicio
    para i<-0 hasta len(rows) hacer
        rinit<-rows[i]
        rend <-rows[i+1]
        
        para k<-rinit hasta rend hacer
            j<-cols[rinit+k]
            d<-data[rinit+k]
            // Código...
        fin_para
    fin_para
Fin
```

