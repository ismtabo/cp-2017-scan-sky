# cp-2017-scan-sky
Scan sky: parallel computation practice

## Enunciado

- Dada una imagen contar el número de cuerpos celestes que aparecen en ella.
- Implementar una solución paralela con OpenMP basada en el código secuencial proporcio-
nado por los profesores.
- Respetar las estructuras de datos y el algoritmo presente en el código secuencial.

## Descripción del problema

El problema planteado en la práctica trata sobre el **etiquetado de componentes conectados** (o
**descubrimiento de regiones**). Este problema supone el etiquetado de nodos o componentes conectados
dentro de un grafo dada una heurística.
El etiquetado de componentes conectados es frecuente en campos como la visión artificial,
donde se utiliza para detectar la aparición de regiones en imágenes, así como en el reconocimiento
de imágenes o en la extracción de blobs, para su posterior conteo, filtrado o seguimiento.

## Descripción del algoritmo
Este algoritmo se compone de las siguientes fases:

1. Leer matriz de códigos de colores de la imagen
2. Etiquetar bloques de colores no vacíos con el índice del pixel
3. Copiar etiquetas a matriz auxiliar
4. Por cada celda, buscar el mínimo valor entre la etiqueta actual y las vecinas superior, inferior,
izquierda y derecha, y asignarlo como etiqueta de la celda.
5. Si se han producido cambios volver al paso 3, en caso contrario, pasar al siguiente paso.
6. Contar número de celdas cuyo valor es igual al número de posición del “pixel”.(Identificador
de región/blob)
