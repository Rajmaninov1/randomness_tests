# Randomness Tests

### Test Monobit

Uno de los Test utilizados es el Monobit el cual es una función que recibe de parámetros una semilla, el número de pasos y el número de bits. La función crea una lista de ceros y unos aleatorios hasta igualar a la semilla ingresada, dependiendo de ciertos parámetros, se obtiene un cero y un uno.


### Test Frecuency Block

El enfoque de la prueba es la proporción de unos dentro de los bloques de M bits. El propósito de esta prueba es determinar si la frecuencia de unos en un bloque de M bits es aproximadamente M/2, como se esperaría bajo la suposición de aleatoriedad. Para el tamaño de bloque M=1, esta prueba degenera a la prueba 1, la prueba de Frecuencia (Monobit).

Primero Calculamos el número de bloques y se desecha lo que no sea dichos números,
Se usa un for para llevar el registro de la proporción de cada bloque, la cadenas binarias de los valores se cortan en bloques para así  mantener un registros de los uno, después actualizamos la ubicación de los cortes por último calculamos el valor de P


### Test Serial

Uno de Los Test utilizados es el Serial el cual representa una cantidad repetida del procedimiento basado en probar la uniformidad de la distribución de patrones con las longitudes dadas.

La prueba seria utiliza una expresión. Tiene la forma de una suma de cuadrados y anteriormente se suponía que tenía asintóticamente un χ2 (γ-variable) de distribución. 
El propósito de esta prueba es determinar si el número de ocurrencias de los patrones superpuestos de 2 mm de bits es aproximadamente el mismo que se esperaría para una secuencia aleatoria.

Se crea un length de tamaño de los datos, a continuación se crea varios array de zeros con numpy  los cuales después se recorren a continuación se calculan
finalmente imprime los resultados y si da menor a 1 si es un buen algoritmo random
