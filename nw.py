import numpy as np
from collections import deque
# para correrlo vayan al directorio donde tienen el archivo y con python3 nw.py <3
def needleman_wunsch(str1, str2, match, mismatch, indel):
    print("Prueba con las cadenas", str1, str2)
    lstr1, lstr2 = len(str1), len(str2)
    # Paso 1. 
    # Inicializamos la matriz de (m+1) x (n+1), en ceros.
    matrix = np.zeros((lstr1+1, lstr2+1), dtype=int)

    # Paso 2. 
    # Con slicing la columna 1 y la fila 1 con -i
    matrix[:,0] = [-x for x in range(lstr1+1)]
    matrix[0,:] = [-x for x in range(lstr2+1)]

    # Paso 3.
    # (Este amigo es para simplificar el trace-back)
    queue = deque()

    for i in range(1, lstr1+1):
        for j in range(1, lstr2+1):
            puntuacion  = match if str1[i-1] == str2[j-1] else mismatch
            diagonal = matrix[i-1][j-1] + puntuacion
            arriba = matrix[i-1][j] + indel
            izquierda = matrix[i][j-1] + indel

            maximo = diagonal
            queue.append("d")

            if izquierda > maximo:
                maximo = izquierda
                queue.pop()
                queue.append("i")
            
            if arriba > maximo:
                maximo = arriba
                queue.pop()
                queue.append("a")

            matrix[i][j] = maximo

    # Paso 4. (Trace-back)
    c1 = ""
    c2 = ""
    a,b = deque(str1), deque(str2)
    while len(queue) > 0:
        i = queue.pop()
        match i:
            case "d":
                c1 += a.popleft()
                c2 += b.popleft()
                # Como es diagonal popeamos lstr1 que no nos sirven
                for j in range(len(str2)):
                    if (len(queue) > 0): queue.pop()
            case "i":
                c1 = "-" + c1
                c2 += b.popleft()
                # No popeamos nada
            case "a":
                c1 += a.popleft()
                c2 = "-" + c2
                # Popeamos sólo lstr1-1
                for j in range(len(str2)-1):
                    if(len (queue) > 0): queue.pop()
        print(c1,c2, sep="\n")


    #Existe el caso de que alguno tenga algo aún
    while len(a) > 0:
        c1 += a.popleft()
        c2 = "-" + c2

    while len(b) > 0:
        c2 += b.popleft()
        c1 = "-" + c1
    #Este print lo borraremos después 
    print(c1,c2, sep="\n")
    return(c1,c2)

# sólo pa pruebas xc
# cadena 2 más grande
needleman_wunsch("agta", "gagta", 1,-1,-1)
# cadena 1 más grande
needleman_wunsch("agta", "agtaa", 1,-1,-1)
# ambas del mismo tamaño
needleman_wunsch("agta", "gagt", 1,-1,-1)