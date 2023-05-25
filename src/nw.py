import numpy as np
from collections import deque

def needleman_wunsch(str1, str2, match, mismatch, indel):
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
    while queue:
        i = queue.pop()
        match i:
            case "d":
                if a: c1 += a.pop()
                if b: c2 += b.pop()
                # Como es diagonal popeamos lstr2 que no nos sirven
                for j in range(lstr2):
                    if queue: queue.pop()
            case "i":
                c1 += "-"
                if b: c2 += b.pop()
                # No popeamos nada
            case "a":
                if a: c1 += a.pop()
                c2 += "-"
                # Popeamos sólo lstr2-1
                for j in range(lstr2-1):
                    if queue: queue.pop()

    #Existe el caso de que alguno tenga algo aún
    while len(a) > 0:
        c1 += a.pop()
        c2 += "-"

    while len(b) > 0:
        c2 += b.pop()
        c1 += "-"
    #originalmente estaban al revés. quiero corregirlo pero mejor lo vemos al final xd
    c1 = c1 [::-1]
    c2 = c2 [::-1]

    return(c1,c2)