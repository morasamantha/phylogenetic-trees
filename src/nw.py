import numpy as np
from collections import deque

def needleman_wunsch(str1, str2, match, mismatch, indel):
    lstr1, lstr2 = len(str1), len(str2)
    # Paso 1. 
    # Inicializamos la matriz de (m+1) x (n+1), en ceros.
    matrix = np.zeros((lstr1+1, lstr2+1), dtype=int)
    # -1 (negativos) Borde, 0 diagonal, 1 izquierda, 2 arriba
    traceback = np.zeros((lstr1+1, lstr2+1), dtype=int)

    # Paso 2. 
    # Con slicing la columna 1 y la fila 1 con -i
    matrix[:,0] = [-x for x in range(lstr1+1)]
    matrix[0,:] = [-x for x in range(lstr2+1)]

    traceback[:,0] = [-x for x in range(lstr1+1)]
    traceback[0,:] = [-x for x in range(lstr2+1)]

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
            traceback[i][j] = 0

            if izquierda > maximo:
                maximo = izquierda
                traceback[i][j] = 1
            
            if arriba > maximo:
                maximo = arriba
                traceback[i][j] = 2

            matrix[i][j] = maximo

    # Paso 4. (Trace-back)
    c1 = ""
    c2 = ""
    a,b = deque(str1), deque(str2)
    i, j = lstr1, lstr2

    while traceback[i][j] > -1 and a:
        actual = traceback[i][j]
        match actual:
            case 0:
                c1 += a.pop()
                c2 += b.pop()
                i -= 1
                j -= 1
            case 1:
                c1 += "-"
                c2 += b.pop()
                j -= 1
            case 2:
                c1 += a.pop()
                c2 += "-"
                i -= 1
            case _:
                pass 
    c1 = c1 [::-1]
    c2 = c2 [::-1]
    return(c1,c2)