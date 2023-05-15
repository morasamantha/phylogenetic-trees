from collections import deque
# para correrlo vayan al directorio donde tienen el archivo y con python3 nw.py <3
def needleman_wunsch(str1, str2, match, mismatch, indel):
    lstr1, lstr2 = len(str1), len(str2)
    # Paso 1. 
    # Inicializamos la matriz de (m+1) x (n+1), en ceros.
    matrix = [[0 for i in range(lstr1+1)] for j in range(lstr2+1)]

    # Paso 2. 
    # Llenamos la columna 1 y la fila 1 con -i
    for i in range(len(matrix)):
        matrix[i][0] -= i

    matrix[0][:] = [-x for x in range(lstr1+1)]

    # Paso 3.
    # (Este amigo es para simplificar el trace-back)
    queue = deque()

    # (En este for es donde está el error u_u)
    for i in range(1, lstr1+2):
        for j in range(1, lstr2):
            puntuacion  = match if str1[j-1] == str2[i-1] else mismatch
            diagonal = matrix[i-1][j-1] + puntuacion
            izquierda = matrix[i-1][j] + indel
            arriba = matrix[i][j-1] + indel

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
    while (queue):
        i = queue.pop()
        match i:
            case "d":
                c1 += a.popleft()
                c2 += b.popleft()
                # Como es diagonal popeamos lstr1 que no nos sirven
                for j in range(lstr1):
                    queue.pop()
            case "i":
                c1 += "-"
                c2 += b.popleft()
                # No popeamos nada
            case "a":
                c1 += a.popleft()
                c2 += "-"
                # Popeamos sólo lstr1-1
                for j in range(lstr1-1):
                    queue.pop()
    #Existe el caso de que alguno tenga algo aún
    while len(a):
        c1 += a.popleft()
        c2 = "-" + c2

    while len(b):
        c2 += b.popleft()
        c1 = "-" + c1
    print(c1,c2, sep="\n")

# sólo pa pruebas xc
# cadena 2 más grande
print("Prueba 1")
needleman_wunsch("agta", "gagta", 1,-1,-1)
# cadena 1 más grande
print("Prueba 2")
needleman_wunsch("agtaa", "gagta", 1,-1,-1)
# ambas del mismo tamaño
print("Prueba 3")
needleman_wunsch("agta", "gagt", 1,-1,-1)