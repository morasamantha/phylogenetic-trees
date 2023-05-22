# Para calcular la distancia genética.
import numpy as np

def jukes_cantor(str1, str2):
    difs = count_diffs(str1,str2)
    longitud = len(str1)
    d = difs/longitud
    distancia = -3/4 * np.log(1 - (4/3 * d))
    print(distancia)
    return distancia

def motoo_kimura(str1, str2):
    longitud = len(str1)
    P, Q = count(str1, str2)
    P /= longitud
    Q /= longitud
    distancia = -1/2 * np.log(1-(2*P)-Q) - 1/4 * np.log(1-(2*Q))
    print(distancia)
    return distancia

#
# Hasegawa, Kishino y Kano (HKY85)
#
def hasegawa_kishino_kano(str1, str2):
    print("Sec 1", str1, "Sec 2", str2)

#
# GTR
#
def gtr(str1, str2):
    print("Sec 1", str1, "Sec 2", str2)

# 
# Funciones auxiliares
#

# Dif. entre dos cadenas del mismo tamaño. 
def count_diffs(str1, str2):
    count = 0
    for i in range(len(str1)):
        if str1[i] != str2[i]: count += 1
    return count

def count(str1, str2):
    transitions, transversions = 0,0
    for i in range(len(str1)):
        a,b = str1[i], str2[i]
        if a != b:
            if a == '-' or b == '-':
                pass
            elif (a == 'A' and b == 'G' or a == 'G' and b == 'A'):
                transitions += 1
            elif (a == 'T' and b == 'C' or a == 'C' and b == 'T'):
                transitions += 1
            else:
                transversions += 1
    return(transitions, transversions)