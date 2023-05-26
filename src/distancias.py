# Para calcular la distancia genética.
import numpy as np

def jukes_cantor(str1, str2):
    difs = 0
    l = len(str1)
    for i in range(l):
        if str1[i] != str2[i]: difs += 1
    d = difs/l
    print(difs, l, d*4/3)
    distancia = -.75 * np.log(1 - (4/3 * d)) if d else 0 
    return distancia

def motoo_kimura(str1, str2):
    longitud = len(str1)
    P, Q = 0,0
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