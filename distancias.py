# Para calcular la distancia genética.

#import array
import numpy as np
#
# Jukes-Cantor (JK69)
#
def jukes_cantor(str1, str2):
    print("Sec 1", str1, "Sec 2", str2)
    difs = count_diffs(str1,str2)
    longitud = len(str1)
    d = difs/longitud
    distancia = -3/4 * np.log(1 - (4/3 * d))
    print(distancia)
    return distancia

#
# Motoo Kimura (k71)
#


#
# Hasegawa, Kishino y Kano (HKY85)
#


#
# GTR
#

# 
# Funciones auxiliares
#

# Dif. entre dos cadenas del mismo tamaño. 
def count_diffs(str1, str2):
    count = 0
    for i in range(len(str1)):
        if str1[i] != str2[i]: count = count +1

    return count

jukes_cantor("-agta", "gagta")
jukes_cantor("ag-ta", "agata")
jukes_cantor("-agta", "gagt-")
jukes_cantor("ag-ta", "-gatt")