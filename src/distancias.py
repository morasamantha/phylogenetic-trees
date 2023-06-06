import numpy as np
from decimal import Decimal

def distancia(str1, str2):
    difs, l = 0, len(str1)
    for i in range(l):
        if str1[i]!= str2[i]:
            difs += 1
    distancia = difs / l 
    return distancia

def jukes_cantor(str1, str2):
    difs = 0
    l = len(str1)
    for i in range(l):
        if str1[i] != str2[i]: difs += 1
    d = difs/l
    distancia = -.75 * np.log(1 - (4/3 * d))
    return distancia

def motoo_kimura(str1, str2):
    longitud = len(str1)
    P, Q = 0,0
    for i in range(longitud):
        a,b = str1[i], str2[i]
        if a != b:
            if a == '-' or b == '-':
                pass
            elif (a == 'A' and b == 'G' or a == 'G' and b == 'A'):
                P += 1
            elif (a == 'T' and b == 'C' or a == 'C' and b == 'T'):
                P += 1
            else:
                Q += 1
    P /= longitud
    Q /= longitud
    distancia = -1/2 * np.log(1-(2*P)-Q) - 1/4 * np.log(1-(2*Q))
    return distancia

def tamura(str1, str2):
    l = len(str1)
    transitions, transversions = 0, 0
    for i in range(l):
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
    sec1 = cuenta_sitio(str1)   
    sec2 = cuenta_sitio(str2)
    theta1 = (sec1[0] + sec1[1]) / (l - sec1[2]) 
    theta2 = (sec2[0] + sec2[1]) / (l - sec2[2])
    c = theta1 + theta2 - 2 * theta1 * theta2
    P = transitions / l
    Q = transversions / l
    distancia = - c * np.log(1 - (P/c) - Q) - 0.5 * (1-c) * np.log(1 - 2*Q)
    return distancia

def cuenta_sitio(seq):
    g, c, p, a, t = 0,0,0,0,0
    for i in range(len(seq)):
        match seq[i]:
            case 'G':
                g += 1
            case 'C':
                c += 1
            case 'A':
                a += 1
            case 'T':
                t += 1
            case '-':
                p += 1
            case _:
                pass 
    return(g,c,p,a,t)

#def gtr(str1, str2):
#    print("Sec 1", str1, "Sec 2", str2)


def calcula_modelo(str1, str2, modelo):
    if (len(str1) != len(str2)):
        return 1
    #    raise Exception("Deben de ser del mismo tamaño.")
    d = 0
    match modelo:
        case 'naive':
            d = distancia(str1, str2)
        case 'jk':
            d = jukes_cantor(str1, str2)
        case 'mk':
            d = motoo_kimura(str1, str2)
        case 'tamura':
            d = tamura(str1, str2)
    return (round(d, 3))