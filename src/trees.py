import os
import numpy as np
import networkx as nx
import matplotlib.pyplot as plt
from decimal import Decimal
from pylocluster import *
from newick import loads
from tabulate import tabulate

def nj(dmatrix, g, taxas, nombre):
    n = len(taxas)
    if n > 2:
        # PASO 1
        # Calculamos R_i
        r = np.zeros(n, dtype = float)
        for i in range(n):
            for j in range(n):
                r[i] += dmatrix[i][j]

        # 'Neighborliness' xp
        # estamos desperdiciando algo de espacio y tiempo xd    
        m = np.zeros((n,n), dtype=float)
        for i in range(n):
            for j in range(n):
                if i != j : m[i][j] = (n-2) * dmatrix[i][j] - r[i] - r[j]
                # Si i==j se queda en 0'ros
        # Taxas a unir:
        minimo = np.unravel_index(np.argmin(m, axis=None), m.shape)
        i = minimo[0]
        j = minimo[1]
        taxa_i = taxas.pop(i)
        if i < j:
            taxa_j = taxas.pop(j-1)
        else:
            taxa_j = taxas.pop(j)

        nuevo_vertice = taxa_i + taxa_j 
        taxas.append(nuevo_vertice)

        # Paso 2
        # Recalculamos la tabla de distancias xp
        dist_update = np.zeros((n-1,n-1), dtype = float)

        nuevas_distancias = []
        for v in range(n):
            if v == i:
                pass
            elif v == j:
                pass
            else :
                nuevas_distancias.append(0.5 * (dmatrix[i][v] + dmatrix[v][j] - dmatrix[i][j]) )

        nuevas_distancias.append(0)

        # aquí creo q me hace falta revisar índices para ver de cuál axis borro primero xc
        dmatrix1 = np.delete(dmatrix, [i,j], 1)
        dmatrix1 = np.delete(dmatrix1, [i,j], 0)
        n_nueva = len(taxas)

        for columna in range(n_nueva-1):
            for fila in range(n_nueva-1):
                dist_update[columna][fila] = dmatrix1[columna][fila]
        
        dist_update[:,-1] = nuevas_distancias
        dist_update[-1,:] = nuevas_distancias
        # Paso 3
        # Calculamos la branch length de ti a V y de tj a V 
        Vi = (.5 * dmatrix[i][j]) + (1/(2*(n - 2)) * (r[i] - r[j]))
        Vj = (.5 * dmatrix[i][j]) + (1/(2*(n - 2)) * (r[j] - r[i]))
        # Los quitamos de la estrella
        if(g.has_edge(0, taxa_i)) : g.remove_edge(0, taxa_i)
        if(g.has_edge(0, taxa_j)) : g.remove_edge(0, taxa_j)
        # Los agregamos a V con todo y sus branch length 
        g.add_edge(nuevo_vertice, taxa_i, weight=(round(Vi, 4)))
        g.add_edge(nuevo_vertice, taxa_j, weight=(round(Vj, 4)))

        # Seguimos calculandoooooo....
        nj(dist_update, g, taxas, nombre)
    # Paso 4
    elif n == 2:
        # Unimos los dos que quedan, by a branch length d(ti,tj)
        g.add_edge(taxas[0], taxas[1], weight=round(dmatrix[0,1], 4))
        g.remove_node(0)

        plt.figure()    
        pos = nx.planar_layout(g)
        weight_labels = nx.get_edge_attributes(g,'weight')
        nx.draw(g,pos,font_color = 'black', node_shape = 's', with_labels = True,)
        nx.draw_networkx_edge_labels(g,pos,edge_labels=weight_labels)
        plt.savefig(nombre + '.png', dpi=300, bbox_inches='tight')
        nx.write_latex(g, nombre + ".tex", caption="{nombre}", edge_label=weight_labels)
    else:
        raise Exception("Algo salió mal. :-) <3")

def neighbor_joining(matrix, headers, nombre):
    G = nx.Graph()
    G.add_node(0)
    for taxa in range(len(headers)):
        G.add_node(headers[taxa])
        G.add_edge(headers[taxa], 0)

    return nj(matrix, G, headers, nombre)

def upgma(matrix, headers, nombre):
    nwk = linkage(matrix, headers, method='upgma')
    with open(nombre + '.txt', 'w') as file:
        file.write(loads(nwk)[0].ascii_art())
        file.write("\n")
        file.write(nwk)
jk = [[0.,    0.114, 0.13,  0.068, 0.034, 0.114, 0.008, 0.172, 0.111, 0.048, 0.306,],
[0.114, 0.,    0.131, 0.117, 0.113, 0.12,  0.113, 0.175, 0.065, 0.11,  0.312],
[0.13,  0.131, 0.,    0.134, 0.13,  0.135, 0.129, 0.184, 0.127, 0.125, 0.312],
[0.068, 0.117, 0.134, 0.,    0.066, 0.12,  0.067, 0.175, 0.116, 0.062, 0.311],
[0.034, 0.113, 0.13,  0.066, 0.,    0.116, 0.033, 0.17,  0.109, 0.049, 0.305],
[0.114, 0.12,  0.135, 0.12,  0.116, 0.,    0.114, 0.171, 0.117, 0.108, 0.309],
[0.008, 0.113, 0.129, 0.067, 0.033, 0.114, 0.,    0.17,  0.11,  0.047, 0.304],
[0.172, 0.175, 0.184, 0.175, 0.17,  0.171, 0.17,  0.,    0.176, 0.169, 0.314],
[0.111, 0.065, 0.127, 0.116, 0.109, 0.117, 0.11,  0.176, 0.,    0.108, 0.309],
[0.048, 0.11,  0.125, 0.062, 0.049, 0.108, 0.047, 0.169, 0.108, 0.,    0.302],
[0.306, 0.312, 0.312, 0.311, 0.305, 0.309, 0.304, 0.314, 0.309, 0.302, 0.   ]]

mk = [[0.,    0.098, 0.11,  0.062, 0.032, 0.097, 0.007, 0.128, 0.095, 0.045, 0.178], [0.098, 0.,    0.112, 0.1,   0.099, 0.105, 0.097, 0.132, 0.061, 0.094, 0.184],
[0.11,  0.112, 0.,    0.111, 0.11,  0.113, 0.109, 0.139, 0.108, 0.105, 0.179],
[0.062, 0.1,   0.111, 0.,    0.061, 0.103, 0.062, 0.13,  0.1,   0.057, 0.183],
[0.032, 0.099, 0.11,  0.061, 0.,    0.1,   0.031, 0.127, 0.094, 0.046, 0.177],
[0.097, 0.105, 0.113, 0.103, 0.1,   0.,    0.096, 0.129, 0.101, 0.094, 0.178],
[0.007, 0.097, 0.109, 0.062, 0.031, 0.096, 0.,    0.126, 0.094, 0.043, 0.177],
[0.128, 0.132, 0.139, 0.13,  0.127, 0.129, 0.126, 0.,    0.131, 0.127, 0.178],
[0.095, 0.061, 0.108, 0.1,   0.094, 0.101, 0.094, 0.131, 0.,    0.091, 0.181],
[0.045, 0.094, 0.105, 0.057, 0.046, 0.094, 0.043, 0.127, 0.091, 0.,    0.176],
[0.178, 0.184, 0.179, 0.183, 0.177, 0.178, 0.177, 0.178, 0.181, 0.176, 0.   ]]

naive = [[0.,0.106, 0.12,  0.065, 0.033, 0.106, 0.008, 0.154, 0.103, 0.047, 0.251],
[0.106, 0.,    0.121, 0.108, 0.105, 0.111, 0.105, 0.156, 0.063, 0.102, 0.256],
[0.12,  0.121, 0.,    0.122, 0.119, 0.124, 0.119, 0.163, 0.117, 0.115, 0.255],
[0.065, 0.108, 0.122, 0.,    0.063, 0.111, 0.064, 0.156, 0.108, 0.06,  0.255],
[0.033, 0.105, 0.119, 0.063, 0.,    0.107, 0.032, 0.152, 0.101, 0.048, 0.251],
[0.106, 0.111, 0.124, 0.111, 0.107, 0.,    0.105, 0.153, 0.109, 0.101, 0.253],
[0.008, 0.105, 0.119, 0.064, 0.032, 0.105, 0.,    0.152, 0.102, 0.046, 0.25 ],
[0.154, 0.156, 0.163, 0.156, 0.152, 0.153, 0.152, 0.,    0.157, 0.151, 0.256],
[0.103, 0.063, 0.117, 0.108, 0.101, 0.109, 0.102, 0.157, 0.,    0.101, 0.253],
[0.047, 0.102, 0.115, 0.06,  0.048, 0.101, 0.046, 0.151, 0.101, 0.,    0.248],
[0.251, 0.256, 0.255, 0.255, 0.251, 0.253, 0.25,  0.256, 0.253, 0.248, 0.   ]]

tamura = [[0., 0.099, 0.111, 0.062, 0.032, 0.098, 0.007, 0.129, 0.096, 0.045, 0.179],
[0.099, 0.,    0.113, 0.1,   0.099, 0.105, 0.097, 0.133, 0.062, 0.094, 0.185],
[0.111, 0.113, 0.,    0.112, 0.111, 0.114, 0.109, 0.14,  0.108, 0.106, 0.18 ],
[0.062, 0.1,   0.112, 0.,    0.061, 0.104, 0.062, 0.131, 0.1,   0.057, 0.184],
[0.032, 0.099, 0.111, 0.061, 0.,    0.101, 0.031, 0.127, 0.094, 0.046, 0.178],
[0.098, 0.105, 0.114, 0.104, 0.101, 0.,    0.097, 0.13,  0.102, 0.094, 0.179],
[0.007, 0.097, 0.109, 0.062, 0.031, 0.097, 0.,    0.127, 0.095, 0.043, 0.178],
[0.129, 0.133, 0.14,  0.131, 0.127, 0.13,  0.127, 0.,    0.132, 0.127, 0.179],
[0.096, 0.062, 0.108, 0.1,   0.094, 0.102, 0.095, 0.132, 0.,    0.092, 0.182],
[0.045, 0.094, 0.106, 0.057, 0.046, 0.094, 0.043, 0.127, 0.092, 0.,    0.178],
[0.179, 0.185, 0.18,  0.184, 0.178, 0.179, 0.178, 0.179, 0.182, 0.178, 0.   ]]

h = ["andersoni", "barbouri", "bishopi", "californiense", "dumerilii", "laterale", "mexicanum", "alpoideum", "texanum", "tigrinum", "lacertina"]
