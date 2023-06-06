import os
import numpy as np
import networkx as nx
import matplotlib.pyplot as plt
from decimal import Decimal
from pylocluster import *
from newick import loads

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

m = [[0.,    0.098, 0.11,  0.062, 0.032, 0.097, 0.007, 0.128, 0.095, 0.045, 0.178], [0.098, 0.,    0.112, 0.1,   0.099, 0.105, 0.097, 0.132, 0.061, 0.094, 0.184],
[0.11,  0.112, 0.,    0.111, 0.11,  0.113, 0.109, 0.139, 0.108, 0.105, 0.179],
[0.062, 0.1,   0.111, 0.,    0.061, 0.103, 0.062, 0.13,  0.1,   0.057, 0.183],
[0.032, 0.099, 0.11,  0.061, 0.,    0.1,   0.031, 0.127, 0.094, 0.046, 0.177],
[0.097, 0.105, 0.113, 0.103, 0.1,   0.,    0.096, 0.129, 0.101, 0.094, 0.178],
[0.007, 0.097, 0.109, 0.062, 0.031, 0.096, 0.,    0.126, 0.094, 0.043, 0.177],
[0.128, 0.132, 0.139, 0.13,  0.127, 0.129, 0.126, 0.,    0.131, 0.127, 0.178],
[0.095, 0.061, 0.108, 0.1,   0.094, 0.101, 0.094, 0.131, 0.,    0.091, 0.181],
[0.045, 0.094, 0.105, 0.057, 0.046, 0.094, 0.043, 0.127, 0.091, 0.,    0.176],
[0.178, 0.184, 0.179, 0.183, 0.177, 0.178, 0.177, 0.178, 0.181, 0.176, 0.   ]]
h = ['a_andersoni_mitgen.fasta', 'a_barbouri_mitgen.fasta', 'a_bishopi_mitgen.fasta','a_californiense_mitgen.fasta', 'a_dumerilii_mitgen.fasta', 'a_laterale_mitgen.fasta',
    'a_mexicanum_mitgen.fasta', 'a_talpoideum_mitgen.fasta','a_texanum_mitgen.fasta', 'a_tigrinum_mitgen.fasta', 's_lacertina_mitgen.fasta']
neighbor_joining(m, h, "prueba")