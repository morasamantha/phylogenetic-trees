import os
import numpy as np
import networkx as nx
import matplotlib.pyplot as plt


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
        g.add_edge(nuevo_vertice, taxa_i, weight=Vi)
        g.add_edge(nuevo_vertice, taxa_j, weight=Vj)

        # Seguimos calculandoooooo....
        nj(dist_update, g, taxas, nombre)
    # Paso 4
    elif n == 2:
        # Unimos los dos que quedan, by a branch length d(ti,tj)
        g.add_edge(taxas[0], taxas[1], weight=dmatrix[0,1])
        g.remove_node(0)
        plt.figure()    
        pos = nx.planar_layout(g)
        weight_labels = nx.get_edge_attributes(g,'weight')
        nx.draw(g,pos,font_color = 'white', node_shape = 's', with_labels = True,)
        nx.draw_networkx_edge_labels(G,pos,edge_labels=weight_labels)
        plt.savefig(os.getcwd() + '/graphs/' + nombre + '.png', dpi=300, bbox_inches='tight')
        #return g De hecho no la quiero devolver xp Vamos a ver si puede crear un png desde acá....
    else:
        raise Exception("Algo salió mal. :-) <3")


distance_matrix = np.array([[0, 5, 9, 9, 8],
                            [5, 0, 10, 10, 9],
                            [9, 10, 0, 8, 7],
                            [9, 10, 8, 0, 3],
                            [8, 9, 7, 3, 0]])
t = ['a', 'b', 'c', 'd', 'e']
# Formando la estreeellaaaa.
G = nx.Graph()
G.add_node(0)
for taxa in range(len(t)):
    G.add_node(t[taxa])
    G.add_edge(0, t[taxa])
    
nj(distance_matrix, G, t, "prueba")
