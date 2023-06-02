import numpy as np
import networkx as nx
import matplotlib.pyplot as plt


def nj(dmatrix, g, taxas):
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

        dmatrix1 = np.delete(dmatrix, j, 0)
        dmatrix1 = np.delete(dmatrix1, j, 1)
        dmatrix1 = np.delete(dmatrix1, i-1, 1)
        dmatrix1 = np.delete(dmatrix1, i-1, 0)

        n_nueva = len(taxas)

        for columna in range(n_nueva-1):
            for fila in range(n_nueva-1):
                dist_update[columna][fila] = dmatrix1[columna][fila]
        
        print(dist_update)
        print(nuevas_distancias)
        dist_update[:,-1] = nuevas_distancias
        dist_update[-1,:] = nuevas_distancias

        # Paso 3
        # Calculamos la branch length de ti a V y de tj a V 
        Vi = (.5 * m[i][j]) + ((2*n -2) * (r[i] - r[j]))
        Vj = (.5 * m[i][j]) + ((2*n -2) * (r[j] - r[i]))
        # Los quitamos de la estrella
        if(g.has_edge(0, taxa_i)) : g.remove_edge(0, taxa_i)
        if(g.has_edge(0, taxa_j)) : g.remove_edge(0, taxa_j)
        # Los agregamos a V con todo y sus branch length 
        g.add_edge(nuevo_vertice, taxa_i, weight=Vi)
        g.add_edge(nuevo_vertice, taxa_j, weight=Vj)

        # Seguimos calculandoooooo....
        nj(dist_update, g, taxas)
    # Paso 4
    elif n == 2:
        # Unimos los dos que quedan, by a branch length d(ti,tj)
        g.add_edge(taxas[0], taxas[1], weight=dmatrix[0,1])

        pos = nx.spring_layout(g, seed=7)  # positions for all nodes - seed for reproducibility
        # nodes
        nx.draw_networkx_nodes(g, pos, node_size=700)

        # edges
        nx.draw_networkx_edges(g, pos, g.edges, width=6)
        # node labels
        nx.draw_networkx_labels(g, pos, font_size=20, font_family="sans-serif")
        # edge weight labels
        edge_labels = nx.get_edge_attributes(g, "weight")
        nx.draw_networkx_edge_labels(g, pos, edge_labels)
        plt.savefig('plotgraph.png', dpi=300, bbox_inches='tight')
        #return g De hecho no la quiero devolver xp Vamos a ver si puede crear un png desde acá....
    else:
        raise Exception("Algo salió mal. :-) <3")


distance_matrix = np.array([[0, 5, 9, 9],
                            [5, 0, 10, 10],
                            [9, 10, 0, 8],
                            [9, 10, 8, 0]])
t = ['a', 'b', 'c', 'd']
# Formando la estreeellaaaa.
G = nx.Graph()
G.add_node(0)
for taxa in range(len(t)):
    G.add_node(t[taxa])
    G.add_edge(0, t[taxa])
    
nj(distance_matrix, G, t)
