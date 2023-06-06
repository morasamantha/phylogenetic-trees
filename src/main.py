import os
import itertools
import nw
import subprocess
import msa
import numpy as np
import distancias
from tabulate import tabulate
import trees

def main():
    # Estos son los fasta que vamos a leer. DEBEN estar en la carpeta de resources
    #archivos =  ['apiscoi.fasta', 'cucarachacoi.fasta', 'dromecoi.fasta', 'shrimpcoi.fasta']
    
    archivos = ['a_andersoni_mitgen.fasta', 'a_barbouri_mitgen.fasta', 'a_bishopi_mitgen.fasta','a_californiense_mitgen.fasta', 'a_dumerilii_mitgen.fasta', 'a_laterale_mitgen.fasta',
    'a_mexicanum_mitgen.fasta', 'a_talpoideum_mitgen.fasta','a_texanum_mitgen.fasta', 'a_tigrinum_mitgen.fasta', 's_lacertina_mitgen.fasta']
    
    taxas = leer_archivos(archivos)

    # Para guardar las tablas
    modelos_distancia = ['naive', 'jk', 'mk', 'tamura']
    tablas_multiples = {} # de hecho son dicts pero weno xp

    # Para la alineación uno a uno
    tablas_pares = alinear_secuencias_pares(taxas, archivos, modelos_distancia)
    
    # Para la alineación múltiple
    #taxas_multiples = alinear_secuencias_multiples(taxas)

    #for modelo in modelos_distancia:
    #   ruta_m = '/tablas/msa/' + modelo  + '.txt'
    #   tabla_m = tabla(taxas, modelo, archivos, ruta_m)
    #   tablas_multiples[modelo] = tabla_m

    # Árboles
    for modelo, tabla in tablas_pares.items():
        ruta = os.getcwd() + '/graphs/p/'
        trees.neighbor_joining(tabla, archivos, ruta+'/nj/'+modelo)
        trees.upgma(tabla, archivos, ruta+'/upgma/'+modelo)

    #for modelo, tabla in tablas_multiples.items():
    #   ruta = '/graphs/m/' + modelo
    #   trees.neighbor_joining(tabla, archivos, ruta+'/nj/'+modelo)
    #   trees.upgma(tabla, archivos, ruta+'/upgma/'+modelo)

def leer_archivos(archivos):
    taxa = []
    for archivo in archivos:
        secuencia = ""
        with open(os.getcwd() + '/resources/' + archivo) as file:
            for linea in itertools.islice(file, 1, None):
                secuencia += linea.rstrip()
        taxa.append(secuencia)
    return taxa
    
def alinear_secuencias_pares(taxas, nombres, modelos):
    tabla_pares = {} # Para que realice las alineaciones una única vez.
    tamanio = len(taxas)
    for modelo in modelos:
        ruta_p = '/tablas/nj/' + modelo + '.txt'
        matrix = np.zeros((tamanio, tamanio), dtype = float)
        tabla_pares[modelo] = matrix

    for i in range(tamanio):
        for j in range(i+1, tamanio):
            alineadas = nw.needleman_wunsch(taxas[i], taxas[j], 1, -1, -1)
            for modelo in modelos: 
                distancia = distancias.calcula_modelo(alineadas[0], alineadas[1], modelo)
                tabla_pares[modelo][i][j] = distancia 
                tabla_pares[modelo][j][i] = distancia

    for modelo, tabla in tabla_pares.items():
        with open(os.getcwd() + '/tablas/nw/' + modelo + '.txt', 'w') as f:
            f.write(str(tabla))

    return tabla_pares

def alinear_secuencias_multiples(taxas):
    args = ['python', 'msa.py', str(len(taxas))] + taxas
    salida = subprocess.check_output(args)
    with open(os.getcwd() + '/tablas/msa/msa.txt', 'w') as f:
        f.write(salida)
    return salida

def tabla(taxas, modelo, nombres, ruta):
    tamanio = len(taxas)
    matrix = np.zeros((tamanio, tamanio), dtype = float)

    for i in range(tamanio):
        for j in range(tamanio):
            if i != j: 
                matrix[i][j] = distancias.calcula_modelo(taxas[i], taxas[j], modelo)
    
    with open(os.getcwd() + ruta , 'w') as new_file:
        new_file.write(tabulate(taxas, headers=nombres, tablefmt="pretty"))
    
    return matrix

main()