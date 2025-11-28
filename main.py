import math

#Funcion para leer los grafos del archivo
def leer_grafo(nombre_archivo):
    with open(nombre_archivo, 'r') as f:
        #Se lee la primera linea
        primer_linea = f.readline().strip().split()

        #Se crean variables para almacenar los datos
        num_total_nodos = int(primer_linea[0])
        num_total_aristas = int(primer_linea[1])

        #Se leen el resto de las lineas 
        lineas = f.readlines()
    #Se guardan los nodos en un diccionario
    nodos = {}
    #Se guarda las aristas en un arreglo como diccionarios
    aristas = []
    #Se guarda la oficina
    oficina = None
    #Se guardan los nuevos nodos
    nuevos_nodos = []

    #Se ubica la seccion en la que se encuentra actualmente
    seccion = None
    for linea in lineas:
        #Se eliminan espacio al inicio y al final
        linea = linea.strip()

        if linea.startswith('[NODES]'):
            seccion = 'NODES'
        elif linea.startswith('[EDGES]'):
            seccion = 'EDGES'
        elif linea.startswith('[OFFICE]'):
            seccion = 'OFFICE'
        elif linea.startswith('[NEW]'):
            seccion = 'NEW'
        #Se verifica que no este vacia y no empiece con corchete
        elif linea and not linea.startswith('['):
            #Se procesan los datos segun la seccion
            if seccion == 'NODES':
                #Se separan los datos
                datos = linea.split()
                #Se crea un id para el nodo
                id_nodo = int(datos[0])
                nodos[id_nodo] = {
                    'x' : float(datos[1]),
                    'y' : float(datos[2]),
                    'es_fuente' : bool(int(datos[3]))
                }
            elif seccion == 'EDGES':
                #Se separa los datos
                datos = linea.split()
                #Se agregan en una lista en formato de diccionario
                aristas.append({
                    'nodo_1' : int(datos[0]),
                    'nodo_2' : int(datos[1]),
                    'diametro' : float(datos[2])
                })
            elif seccion == 'OFFICE':
                #Se agrega el valor a la variable oficina
                oficina = int(linea)
            elif seccion == 'NEW':
                #Se separa los datos
                datos = linea.split()
                #Se crean variable para guardar los nodos y la arista
                x = datos[0]
                y = datos[1]
                diametro = datos[2]
                #Se agregan en una lista en formato de tupla
                nuevos_nodos.append((x,y,diametro))
    #Se retornan los valores agregados
    return num_total_nodos,  num_total_aristas, nodos, aristas, oficina, nuevos_nodos
     
#Longitud de las aristas o tuberias
def calcular_longitud(nodo1, nodo2):
    #Se calcula la distancia euclidiana por el teorema de pitagoras
    dx = (nodo2['x'] - nodo1['x'])
    dy = (nodo2['y'] - nodo1['y'])
    #Se saca la raiz
    distancia = math.sqrt(dx**2 + dy**2)
    #Se retorna la distancia
    return distancia

#Funcion para ir guardando los valores de la longitud en la lista de diccionarios
def guardar_arista(nodos, aristas):
    #Se recorren todas las aristas
    for arista in aristas:
        #Se accede al nodo 1
        nodo1 = nodos[arista['nodo_1']]
        #Se accede al nodo2
        nodo2 = nodos[arista['nodo_2']]
        #Se agrega al diccionario de aristas y se manda a llamar a la funcion de calc longitud
        arista['longitud'] = calcular_longitud(nodo1, nodo2)
    #Se retorna el arreglo de aristas con la nueva longitud agregada
    return aristas

    

#Inicializacion del problema

archivos = [
    "FOS.txt",
    "HAN.txt",
    "NYT.txt",
    "PES.txt"
]

for archivo in archivos:
    num_total_nodos,  num_total_aristas, nodos, aristas, oficina, nuevos_nodos = leer_grafo(archivo)

    aristas_longitud = guardar_arista(nodos,aristas)

    print(f"La longitud de las aristas es {aristas_longitud}")


