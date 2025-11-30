import math #Libreria para facilidad de operaciones matematicas
import heapq #Libreria que se utiliza para la cola de prioridad

#Punto 1 - Lectura grafo
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
    
#Punto 2 - Longitud de tuberias

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



#Punto 3 - Sectorizacion
def grafo_adyacencia(nodos, aristas):
    #Se crea un diccionario grafo que contendra el vecino del nodo actual y su longitud
    grafo = {}
    #Se recorren todos los id del diccionario nodos
    for id_nodo in nodos.keys():
        #Por cada id del nodo se crea una lista vacia
        grafo[id_nodo] = []
    
    #Se recorre el diccionario de aristas
    for arista in aristas:
        #Se extraen los datos de las aristas
        nodo1 = arista['nodo_1']
        nodo2 = arista['nodo_2']
        longitud = arista['longitud']

        #Se guardan los datos en el grafo
        grafo[nodo1].append((nodo2,longitud))
        grafo[nodo2].append((nodo1,longitud))
    #Se regresa el diccionario grafo con los valores
    return grafo

#Algoritmo de dijkstra
def dijkstra(grafo, nodo_origen):
    #Se crea un diccionario para guardar las distancias minimas de cada nodo al nodo origen
    distancias = {}
    #Se recorren las llaves del grafo
    for nodo in grafo.keys():
        #Para cada nodo se le asigna una distancia infinita
        distancias[nodo] = float('inf')
    #El nodo origen tiene distancia de 0 ya que es el origen
    distancias[nodo_origen] = 0
    #Se crea una cola de prioridad para mantener primero los minimos
    #Se coloca inicialmente el nodo origen
    cola = [(0, nodo_origen)]
    #Se crea un conjunto vacio con set para poder almacenar los nodos ya procesados o visitados
    visitados = set()
    #Se recorre hasta que la cola se termine
    while cola:
        #Se extrae de la cola la distancia y el nodo actual
        distancia_actual, nodo_actual = heapq.heappop(cola)
        #Si el nodo actual esta en visitados entonces se continua con la siguiente ejecucion
        if nodo_actual in visitados:
            continue
        #Se agrega el nodo actual al conjunto de visitados
        visitados.add(nodo_actual)
        #Se exploran los vecinos del nodo actual
        for vecino, longitud in grafo[nodo_actual]:
            #Se calcula la distancia del nodo origen al nodo actual y se almacena en una nueva variable
            distancia_nueva = distancia_actual + longitud
            #Si la distancia es menor a la distancia del vecino actual
            #entonces se asigna su nueva distancia y se mete a la cola
            if distancia_nueva < distancias[vecino]:
                #Se asigna la nueva distancia  
                distancias[vecino] = distancia_nueva
                #Se mete a la cola
                heapq.heappush(cola, (distancia_nueva, vecino))
    #Se retorna el arreglo de distancias
    return distancias

#Funcion para identificar fuentes
def identificar_fuentes(nodos):
    #Se crea una lista para guardar las fuentes
    fuentes = []
    #Se recorre los nodos, items nos ayuda a separar la clave y el valor de lo que contiene la clave
    for id_nodo, info in nodos.items():
        #Se checa el valor para ver si es fuente es true
        if info['es_fuente']:
            #Si es_fuente es true entonces se agrega el nodo a la lista de fuentes
            fuentes.append(id_nodo)
    #Se retornan los id de los nodos que son fuentes
    return fuentes

#Funciona que asigna los sectores a los que pertenece cada nodo
def asignar_sectores(nodos, grafo):
    #Se buscan los nodos que son fuentes
    fuentes = identificar_fuentes(nodos)

    #Se verifica que haya fuentes en la red
    if not fuentes:
        #Se muestra mensaje de error y no se retorna nada
        print("Hay un error no hay fuentes en la red")
        return {}

    #Se calcula la distancia desde cada fuente a todos los demas nodos
    distancias_fuentes = {}
    #Se recorre las fuentes
    for fuente in fuentes:
        #Se asigna a la clave fuente los valores de la distancia minima entre los nodos
        distancias_fuentes[fuente] = dijkstra(grafo, fuente)

    #Se asigna cada nodo a la fuente mas cercana
    sectores = {}
    #Se asigna cada nodo a la fuente mas cercana
    for id_nodo in nodos.keys():
        #Se crea una variable para almacenar la fuente mas cercana
        fuente_mas_cercana = None
        #Se crea una variable infinito para la distancia minima
        distancia_minima = float('inf')

        #Se recooren todas las fuentes
        for fuente in fuentes:
            #Se consulta la distancia de la fuente actual con el id que se tiene
            distancia = distancias_fuentes[fuente][id_nodo]
            #Si la distancia es menor a la minima ya establecida
            if distancia < distancia_minima:
                #La distancia minima toma el valor de la distancia
                distancia_minima = distancia
                #Se asigna a la fuente mas cercana
                fuente_mas_cercana = fuente
        #Se asignan los sectores a la fuente mas cercana
        sectores[id_nodo] = fuente_mas_cercana
    #Se retornan los sectores
    return sectores

#Se identifican las tuberias cerradas
def identificar_tuberias_cerradas(sectores, aristas):
    #Se crea un arreglo para almacenar las tuberias cerradas
    tuberias_cerradas = []
    
    #Se recorren las aristas
    for arista in aristas:
        nodo1 = arista['nodo_1']
        nodo2 = arista['nodo_2']

        #Se checa si los nodos pertenecen a diferentes sectores, entonces se debe de cerrar la tuberia
        if sectores[nodo1] != sectores[nodo2]:
            #Se agregan al arreglo de tuberias_cerradas
            tuberias_cerradas.append(arista)
    #Se retorna el arreglo de las tuberias cerradas
    return tuberias_cerradas

#Funciona de sectorizacion que manda a llamar a las minifunciones
def sectorizacion(nodos, aristas):
    #Se crea el grafo de adyacencia, que nos permitira facilitar la navegacion por las aristas
    grafo = grafo_adyacencia(nodos,aristas)

    #Se asignan los sectores
    sectores = asignar_sectores(nodos,grafo)

    #Se identifican las tuberias a cerrar
    tuberias_cerradas = identificar_tuberias_cerradas(sectores, aristas)
    #Se retornan los sectores y las tuberias cerradas
    return sectores, tuberias_cerradas, grafo

#Punto 4 - Frescura del agua en función de la distancia del nodo a la fuente. 

#Funcion para calcular las distancias de los nodos a las fuentes
def calcular_distancias_fuentes(nodos, grafo):
    #Se identifican las fuentes
    fuentes = identificar_fuentes(nodos)
    #Se crea un diccionario para almacenar las distancias de los nodos
    distancias_fuentes = {}
    #Por cada fuente
    for fuente in fuentes:
        #Se calculan las distancia con dijkstra
        distancias_fuentes[fuente] = dijkstra(grafo, fuente)

    #Se retornan las distancias a las fuentes
    return distancias_fuentes

#Funcion para encontrar el nodo mas lejando por sector
def encontrar_nodo_mas_lejano(sectores, distancias_fuentes):
    #Se agrupan nodos por sector
    #Se crea un diccionario para poder agruparlos
    nodos_por_sector = {}
    #Se recorre por el numero de sectores, aparte se utiliza items para separar clave y valor
    for nodo, fuente in sectores.items():
        #Si la fuente no esta en el diccionario de nodos por sector entonces 
        if fuente not in nodos_por_sector:
            #Se agrega a la fuente en nodos por sector con una lista vacia
            nodos_por_sector[fuente] = []
        #Se agrega al nodo por sector
        nodos_por_sector[fuente].append(nodo)
    
    #Para cada sector se encuentra el nodo mas lejano
    resultados = {}
    #Se recorre los nodos_por sector, con items se descomponen la clave y el vlaor
    for fuente, nodos_sector in nodos_por_sector.items():
        #Se crean variable para almacenar el nodo_mas_lekano y la distancia maxima de ese nodo lejano
        nodo_mas_lejano = None
        distancia_maxima = -1
        #Se busca el nodo con mayor distancia a la fuente
        for nodo in nodos_sector:
             #Se consulta la distancia de la fuente actual con el id que se tiene
            distancia = distancias_fuentes[fuente][nodo]
            #Si la distancia es mayor a la maxima ya establecida
            if distancia > distancia_maxima:
                #La distancia minima toma el valor de la distancia
                distancia_maxima = distancia
                #Se asigna a la fuente mas cercana
                nodo_mas_lejano = nodo
        #Se asignan los resultados del nodo mas lejano
        resultados[fuente] = {
            'nodo_lejano': nodo_mas_lejano,
            'distancia': distancia_maxima
        }
    #Se retornan los resultados
    return resultados

#Funciona de frescura del agua
def frescura_agua(nodos,sectores,grafo):
    #Se buscan las distancias de las fuentes 
    distancias_fuentes = calcular_distancias_fuentes(nodos,grafo)

    #Se buscan los nodos mas lejanos de cada sector
    nodos_lejanos = encontrar_nodo_mas_lejano(sectores, distancias_fuentes)

    #Se retornan los nodos_lejanos
    return nodos_lejanos        


#Inicializacion del problema

archivos = [
    "FOS.txt",
    #"HAN.txt",
    #"NYT.txt",
    #"PES.txt"
]

for archivo in archivos:
    #Se logra hacer el primer punto
    num_total_nodos,  num_total_aristas, nodos, aristas, oficina, nuevos_nodos = leer_grafo(archivo)

    #Se logra hacer el segundo punto
    aristas_longitud = guardar_arista(nodos,aristas)

    # print(aristas_longitud)

    #Se logra hacer el tercer punto
    sectores, tuberias_cerradas, grafo = sectorizacion(nodos,aristas_longitud)

    # print(sectores)
    # print("*" * 80)
    # print(tuberias_cerradas)

    #Se logra hacer el cuarto punto
    nodos_lejanos = frescura_agua(nodos, sectores, grafo)

    print(nodos_lejanos)


