
#Funcion para leer los grafos del archivo
def leer_grafo(nombre_archivo):
    with open(nombre_archivo, 'r') as f:
        
        primer_linea = f.readline().strip().split()

        #Se crean variables para almacenar los datos
        num_total_nodos = int(primer_linea[0])
        num_total_aristas = int(primer_linea[1])

        #Se leen el resto de las lineas 
        lineas = f.readlines()
    #Se guardan los nodos en un diccionario
    nodos = {}
    #Se guarda en una tupla los valores de la arista
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
                #Se crean variable para guardar los nodos y la arista
                nodo1 = datos[0]
                nodo2 = datos[1]
                diametro = datos[2]
                #Se agregan en una lista en formato de tupla
                aristas.append((nodo1,nodo2,diametro))
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

    return num_total_nodos,  num_total_aristas, nodos, aristas, oficina, nuevos_nodos
        
#Inicializacion del problema

archivos = [
    "FOS.txt",
    "HAN.txt",
    "NYT.txt",
    "PES.txt"
]

for archivo in archivos:
    num_total_nodos,  num_total_aristas, nodos, aristas, oficina, nuevos_nodos = leer_grafo(archivo)

    print(num_total_nodos)
    print(num_total_aristas)
    print(nodos)
    print(aristas)
    print(oficina)
    print(nuevos_nodos)