from decimal import *           # Toma posiciones decimales de string
from functools import partial   # Pasa funciones a widgets
from math import erfc, sqrt     # No confundir con cmath
from scipy import stats         # Requiere instalación desde pip
from tkinter import *           # Requiere instalación desde pip (pytk)
from tkinter.messagebox import showinfo  # Mensaje emergente
# Requiere instalación desde pip (PythonTurtle): sirve para graficar
from turtle import ScrolledCanvas, RawTurtle, TurtleScreen


class Estructura:
    """
    Tipo de dato estructurado que agrupa un conjunto de variables,
    necesario para la generación de números pseudoaleatorios a
    través del método de Congruencias Fundamental.
    """
    n: int          # Cantidad de números aleatorios a generar
    a: int = 7      # Semilla factor del último lugar
    c: int = 13     # Semilla factor de k lugares anteriores
    k: int = 920    # Cantidad de semillas generadas de Von Neumann
    m: int = 99991  # Módulo (actúa también como divisor entero)
    x: int = 1115   # Semilla inicial para los k lugares anteriores

    # Controlar a alto nivel que n sea mayor a 1
    def __init__(self, n: int, /) -> None:
        if n < 2000:  # En caso que necesite pocos números aleatorios
            self.k = n // 2
        self.n = n


def von_neumann(n: int, z: int, /) -> tuple[int, ...]:
    """
    Genera n números pseudoaleatorios con el método de
    Von Neumann, cuya semilla inicial es z.
    Debe controlarse desde afuera que 1000 <= z <= 9999
    y que n > 0.
    """
    a: list[int] = []  # Contendrá los números pseudoaleatorios
    x: int = z  # Favorece a entender que z es un parámetro por valor
    y: int      # Actúa como auxiliar
    for i in range(0, n):
        # Primer condicional necesario para mejorar método
        if (x // 100 == 0):     # Si los primeros dos dígitos son 00
            y = 100 - (i % 99)  # y = [2;100]
            y -= 1              # y = [1;99]
            y *= 100            # y = [100;9900]
            x += y              # Los primeros dos dígitos ya no son 00
        if (x % 100 == 0):      # Si los últimos dos dígitos son 00
            y = i % 99          # y = [0;98]
            y += 1              # y = [1;99]
            x += y              # Los últimos dos dígitos ya no son 00
        y = x ** 2      # Resulta más legible que x * x
        y %= 1000000    # Remueve los dos primeros dígitos (de 8)
        y //= 100       # Remueve los dos últimos dígitos (de 8)
        a.append(y)     # Se agrega los 4 dígitos del medio (de 8)
        x = y           # Nueva semilla
    return tuple(a)     # Conviene las tuplas por los índices


def congruencias_fundamental(
    cf: Estructura,
    /) -> tuple[tuple[tuple[int, ...],
                      float], ...] | None:
    """
    Genera n números pseudoaleatorios con el método de
    Congruencias Fundamental, utilizando una estructura por defecto.
    Debe controlarse desde afuera que n > 0.
    Devuelve una tupla conteniendo n elementos, cada uno posee
    una tupla de dos partes; éstas son el conjunto de dígitos
    del número aleatorio y el número en formato flotante.
    En caso de que exista un error, se devuelve None; por lo tanto,
    debe revisarse posteriormente si la estructura corresponde.
    """
    p: int = 1  # Cantidad máxima de dígitos
    v: tuple[int, ...] = von_neumann(cf.k, cf.x)
    if v[0] >= cf.m:
        print("Error: el elemento 1 de la sucesion es",
              "mayor o igual al modulo.")
        return None
    q: int = (cf.a * v[cf.k-1] + cf.c * v[0]) % cf.m
    # p se basa en la semilla más grande, no en el módulo
    while ((10 ** p) <= q):
        p += 1
    r: float = q / cf.m  # Es lo mismo que q, pero como flotante
    u: list[float] = []
    y: list[int] = []
    u.append(r)
    y.append(q)
    for i in range(1, cf.k):
        if (v[i] >= cf.m):
            print("Error: el elemento", i+1,
                  "de la sucesion es mayor o igual al modulo.")
            return None
        q = (cf.a * y[i-1] + cf.c * v[i]) % cf.m
        while ((10 ** p) <= q):
            p += 1
        r = q / cf.m
        u.append(r)
        y.append(q)
    for i in range(cf.k, cf.n):
        if (y[i-cf.k] >= cf.m):
            print("Error: el elemento ", i+1,
                  " de la sucesion es mayor o igual al modulo.")
            return None
        q = (cf.a * y[i-1] + cf.c * y[i-cf.k]) % cf.m
        while ((10 ** p) <= q):
            p += 1
        r = q / cf.m
        u.append(r)
        y.append(q)
    b: int = 0
    e: int
    f: int
    w: list[int]
    z: list[tuple[tuple[int, ...], float]] = []
    for i in range(cf.n):
        f = 10 ** p
        w = []
        while f > 10:
            e = y[i] % f
            e *= 10
            e //= f
            w.append(e)
            b += 1
            f //= 10
        e = y[i] % 10
        w.append(e)
        z.append((tuple(w), u[i]))
    return tuple(z)     # Conviene las tuplas por los índices


def monobits(x: tuple[tuple[tuple[int, ...],
                            float], ...], /) -> bool:
    """
    Prueba que los dígitos obtenidos y los números flotantes
    se distribuyen aleatoriamente (equitativamente).
    Esto se realiza dividiendo el conjunto de números
    pseudoaleatorios en dos partes y revisando si su diferencia
    en cantidades no supera el nivel de tolerancia permitido (alfa),
    basado en una distribución normal.
    Para los dígitos, se utilizan los grupos [0;4] y [5;9].
    Para los flotantes, se utilizan los grupos [0.0;0.5) y [0.5;1.0).
    """
    ALFA = Decimal('0.01')      # Es más preciso que flotante
    L = len(x)                  # Cantidad de números flotantes
    D = len(x[0][0])            # Proporción de dígitos por flotante
    u: int = 0                  # Contador de dígitos
    f: int = 0                  # Contador de flotantes
    ''' En el resto de código, se cuenta los dígitos y los flotantes.
    En caso de pertenecer a la mitad superior del dominio,
    se cuenta; sino, se realiza un descuento.
    Tras obtener las diferencias contadas, se obtiene el cuadrado de
    éstos (para asegurar que sean positivos); dividido por el doble de
    cantidad de números, distinguiéndose entre dígitos y flotantes.
    Ese resultado termina operado por una raíz cuadrada (por separado)
    Esto se debe a que originalmente se opera de la siguiente forma:
    Z = [|S| / sqrt(N)] / sqrt(2)
    donde Z es el estadístico que se someterá a la función error
    complementaria, S es la diferencia encontrada de elementos (se
    halla dentro de una función de absoluto); sqrt() es una raíz
    cuadrada y N es la cantidad de números pseudoaleatorios.
    Por último, se prueba si los estadísticos Z son mayores o iguales
    al nivel de tolerancia alfa; ajustado a 0.01.
    '''
    for i in x:
        for n in i[0]:
            u = u+1 if (n >= 5) else u-1
        f = f+1 if (i[1] >= 0.5) else f-1
    m: float = (f ** 2) / (L * 2)
    b: float = (u ** 2) / (L * D * 2)
    return (erfc(sqrt(m)) >= ALFA and erfc(sqrt(b)) >= ALFA)


def chi_cuadrado(x: tuple[tuple[tuple[int, ...],
                                float], ...], /) -> bool:
    """
    Prueba que los dígitos obtenidos y los números flotantes
    se distribuyen aleatoriamente (equitativamente).
    Esto se realiza dividiendo el conjunto de números
    pseudoaleatorios en diez partes y revisando si la diferencia
    cuadrada entre lo esperado y lo obtenido no supera el nivel de
    tolerancia permitido, basado en una distribución chi cuadrado.
    Para los dígitos, se usa cada uno de ellos del sistema decimal.
    Para los flotantes, se fragmenta en porciones de 0.1; donde
    se incluye el valor inferior y se excluye el valor superior,
    por ejemplo: [0.0;0.1).
    """
    C = 14.6837         # alfa = 0.1 con gl = 9
    L = len(x)          # Cantidad de números flotantes
    D = len(x[0][0])    # Proporción de dígitos por flotante
    x: int              # Auxiliar para fragmentador de flotantes
    v: dict[int, int] = {}  # Diccionario para dígitos decimales
    w: dict[int, int] = {}  # Diccionario para intervalos flotantes
    for i in range(10):  # Inicia cada contador de los dos dict
        v[i]: int = 0
        w[i]: int = 0
    ''' En el resto de código, se cuenta los dígitos y los flotantes.
    Tras obtener las cantidades contadas, se suma los cuadrados de
    las diferencias entre lo esperado y lo observiado; dividido por lo
    esperado, optimizando lo último al realizar tal operación al final.
    Por último, se prueba si los estadísticos no superan el máximo
    permitido por el nivel de tolerancia alfa; ajustado a 0.01, en
    conjunto con los 9 grados de libertad (10-1 porque la probabilidad
    de todo sumado es igual a 1).
    '''
    for i in x:
        for n in i[0]:
            v[n] = v[n] + 1
        x = int(Decimal(i[1])//Decimal("0.1"))
        w[x] = w[x] + 1
    ce: float = 0.0
    cf: float = 0.0
    for i in range(10):
        ce += (v[i] - (L * D / 10)) ** 2
        cf += (w[i] - (L / 10)) ** 2
    ce /= L * D
    cf /= L
    ce *= 10
    cf *= 10
    return (ce < C and cf < C)


def poker(x: tuple[tuple[tuple[int, ...],
                         float], ...], /) -> bool:
    """Prueba grupos de números juntos como una mano de póker y
    compara cada mano con la mano esperada usando Chi-cuadrado.
    La prueba se utiliza para analizar la frecuencia con la
    que se repiten los dígitos en números pseudoaleatorios
    individuales.
    Determina si los números cumplen con las propiedades de
    uniformidad e independencia.
    """
    L = len(x)          # Se utiliza como constante por repitencia
    UNI = [1, 1, 1, 1]  # Todas diferentes
    PAR = [2, 1, 1, 1]  # Un par (resto diferentes)
    DUP = [2, 2, 1, 1]  # Dos pares (el restante único)
    TER = [3, 2, 1, 1]  # Una tercia (resto diferentes)
    FUL = [3, 2, 2, 1]  # Full house (tercia y par)
    CUA = [4, 3, 2, 1]  # Póker (4 iguales)
    QUI = [5, 4, 3, 2]  # Quintilla (todas iguales)
    v: list[int]        # Cantidad de repeticiones de elemento
    y: dict[Decimal, int] = {}  # Diccionario de esperado vs observado
    cu: Decimal = Decimal("0.3024")*L   # Todas diferentes
    cp: Decimal = Decimal("0.504")*L    # Un par (resto diferentes)
    cd: Decimal = Decimal("0.108")*L    # Dos pares (el restante único)
    ct: Decimal = Decimal("0.072")*L    # Una tercia (resto diferentes)
    cf: Decimal = Decimal("0.009")*L    # Full house (tercia y par)
    cc: Decimal = Decimal("0.0045")*L   # Póker (4 iguales)
    cq: Decimal = Decimal("0.0001")*L   # Quintilla (todas iguales)
    y[cu] = 0   # Todas diferentes
    y[cp] = 0   # Un par (resto diferentes)
    y[cd] = 0   # Dos pares (el restante único)
    y[ct] = 0   # Una tercia (resto diferentes)
    y[cf] = 0   # Full house (tercia y par)
    y[cc] = 0   # Póker (4 iguales)
    y[cq] = 0   # Quintilla (todas iguales)
    for i in x:
        v = []                      # Debe limpiarse en cada iteración
        for n in range(len(i[0])-1):  # El último siempre será 1(uno)
            v.append(i[0][n:len(i[0])].count(i[0][n]))
        v = sorted(v, reverse=True)  # Necesario para comparar
        ''' En este bloque condicional, se cuenta cada mano de
        póker según el patrón que sigue. Se ha ordenado priorizando
        a los que tengan mayor probabilidad, esperando mayor
        frecuencia de éstos; mientras que los menos recurrentes
        se consultan más abajo.
        '''
        if v == PAR:        # Un par (resto diferentes)
            y[cp] = y[cp]+1
        elif v == UNI:      # Todas diferentes
            y[cu] = y[cu]+1
        elif v == DUP:      # Dos pares (el restante único)
            y[cd] = y[cd]+1
        elif v == TER:      # Una tercia (resto diferentes)
            y[ct] = y[ct]+1
        elif v == FUL:      # Full house (tercia y par)
            y[cf] = y[cf]+1
        elif v == CUA:      # Póker (4 iguales)
            y[cc] = y[cc]+1
        elif v == QUI:      # Quintilla (todas iguales)
            y[cq] = y[cq]+1
        else:  # Si no entra en ninguna categoría
            return False  # Tiene que ser un error
    w: Decimal      # Variable auxiliar para guardar temporalmente
    '''Mientras la frecuencia esperada sea menor de 5, se agrupa
    la fila inferior con la inmediata superior hasta que la suma
    sea al menos 5. Si al quedar sólo dos grupos, sigue siendo
    menor a 5; se considera que no pasa la prueba, a causa de
    no ser lo suficientemente contundente.
    '''
    if cq < 5:
        w = Decimal("0.0046")*L
        y[w] = y[cc] + y[cq]
        y.pop(cq)
        y.pop(cc)
        cc = w
        if cc < 5:
            w = Decimal("0.0136")*L
            y[w] = y[cf] + y[cc]
            y.pop(cc)
            y.pop(cf)
            cf = w
            if cf < 5:
                w = Decimal("0.0856")*L
                y[w] = y[ct] + y[cf]
                y.pop(cf)
                y.pop(ct)
                ct = w
                if ct < 5:
                    w = Decimal("0.1936")*L
                    y[w] = y[cd] + y[ct]
                    y.pop(ct)
                    y.pop(cd)
                    cd = w
                    if cd < 5:
                        w = Decimal("0.496")*L
                        y[w] = y[cu] + y[cd]
                        y.pop(cd)
                        y.pop(cu)
                        cu = w
                        if cu < 5:
                            return False  # No hay más de 10 elementos
    z: Decimal = Decimal(0)
    '''En el bucle, se suma el cuadrado de la diferencia entre
    lo observado y lo esperado; dividiéndose por lo esperado
    '''
    for w in y.items():
        z += (((Decimal(w[1])) - w[0]) ** 2) / w[0]
    '''Según la cantidad de categorías, se comprueba el estadístico;
    para alfa=0.1 y grados de libertad iguales a la cantidad de
    categorías menos uno (la suma de probabilidades de todas las
    categorías es 1).
    '''
    T = len(y)
    if T == 7:
        return (z < 10.6446)    # alfa = 0.1 con gl = 6
    elif T == 6:
        return (z < 9.2363)     # alfa = 0.1 con gl = 5
    elif T == 5:
        return (z < 7.7794)     # alfa = 0.1 con gl = 4
    elif T == 4:
        return (z < 6.2514)     # alfa = 0.1 con gl = 3
    elif T == 3:
        return (z < 4.6052)     # alfa = 0.1 con gl = 2
    elif T == 2:
        return (z < 2.7055)     # alfa = 0.1 con gl = 1
    else:                       # Para llegar hasta aquí
        return False            # Debe ser un error


def rachas(x: tuple[tuple[tuple[int, ...],
                          float], ...], /) -> bool:
    """
    Prueba si los números flotantes siguen algún patrón para hallarse
    por debajo o arriba de la media. Utiliza la distribución normal
    para realizar esta comparación, basándose en el siguiente
    programa:
    https://www.geeksforgeeks.org/runs-test-of-randomness-in-python/
    """
    ALFA = Decimal('0.01')
    mediana: float = 0.5
    observado: int = 0      # Cuenta cada racha (cambio)
    pos: int = 0            # El contador no descontará
    neg: int = 0            # Debe distinguirse de pos
    l: list[float] = []
    for i in x:                 # Por cada número flotante
        l.append(i[1])          # Guardarlo en una lista
    t: tuple[float, ...] = tuple(l)
    # El recorrido del bucle termina siendo circular
    for i in range(len(t)):
        # Acceder con índice -1 en Python es seguro (último elemento)
        if (t[i] >= mediana and t[i-1] < mediana) or \
                (t[i] < mediana and t[i-1] >= mediana):
            observado += 1
        if(t[i]) >= mediana:
            pos += 1
        else:
            neg += 1
    esperado: float = ((2*pos*neg)/(pos+neg))+1
    desvio_estandar: float = sqrt((2*pos*neg*(2*pos*neg-pos-neg)) /
                                  (((pos+neg)**2)*(pos+neg-1)))
    z: float = (observado-esperado)/desvio_estandar
    return (stats.norm.cdf(-abs(z)) >= ALFA)


def alta(e_peso: Entry, c_peso: Label, e_precio: Entry, c_precio: Label,
         lista: dict, listado: Listbox, /) -> None:
    peso: int
    precio: Decimal
    invalido: bool = False
    c_peso.config(text="")
    c_precio.config(text="")
    s_peso: str = e_peso.get()
    try:
        peso = int(s_peso)
        if peso <= 0:
            c_peso.config(text="No es un número natural")
            invalido = True
    except:
        c_peso.config(text="No es un número natural")
        invalido = True
    s_precio: str = e_precio.get()
    try:
        precio = Decimal(s_precio)
        if precio <= 0:
            c_precio.config(text="Valor inválido como precio")
            invalido = True
    except:
        c_precio.config(text="No es un número")
        invalido = True
    if invalido:
        return
    pos: int = len(lista)
    for indice, clave in enumerate(sorted(lista.keys())):
        if peso < clave:
            pos = indice
            break
        elif peso == clave:
            showinfo("Alta", "Un camión con tal capacidad "
                     + "de peso ya existe.")
    e_peso.delete(0, len(s_peso))
    e_precio.delete(0, len(s_precio))
    precio = precio.quantize(Decimal("1.00"))
    lista[peso] = precio
    listado.insert(pos, str(peso) + "kg: $" + str(precio))


def baja(lista: dict, listado: Listbox, /) -> None:
    tupla: tuple = listado.curselection()
    if tupla:
        elemento: str = listado.get(tupla[0])
        indice: int = 0
        while elemento[indice] != "k":
            indice += 1
        lista.pop(int(elemento[0:indice]))
        listado.delete("active", "active")
        listado.selection_clear(0, "end")
    else:
        showinfo("Baja", "Para borrar un tipo de camión, "
                 + "primero debe seleccionar uno.")


lista: dict[int, Decimal] = {}
inicio = Tk()
inicio.title("Integrador de Modelo y Simulación")
inicio.resizable(False, False)
Label(inicio, text="Datos",
      font=("Times New Roman", 28,
            "bold")).grid(row=0, column=0, columnspan=12)
Label(inicio, text="Camiones",
      font=("Unicode", 14,
            "roman")).grid(row=1, column=0, columnspan=2)
Label(inicio, text="Capacidad (kg)",
      font=("Arial", 12, "italic")).grid(row=2, column=0)
Label(inicio, text="Consumo ($/l)",
      font=("Arial", 12, "italic")).grid(row=4, column=0)
e_peso = Entry(inicio, highlightthickness=1,
               highlightbackground="black", highlightcolor="black")
e_precio = Entry(inicio, highlightthickness=1,
                 highlightbackground="black", highlightcolor="black")
c_peso = Label(inicio, text="", font=("Helvetica", 10),
               foreground="#ff0000", padx=0, pady=0)
c_precio = Label(inicio, text="", font=("Helvetica", 10),
                 foreground="#ff0000", padx=0, pady=0)
e_peso.grid(row=2, column=1)
e_precio.grid(row=4, column=1)
c_peso.grid(row=3, column=0, columnspan=2)
c_precio.grid(row=5, column=0, columnspan=2)
transportes = Frame(inicio)
barra_uno = Scrollbar(transportes)
listado = Listbox(transportes, height=9,
                  yscrollcommand=barra_uno.set)
barra_uno.config(command=listado.yview)
barra_uno.pack(side="right")
listado.pack(side="left", fill="x")
transportes.grid(row=7, column=0, rowspan=3, columnspan=2)
agregar = Button(inicio, text='Agregar', background="#eeeeee",
                 foreground="#000000", font=("Unicode", 14, "roman"),
                 command=partial(alta, e_peso, c_peso, e_precio,
                                 c_precio, lista, listado), pady=0)
eliminar = Button(inicio, text='Eliminar', background="#eeeeee",
                  foreground="#000000", font=("Unicode", 14, "roman"),
                  command=partial(baja, lista, listado), pady=0)
agregar.grid(row=6, column=0)
eliminar.grid(row=6, column=1)
Label(inicio, text="Cantidad de vacas",
      font=("Arial", 12, "italic")).grid(row=4, column=3, sticky="w")
Label(inicio, text="Sueldo de conductor ($)",
      font=("Arial", 12,
            "italic")).grid(row=1, column=2, columnspan=2, padx=20)
Label(inicio, text="Cantidad de marcas",
      font=("Arial", 12, "italic")).grid(row=1, column=4,
                                         columnspan=2, padx=20)
Label(inicio, text="Peso moda de vaca (kg)",
      font=("Arial", 12, "italic")).grid(row=4, column=5, sticky="w")
Label(inicio, text="Peso mínimo de vaca", pady=0,
      font=("Arial", 12, "italic")).grid(row=6, column=3, sticky="sw",
                                         ipady=0)
Label(inicio, text="Peso máximo de vaca", pady=0,
      font=("Arial", 12, "italic")).grid(row=6, column=5, sticky="s",
                                         ipady=0)
Label(inicio, text="Distancia a recorrer (km)",
      font=("Arial", 12,
            "italic")).grid(row=1, column=6, columnspan=2)
e1 = Entry(inicio, highlightthickness=1,
           highlightbackground="black", highlightcolor="black")
e2 = Entry(inicio, highlightthickness=1,
           highlightbackground="black", highlightcolor="black")
e3 = Entry(inicio, highlightthickness=1,
           highlightbackground="black", highlightcolor="black")
e4 = Entry(inicio, highlightthickness=1,
           highlightbackground="black", highlightcolor="black")
e5 = Entry(inicio, highlightthickness=1,
           highlightbackground="black", highlightcolor="black")
e6 = Entry(inicio, highlightthickness=1,
           highlightbackground="black", highlightcolor="black")
e7 = Entry(inicio, highlightthickness=1,
           highlightbackground="black", highlightcolor="black")
e1.grid(row=2, column=2, columnspan=2)
e2.grid(row=2, column=4, columnspan=2)
e3.grid(row=5, column=2, columnspan=2)
e4.grid(row=5, column=4, columnspan=2)
e5.grid(row=7, column=2, columnspan=2, sticky='n')
e6.grid(row=7, column=4, columnspan=2, sticky='n')
e7.grid(row=2, column=6, columnspan=2)
c1 = Label(inicio, text="", font=("Helvetica", 10),
           foreground="#ff0000", padx=0, pady=0)
c2 = Label(inicio, text="", font=("Helvetica", 10),
           foreground="#ff0000", padx=0, pady=0)
c3 = Label(inicio, text="", font=("Helvetica", 10),
           foreground="#ff0000", padx=0, pady=0)
c4 = Label(inicio, text="", font=("Helvetica", 10),
           foreground="#ff0000", padx=0, pady=0)
c5 = Label(inicio, text="", font=("Helvetica", 10),
           foreground="#ff0000", padx=0, pady=0)
c6 = Label(inicio, text="", font=("Helvetica", 10),
           foreground="#ff0000", padx=0, pady=0)
c7 = Label(inicio, text="", font=("Helvetica", 10),
           foreground="#ff0000", padx=0, pady=0)
c1.grid(row=3, column=2, columnspan=2)
c2.grid(row=3, column=4, columnspan=2)
c3.grid(row=6, column=2, columnspan=2, sticky="n")
c4.grid(row=6, column=4, columnspan=2, sticky="n")
c5.grid(row=7, column=2, columnspan=2, sticky="s")
c6.grid(row=7, column=4, columnspan=2, sticky="s")
c7.grid(row=3, column=6, columnspan=2)
conf = Label(inicio, text="", font=("Verdana", 12),
             foreground="#000000")
conf.grid(row=9, column=2, columnspan=4, sticky="s")
res = Label(inicio, text="Peso total de vacas:",
            font=("Cambria", 14))
res.grid(row=9, column=2, columnspan=4, sticky="n")
Label(inicio, text="Pesos (kg)",
      font=("Unicode", 14,
            "roman")).grid(row=4, column=6, columnspan=2)
ganado = Frame(inicio)
barra_dos = Scrollbar(ganado)
vacas = Listbox(ganado, height=12,
                yscrollcommand=barra_dos.set)
barra_dos.config(command=vacas.yview)
barra_dos.pack(side="right", expand=True, fill="y")
vacas.pack(side="left", expand=True, fill="both")
ganado.grid(row=5, column=6, rowspan=5, columnspan=2)


def calculo(ventana: Tk, dicc: dict[int, Decimal], e_sueldo: Entry,
            c_sueldo: Label, e_marcas: Entry, c_marcas: Label,
            e_vacas: Entry, c_vacas: Label, e_moda: Entry,
            c_moda: Label, e_minimo: Entry, c_minimo: Label,
            e_maximo: Entry, c_maximo: Label, e_distancia: Entry,
            c_distancia: Label, resultado: Label, confianza: Label,
            l_vacas: Listbox, /) -> None:
    c_sueldo.config(text="")
    c_marcas.config(text="")
    c_vacas.config(text="")
    c_moda.config(text="")
    c_minimo.config(text="")
    c_maximo.config(text="")
    c_distancia.config(text="")
    resultado.config(text="Peso total de vacas:")
    confianza.config(text="")
    for widget in ventana.winfo_children():
        if isinstance(widget, Toplevel):
            widget.destroy()
    l_vacas.delete(0, "end")
    sueldo: Decimal
    cant_marcas: int
    cant_vacas: int
    moda: Decimal
    minimo: Decimal
    maximo: Decimal
    distancia: Decimal
    invalido: bool = False
    s_sueldo: str = e_sueldo.get()
    try:
        sueldo = Decimal(s_sueldo)
        if sueldo <= 0:
            c_sueldo.config(text="Valor inválido como sueldo")
            invalido = True
    except:
        c_sueldo.config(text="No es un número")
        invalido = True
    s_marcas: str = e_marcas.get()
    try:
        cant_marcas = int(s_marcas)
        if cant_marcas <= 0:
            c_marcas.config(text="No es un número natural")
            invalido = True
    except:
        c_marcas.config(text="No es un número natural")
        invalido = True
    s_vacas: str = e_vacas.get()
    try:
        cant_vacas = int(s_vacas)
        if cant_vacas <= 0:
            c_vacas.config(text="No es un número natural")
            invalido = True
    except:
        c_vacas.config(text="No es un número natural")
        invalido = True
    s_moda: str = e_moda.get()
    try:
        moda = Decimal(s_moda)
        if moda <= 0:
            c_moda.config(text="Valor inválido como moda")
            invalido = True
    except:
        c_moda.config(text="No es un número")
        invalido = True
    s_minimo: str = e_minimo.get()
    try:
        minimo = Decimal(s_minimo)
        if minimo <= 0:
            c_minimo.config(text="Valor inválido como minimo")
            invalido = True
    except:
        c_minimo.config(text="No es un número")
        invalido = True
    s_maximo: str = e_maximo.get()
    try:
        maximo = Decimal(s_maximo)
        if maximo <= 0:
            c_maximo.config(text="Valor inválido como maximo")
            invalido = True
    except:
        c_maximo.config(text="No es un número")
        invalido = True
    s_distancia: str = e_distancia.get()
    try:
        distancia = Decimal(s_distancia)
        if distancia <= 0:
            c_distancia.config(text="Valor inválido como distancia")
            invalido = True
    except:
        c_distancia.config(text="No es un número")
        invalido = True
    if dicc == {}:
        showinfo("Cálculo", "Para calcular la cantidad de camiones, "
                 + "primero debe cargar al menos un tipo.")
        invalido = True
    if invalido:
        return
    sueldo = sueldo.quantize(Decimal("1.00"))
    if minimo < moda < maximo:
        estruct = Estructura(cant_vacas)  # 60383
        cf = congruencias_fundamental(estruct)
        if cf != None:
            if monobits(cf) and chi_cuadrado(cf) \
                    and poker(cf) and rachas(cf):
                confianza.config(text="La muestra de vacas es "
                                 + "suficientemente aleatoria")
            else:
                confianza.config(text="La muestra de vacas no es "
                                 + "suficientemente aleatoria")
            dict_marcas: dict[int, int] = {}
            clave: int
            for i in range(cant_marcas):
                dict_marcas[i] = 0
            lis_vacas: list[float] = []
            fc: Decimal = (moda-minimo)/(maximo-minimo)
            aux: float
            suma: float = 0.0
            for var in cf:
                if var[1] < fc:
                    aux = float(minimo) + sqrt(var[1]
                                               * float(maximo-minimo)
                                               * float(moda-minimo))
                else:
                    aux = float(maximo) - sqrt((1-var[1])
                                               * float(maximo-minimo)
                                               * float(maximo-moda))
                clave = round((((aux-float(minimo))
                                / float(maximo-minimo))
                               * cant_marcas)-0.5)
                dict_marcas[clave] = dict_marcas[clave] + 1
                lis_vacas.append(aux)
                l_vacas.insert("end", aux)
                suma += aux
            resultado.config(text="Peso total de vacas: "
                             + str(suma) + " kg.")
            salida: Toplevel = Toplevel(ventana, width=800,
                                        height=600)
            salida.title("Gráfico")
            salida.resizable(False, False)
            Label(salida, text="Resultados",
                  font=("Times New Roman", 20,
                        "bold")).pack(side="top")
            marcador = Frame(salida)
            barra_marcas = Scrollbar(marcador)
            lista_marcas = Listbox(marcador, height=15, width=35,
                                   yscrollcommand=barra_marcas.set)
            barra_marcas.config(command=lista_marcas.yview)
            barra_marcas.grid(row=0, column=1)
            lista_marcas.grid(row=0, column=0)
            Label(marcador, text="Probabilidad m:",
                  foreground="#000000").grid(row=1, column=0,
                                             columnspan=2)
            Label(marcador, text=str(2/(maximo-minimo)),
                  foreground="#00ff00").grid(row=2, column=0,
                                             columnspan=2)
            Label(marcador, text="a: Mínimo",
                  foreground="#000000").grid(row=3, column=0,
                                             columnspan=2)
            Label(marcador, text="b: Máximo",
                  foreground="#000000").grid(row=4, column=0,
                                             columnspan=2)
            Label(marcador, text="c: Moda",
                  foreground="#000000").grid(row=5, column=0,
                                             columnspan=2)
            Label(marcador, text="Dist. Teórica",
                  foreground="#0000ff").grid(row=6, column=0)
            Label(marcador, text="Dist. p/ Marca",
                  foreground="#ff0000").grid(row=6, column=1)
            Label(marcador, text="Eje X: Peso (kg)",
                  foreground="#000000").grid(row=7, column=0)
            Label(marcador, text="Eje Y: P(X)",
                  foreground="#000000").grid(row=7, column=1)
            tupla: tuple[tuple[int, Decimal], ...] = tuple(dicc.items())
            cantidad: int
            precio: Decimal
            peso: int = tupla[0][0]
            camiones: int = (Decimal(suma)
                             + (moda / 2)) // tupla[0][0]
            ideal: Decimal = camiones * (tupla[0][1] * distancia
                                         + sueldo)
            for i in range(1, len(tupla)):
                cantidad = (Decimal(suma) + (moda / 2)) // tupla[i][0]
                precio = cantidad * (tupla[i][1] * distancia + sueldo)
                if precio < ideal:
                    camiones = cantidad
                    peso = tupla[i][0]
                    ideal = precio
            Label(marcador, text="Tipo de camión ideal: " + str(peso)
                                 + "kg.",
                  foreground="#000000").grid(row=8, column=0,
                                             columnspan=2)
            Label(marcador, text="Cantidad de camiones: "
                                 + str(camiones),
                  foreground="#000000").grid(row=9, column=0,
                                             columnspan=2)
            Label(marcador, text="Costo: $" + str(ideal),
                  foreground="#000000").grid(row=10, column=0,
                                             columnspan=2)
            marcador.pack(side="right")
            for marca in dict_marcas.items():
                aux = float(minimo) + ((marca[0] + 0.5)
                                       * float(maximo - minimo)
                                       / cant_marcas)
                lista_marcas.insert("end", str(aux) + "kg: "
                                    + str(marca[1]) + " vacas")
            canvas = ScrolledCanvas(salida, width=600,
                                    height=400)
            canvas.pack(side="bottom")
            pantalla = TurtleScreen(canvas, delay=0)
            pantalla.screensize(canvwidth=600, canvheight=400)
            tortuga = RawTurtle(pantalla, visible=False)
            tortuga.resizemode("noresize")
            tortuga.up()
            tortuga.goto(-280, -180)
            tortuga.down()
            # Eje X
            tortuga.goto(280, -180)
            # Parte inferior de flecha X
            tortuga.goto(275, -185)
            tortuga.up()
            tortuga.goto(280, -180)
            tortuga.down()
            # Parte superior de flecha X
            tortuga.goto(275, -175)
            tortuga.up()
            tortuga.goto(-280, -180)
            tortuga.down()
            # Eje Y
            tortuga.goto(-280, 180)
            # Parte izquierda de flecha Y
            tortuga.goto(-285, 175)
            tortuga.up()
            tortuga.goto(-280, 180)
            tortuga.down()
            # Parte derecha de flecha Y
            tortuga.goto(-275, 175)
            tortuga.up()
            tortuga.goto(-289, -168)
            tortuga.down()
            # Etiqueta 0 en eje Y
            tortuga.goto(-287, -168)
            tortuga.goto(-286, -169)
            tortuga.goto(-286, -174)
            tortuga.goto(-287, -175)
            tortuga.goto(-289, -175)
            tortuga.goto(-290, -174)
            tortuga.goto(-290, -169)
            tortuga.goto(-289, -168)
            tortuga.up()
            tortuga.goto(-290, -180)
            tortuga.down()
            # Eje de punto 0 en Y
            tortuga.goto(-280, -180)
            tortuga.up()
            tortuga.goto(-290, 116)
            tortuga.down()
            # Etiqueta m (prob. max.) en eje Y
            tortuga.goto(-290, 110)
            tortuga.up()
            tortuga.goto(-290, 115)
            tortuga.down()
            tortuga.goto(-289, 116)
            tortuga.goto(-288, 116)
            tortuga.goto(-287, 115)
            tortuga.goto(-287, 110)
            tortuga.up()
            tortuga.goto(-287, 115)
            tortuga.down()
            tortuga.goto(-286, 116)
            tortuga.goto(-285, 116)
            tortuga.goto(-284, 115)
            tortuga.goto(-284, 110)
            tortuga.up()
            tortuga.goto(-290, 120)
            tortuga.down()
            # Eje de punto m en Y
            tortuga.goto(-270, 120)
            tortuga.up()
            tortuga.goto(-277, -187)
            tortuga.down()
            # Etiqueta a en eje X
            tortuga.goto(-275, -187)
            tortuga.goto(-274, -188)
            tortuga.goto(-274, -189)
            tortuga.goto(-276, -189)
            tortuga.goto(-277, -190)
            tortuga.goto(-276, -191)
            tortuga.goto(-274, -191)
            tortuga.goto(-274, -189)
            tortuga.up()
            tortuga.goto(-270, -190)
            tortuga.down()
            # Eje de punto a (min) en X
            tortuga.goto(-270, -170)
            tortuga.up()
            tortuga.goto(263, -184)
            tortuga.down()
            # Etiqueta b en eje X
            tortuga.goto(263, -191)
            tortuga.goto(265, -191)
            tortuga.goto(266, -190)
            tortuga.goto(266, -188)
            tortuga.goto(265, -187)
            tortuga.goto(263, -187)
            tortuga.up()
            tortuga.goto(270, -190)
            tortuga.down()
            # Eje de punto b (max) en X
            tortuga.goto(270, -170)
            tortuga.up()
            aux = -270 + (float(fc)*540)
            tortuga.goto(aux-4, -187)
            tortuga.down()
            # Etiqueta c en eje X
            tortuga.goto(aux-7, -187)
            tortuga.goto(aux-8, -188)
            tortuga.goto(aux-8, -190)
            tortuga.goto(aux-7, -191)
            tortuga.goto(aux-4, -191)
            tortuga.up()
            tortuga.goto(aux, -190)
            tortuga.down()
            # Eje de punto c (moda) en X
            tortuga.goto(aux, -170)
            tortuga.up()
            tortuga.color("#0000ff")
            tortuga.goto(-270, -180)
            tortuga.down()
            # Distribución teórica
            tortuga.goto(aux, 120)
            tortuga.goto(270, -180)
            tortuga.up()
            tortuga.color("#ff0000")
            tortuga.goto(-270, -180)
            fx: float
            tortuga.down()
            # Distribución por marca
            for marca in dict_marcas.items():
                aux = -270 + ((marca[0] + 0.5) * 540 / cant_marcas)
                fx = -180 + (marca[1] * cant_marcas * 300
                             / float(2 * cant_vacas))
                tortuga.goto(aux, fx)
            tortuga.goto(270, -180)
            tortuga.up()
        else:
            showinfo("Cálculo", "No se puede generar la cantidad de "
                     + "vacas aleatorias deseadas.")
    else:
        showinfo("Cálculo", "Los valores de los pesos de vacas "
                 + "no tienen sentido, se solapan.")


calcular = Button(inicio, text="Calcular", background="#eeeeee",
                  foreground="#000000", font=("Unicode", 14, "roman"),
                  command=partial(calculo, inicio, lista, e1, c1, e2,
                                  c2, e3, c3, e4, c4, e5, c5, e6, c6,
                                  e7, c7, res, conf, vacas))
calcular.grid(row=8, column=3, columnspan=3)
inicio.mainloop()
