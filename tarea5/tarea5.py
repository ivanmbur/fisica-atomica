import numpy as np
import matplotlib.pyplot as plt

#Entrega la solucion a la parte radial de la ecuacion de Schrodinger en coordenadas naturales del problema. Se integra en el intervalo [a, b] con condiciones iniciales (R0, DR0) con un salto dx. Se tiene como entrada tambien el valor de l y E. La salida son tres arrays, el primero con los valores de u, el segundo con los valores de la solucion y el ultimo con los de su derivada. La integracion se hace con el metodo de Euler. 

def integrador(a, b, dx, R0, DR0, l, E):
	R = [R0]
	DR = [DR0]
	res = int((b - a)/dx + 1)
	U = np.linspace(a, b, res)
	for n in range(1, res):
		DR.append(DR[n - 1] - ((2/U[n])*DR[n - 1] + (E + (2/U[n]) - (l*(l+1)/(U[n]**2)))*R[n - 1])*dx)
		R.append(R[n - 1] + DR[n]*dx)
	return U, np.array(R), np.array(DR)

#Se va a definir de una vez por todas el intervalo de integracion y la logitud de paso.

a = 0
b = 50
dx = 0.01

#El codigo inferior (que esta comentado para el funcionamiento del codigo posterior) muestra los primero diez valores de la derivada para una solucion con R0=1, l=0 y E=-1. Jugando con distintos de DR0 se encontro que para valores encima de cierto numero la primera derivada tenia un pico hacia abajo pronunciado en su segundo valor y para los inferiores a cierto numero este se volvia un pico hacia abajo. Para obtener una solucion suave es necesario examinar cual es este cierto numero. Ensayando distintos valores es facil notar que es cercano a -1. Sin embargo, es necesario notar que este valor si depende de la energia. Por lo tanto se deberia desarrollar un algoritmo que haga la busqueda de DR0 para cada valor de energia. Esto es justamente lo que se desarrolla abajo 
#U, R, DR = integrador(a, b, dx, 1, -5, 0, -1)
#plt.scatter(U[range(0,10)], DR[range(0,10)]) 
#plt.show()

#Basados en las observaciones superiores, este codigo busca el mejor valor de DR0 para una solucion que empieza en a, tiene longitud de paso dx, tiene una condicion inicial R0 y una energia E. El codigo ademas necesita unos valores iniciales entre los cuales se espera que este el valor buscado DR0_min y DR0_max, y la resolucion (el numero de pasos) res con la que se va a buscar. Finalmente se necesita el criterio despues del cual el buscador puede dejar de buscar. 
def DR0(a, dx, R0, E, DR0_min, DR0_max, res, criterio):
	Encontrado = False
	while Encontrado == False:
		possible_DR0 = np.linspace(DR0_min, DR0_max, res)
		encontrado = False
		n = 1
		while (encontrado == False):
			U, R, DR_1 = integrador(a, a + 2*dx, dx, R0, possible_DR0[n], 0, E)
			U, R, DR_2 = integrador(a, a + 2*dx, dx, R0, possible_DR0[n+1], 0, E)
			if (np.sign((DR_1[2] - DR_1[1])-(DR_1[1]-DR_1[0])) != np.sign((DR_2[2] - DR_2[1])-(DR_2[1]-DR_2[0]))):
				encontrado = True
				DR0_min = possible_DR0[n]
				DR0_max = possible_DR0[n + 1]
			else:
				n += 1
		if(abs(DR0_min - DR0_max) < criterio):
			Encontrado = True	
	return (DR0_min + DR0_max)/2

#Ya que lo siguiente es para el caso l=0, se definen R0, el criterio y la resolucion.
R0 = 1
res = 10000
criterio = 0.001

#Se ha encontrado que los valores de DR0 nunca salieron de este rango
DR0_min = -2
DR0_max = -0.5

#Se hizo una pequena prueba con resultados correctos
#print DR0(a, dx, R0, -1, DR0_min, DR0_max, res, criterio)

E = -0.5
testigo = False
previo = 100000000
while testigo == False and E > -1.2:
	U, R, DR = integrador(a, b, dx, R0, DR0(a, dx, R0, E, DR0_min, DR0_max, res, criterio), 0, E)
	if(R[len(R) - 1] < previo):
		testigo == True
	else:
		E += -0.1
	previo = R[len(R) - 1]
		
print E

