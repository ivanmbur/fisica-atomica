import numpy as np
import matplotlib.pyplot as plt

#El intervalo en el que la ecuacion x'' + k*x=0  se va a resolver es [a, b].

a = 0.0
b = 40.0
t0 = 0.0

#x0 y v0 son las condiciones de Cauchy en el instante t0.

def solve(x0, v0, k, dt):
	n = 1
	T = [t0]
	X = [x0]
	V = [v0]
	while (T[len(T) - 1] < b):
		V.append(V[n - 1] - k*X[n - 1]*dt)
		X.append(X[n - 1] + V[n]*dt)
		T.append(T[n - 1] + dt)
		n = n + 1
	return T, X

T1, X1 = solve(0, 1, 1, 0.01)
k1 = 1.0
T2, X2 = solve(0, 1, 50, 0.01)
k2 = 10.0
T3, X3 = solve(0, 1, 0.3, 0.01)
k3 = 0.5
T4, X4 = solve(0, 1, 10, 0.01)
k4 = 5.0

T1 = np.array(T1)
T2 = np.array(T2)
T3 = np.array(T3)
T4 = np.array(T4)

plt.figure()
plt.title("Solution of x''(t)+kx(t)=0 with x(0)=0 and v(0)=1")
plt.scatter(T3, X3, s = 1, label = r"$k=0.5$")
plt.scatter(T1, X1, s = 1, label = r"$k=1$")
plt.scatter(T4, X4, s = 1, label = r"$k=5$")
plt.scatter(T2, X2, s = 1, label = r"$k=10$")
plt.xlabel('$t$')
plt.ylabel("x")
plt.legend(fontsize = 10)
plt.savefig("tarea1.png")

#Observaciones: decidi graficar contra tiempo ya que facilita la interpretacion en terminos de cantidades fisicas. Jugar con el valor de k nos muestra que a medida que este aumenta tambien lo hace la frecuencia y disminuye la amplitud. Para entender esto, notemos que la ecuacion que queremos resolver es la de una masa de 1kg bajo la accion de un resorte de rigidez k. En efecto, cuando la rigidez aumenta se espera que la frecuencia de oscilacion aumente y la amplitud disminuya ya que es mas dificil deformar el resorte. Por otra parte, es importante notar que la solucion numerica entrega el comportamiento armonico esperado de la solucion analitica.
