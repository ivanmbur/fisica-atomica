import numpy as np
import matplotlib.pyplot as plt

plt.rc("text", usetex = True)

dx = 100000
res = 100
E_ligado = np.linspace(0, 18, res)
inicial = 1
criterio = 0.03
Epar = []
Eimpar = []

def U(x):
	return x**2

def psi(E, paridad, inicial):
	x = np.linspace(0, np.sqrt(E) + 3, dx)
	if (paridad == 0):
		phi = [inicial]
		Dphi = [0]
	else:
		phi = [0]
		Dphi = [inicial]
	for n in range(1,len(x)):
			Dphi.append(Dphi[n - 1] + ((U(x[n - 1]) - E) * phi[n - 1] * (x[n] - x[n - 1])))
			phi.append(phi[n - 1] + (Dphi[n - 1] * (x[n] - x[n - 1])))
	return [phi, Dphi]

def buscador(E_iniciales, paridad):
	E = []
	for n in range(0, len(E_iniciales) - 1):
		Dsol1 =  psi(E_iniciales[n], paridad, inicial)[1]
		Dsol2 =  psi(E_iniciales[n + 1], paridad, inicial)[1]
		if(np.sign(Dsol1[dx-1]) != np.sign(Dsol2[dx-1])):
			E.append(E_iniciales[n])	
			E.append(E_iniciales[n + 1])	
	return E

#plt.plot(x, psi(x, 0.91, 0, inicial)[0])
#plt.show()

Epar = buscador(E_ligado, 0)
Eimpar = buscador(E_ligado, 1)
print Epar
print Eimpar


for n in range(0, len(Epar[::2])):
	valor = criterio + 1
	while (abs(valor) > criterio):
		Energias = np.linspace(Epar[2*n], Epar[2*n + 1], res)
		[Epar[2*n], Epar[2*n + 1]] = buscador(Energias, 0)
		valorprevio = valor
		valor = psi((Epar[2*n] + Epar[2*n + 1])/2, 0, inicial)[0][dx-1]
		print valor
		if valorprevio == valor:
			print n, Epar[2*n], Epar[2*n+1]
			break

for n in range(0, len(Eimpar[::2])):
	valor = criterio + 1
	while (abs(valor) > criterio):
		Energias = np.linspace(Eimpar[2*n], Eimpar[2*n + 1], res)
		[Eimpar[2*n], Eimpar[2*n + 1]] = buscador(Energias, 1)
		valorprevio = valor
		valor = psi((Eimpar[2*n] + Eimpar[2*n + 1])/2, 1, inicial)[0][dx-1]
		print valor
		if valorprevio == valor:
			print n, Epar[2*n], Epar[2*n+1]
			break

plt.figure()
plt.title(r"Comportamiento en estados ligados de part\'iculas pares")
for n in range(0, len(Epar) - 1):
	if ((2*n + 1) < len(Epar)):
		E = (Epar[2*n + 1] + Epar[2*n])/2
		phi = psi(E, 0, inicial)[0]
		plt.plot(np.linspace(0, np.sqrt(E) + 3, dx), phi, label = r"$%f \leq E \leq %f$" % (Epar[2*n], Epar[2*n + 1]))
	else:
		break
plt.ylabel(r"$\psi$", fontsize = 15)
plt.xlabel(r"$x$", fontsize = 15)
plt.legend()
plt.savefig("comportamiento_ligado_par.pdf", format = "pdf") 
plt.figure()
plt.title(r"Comportamiento en estados ligados de part\'iculas impares")
for n in range(0, len(Eimpar) - 1):
	if ((2*n + 1) < len(Eimpar)):
		E = (Eimpar[2*n + 1] + Eimpar[2*n])/2
		phi = psi(E, 1, inicial)[0]
		plt.plot(np.linspace(0, np.sqrt(E) + 3, dx), phi, label = r"$%f \leq E \leq %f$" % (Eimpar[2*n], Eimpar[2*n + 1]))
	else:
		break
plt.ylabel(r"$\psi$", fontsize = 15)
plt.xlabel(r"$x$", fontsize = 15)
plt.legend()
plt.savefig("comportamiento_ligado_impar.pdf", format = "pdf") 
