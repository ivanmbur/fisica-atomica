import numpy as np
import matplotlib.pyplot as plt

plt.rc("text", usetex = True)

x = np.linspace(0, 3, 300 + 1)
U0 = 100
res = 100
E_ligado = np.linspace (0, 100, res)
inicial = 1
criterio = 0.014
Epar = []
Eimpar = []

def U(x):
	if(x <= 1):
		return 0
	else:
		return U0

def psi(x, U0, E, paridad, inicial):
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
		Dsol1 =  psi(x, U0, E_iniciales[n], paridad, inicial)[1]
		Dsol2 =  psi(x, U0, E_iniciales[n + 1], paridad, inicial)[1]
		if(np.sign(Dsol1[300]) != np.sign(Dsol2[300])):
			E.append(E_iniciales[n])	
			E.append(E_iniciales[n + 1])	
	return E

Epar = buscador(E_ligado, 0)
Eimpar = buscador(E_ligado, 1)

for n in range(0, len(Epar[::2])):
	valor = criterio + 1
	while (abs(valor) > criterio):
		Energias = np.linspace(Epar[2*n], Epar[2*n + 1], res)
		[Epar[2*n], Epar[2*n + 1]] = buscador(Energias, 0)
		valor = psi(x, U0, (Epar[2*n] + Epar[2*n + 1])/2, 0, inicial)[0][300]
		print valor

for n in range(0, len(Eimpar[::2])):
	valor = criterio + 1
	while (abs(valor) > criterio):
		Energias = np.linspace(Eimpar[2*n], Eimpar[2*n + 1], res)
		[Eimpar[2*n], Eimpar[2*n + 1]] = buscador(Energias, 1)
		valor = psi(x, U0, (Eimpar[2*n] + Eimpar[2*n + 1])/2, 1, inicial)[0][300]
		print valor

plt.figure()
plt.title(r"Comportamiento en estados ligados de part\'iculas pares")
for n in range(0, len(Epar) - 1):
	if ((2*n + 1) < len(Epar)):
		E = (Epar[2*n + 1] + Epar[2*n])/2
		phi = psi(x, U0, E, 0, inicial)[0]
		plt.plot(x, phi, label = r"$%f \leq E \leq %f$" % (Epar[2*n], Epar[2*n + 1]))
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
		phi = psi(x, U0, E, 1, inicial)[0]
		plt.plot(x, phi, label = r"$%f \leq E \leq %f$" % (Eimpar[2*n], Eimpar[2*n + 1]))
	else:
		break
plt.ylabel(r"$\psi$", fontsize = 15)
plt.xlabel(r"$x$", fontsize = 15)
plt.legend()
plt.savefig("comportamiento_ligado_impar.pdf", format = "pdf") 

f, axarr = plt.subplots(2, 2)
axarr[0, 0].set_title(r"$E = 100$")
axarr[0, 0].plot(x, psi(x, U0, 100, 0, inicial)[0], label = "par")
axarr[0, 0].plot(x, psi(x, U0, 100, 1, inicial)[0], label = "impar")
axarr[0, 0].legend()
axarr[0, 0].set_xlabel(r"Posici\'on")
axarr[0, 0].set_ylabel(r"Funci\'on de onda")
axarr[0, 1].set_title(r"$E = 200$")
axarr[0, 1].plot(x, psi(x, U0, 200, 0, inicial)[0], label = "par")
axarr[0, 1].plot(x, psi(x, U0, 200, 1, inicial)[0], label = "impar")
axarr[0, 1].set_xlabel(r"Posici\'on")
axarr[0, 1].set_ylabel(r"Funci\'on de onda")
axarr[1, 0].set_title(r"$E = 700$")
axarr[1, 0].plot(x, psi(x, U0, 700, 0, inicial)[0], label = "par")
axarr[1, 0].plot(x, psi(x, U0, 700, 1, inicial)[0], label = "impar")
axarr[1, 0].set_xlabel(r"Posici\'on")
axarr[1, 0].set_ylabel(r"Funci\'on de onda")
axarr[1, 1].set_title(r"$E = 1000$")
axarr[1, 1].plot(x, psi(x, U0, 1000, 0, inicial)[0], label = "par")
axarr[1, 1].plot(x, psi(x, U0, 1000, 1, inicial)[0], label = "impar")
axarr[1, 1].set_xlabel(r"Posici\'on")
axarr[1, 1].set_ylabel(r"Funci\'on de onda")
f.tight_layout()
plt.savefig("comportamiento_libre.pdf")

