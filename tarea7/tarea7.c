#include <stdio.h>
#include <math.h>
#include <stdlib.h>

/*Las funciones que necesita el código*/
double D2psi(double *u, double psi_0, double Dpsi_0, double *V, double E);
double *integrador(double *u, double psi_0, double Dpsi_0, double *V, double E, int n);
double Dpsi(double *V, double E, double epsilon, double minimo, double maximo);
double *integrador1(double *V, double E);
double *energias(double *V, double minimo, double maximo, double deltas);
double refinador(double *V, double epsilon, double minimo, double maximo);
double *integrador_final(double *V);

int main()
{
	/*Se definen las constantes principales*/
	
	/*Numero atómico del átomo de dos electrones*/
	int Z = 2;

	/*Numero de pasos en la integración*/
	int n_pasos = 100000; 

	/*Punto de inicio*/
	float u_0 = 0.00001;

	/*Punto final*/
	float u_n = 40;

	/*Numero de veces que se va a iterar el metodo de Hartree*/
	int n_iteraciones = 10;

	/*Minima derivada inicial de las funciones de onda*/
	double minimo_derivada = -5;

	/*Maxima derivada inicial de las funciones de onda*/
	double maximo_derivada = 0;

	/*Criterio derivadas*/
	double epsilon_derivada = 0.01;
	
	/*Minima energia que se considera*/
	double minimo_energia = -1;

	/*Maxima energia que se considera*/
	double maximo_energia = -0.3;

	/*Finura de la busqueda de energias*/
	double delta_energia = 0.0001;
	
	/*Preción energias*/
	double epsilon_energia = 0.0000000001;

	/*Valor inicial de la funcion de onda*/
	double psi_0 = 1;

	/*Se discretiza la posición*/
	double *u = malloc(n_pasos*sizeof(double));
	int i;
	for(i=0;i<n_pasos;i++)
	{
		double u[i] = pow(10,log10(u_0) + ((log10(u_n)-log10(u_0))*i/(n_pasos-1)));
	} 

	/*Se inicializa el potencial*/
	double *V = malloc(n_pasos*sizeof(double));
	int j;
	for(j=0;j<n_pasos;j++)
	{
		V[j] = -2.0*(1.0+((Z-1.0)*exp(-u[j])))/(Z*u[j]);
	}	
	
	/*double *psi = integrador_final(V);*/
	int m;
	for(m=0;m<n_pasos;m++)
	{
		printf("%f %f\n", u[m], V[m]);
	} 
	
	return 0;
}

/*Calcula la segunda derivada de la función de onda dado el punto en el que se quiere calcular u, su valor en el punto anterior psi, su derivada en el punto anterior Dpsi, el potencial en el punto V y la energía de la solucion E*/
double D2psi(u, psi_0, Dpsi_0, V, E)
{
	return ((psi_0*(V-E)) - (2*Dpsi_0/u));
} 

/*Calcula una solucion a la ecuacion de Schrödinger radial dadas la posiciones u, las condiciones iniciales sobre la función psi y su derivada Dpsi, el potencial V, la energía de la solucion E y el numero de valores que se quiere calcular n*/	
double *integrador(u, psi_0, Dpsi_0, V, E, n)
{
	double *solucion = malloc(n*sizeof(float));
	double Dsolucion = Dpsi_0;
	solucion[0] = psi_0;
	int k;
	for(k=1;k<n;k++)
	{	
		double delta_u = u[k] - u[k-1];
		double D2solucion = D2psi(u[k], psi[k-1], Dsolucion, V, E);
		Dsolucion = Dsolucion + (D2solucion*delta_u);
		solucion[k] = solucion[k-1] + (Dsolucion*delta_u);
	}
	return solucion;
}

/*Calcula la primera derivada de la función de onda que la hace continua dadas el potencial V y la energía E, la distancia entre las derivadas a la cual se considera aceptable el resultado epsilon, el minimo valor que se considera minimo y el maximo valor que se considera maximo. El metodo que se utiliza es de biseccion. Vale la pena decir que si epsilon es muy pequeño este metodo no va a converger pues en efecto la solucion tiene un cambio en las derivadas*/
double Dpsi(V, E, epsilon, minimo, maximo)
{
	double inferior = minimo;
	double superior = maximo;
	double *resultado = integrador(u, psi_0, (superior + inferior)/2, V, E, 3); 
	double distancia = sqrt((((resultado[2]-resultado[1])/(u[2] - u[1]))-((resultado[1]-resultado[0])/(u[1] - u[0])))*(((resultado[2]-resultado[1])/(u[2] - u[1]))-((resultado[1]-resultado[0])/(u[1] - u[0])))); 
	do
	{
		double mitad = (superior + inferior)/2;
		resultado = integrador(u, psi_0, (superior + mitad)/2, V, E, 3); 
		distancia1 = sqrt((((resultado[2]-resultado[1])/(u[2] - u[1]))-((resultado[1]-resultado[0])/(u[1] - u[0])))*(((resultado[2]-resultado[1])/(u[2] - u[1]))-((resultado[1]-resultado[0])/(u[1] - u[0]))));
		resultado = integrador(u, psi_0, (inferior + mitad)/2, V, E, 3); 
		distancia2 = sqrt((((resultado[2]-resultado[1])/(u[2] - u[1]))-((resultado[1]-resultado[0])/(u[1] - u[0])))*(((resultado[2]-resultado[1])/(u[2] - u[1]))-((resultado[1]-resultado[0])/(u[1] - u[0])))); 
		if(distancia1<distacia2)
		{
			inferior = mitad;
			distancia = distancia1;
		}	
		else
		{
			superior = mitad;
			distancia = distancia2;
		}
	}
	while(distancia>epsilon)
	return (superior + inferior)/2;
}		

/*Calcula una solucion a la ecuacion de Schrödinger radial dadas el potencial V y la energía de la solucion E*/	
double *integrador1(V, E)
{
	return integrador(u, psi, Dpsi(V, E, epsilon_derivada, minimo_derivada, maximo_derivada), V, E, n_points);  
}

/*Encuentra energias dado el potencial V, la energia minima a considerar minimo, la maxima maximo y la finura con la cual se van a buscar energias delta*/
double *energias(V, minimo, maximo, delta)
{
	double *encontradas = malloc(50*sizeof(double));
	int posicion = 0;
	double pasado = minimo;
	double *psi_pasado = integrador1(V, pasado);
	double paso_pasado = psi_pasado[n_pasos - 1] - psi_pasado[n_pasos - 2];
	for(E=pasado+delta;E<maximo;E = E+delta)
	{	
		double *psi_presente = integrador1(V, E);
		paso_presente = psi_presente[n_pasos - 1] - psi_presente[n_pasos - 2];
		if(paso_presente/sqrt(paso_presente*paso_presente) != paso_pasado/sqrt(paso_pasado*paso_pasado))
		{
			encontradas[posicion] = pasado;
			encontradas[psicion + 1] = presente;
			posicion = posicion + 2;
		}
		pasado = presente;
		paso_pasado = paso_presente;
	}
	return energias;
}
		

/*Refina la energia dado el potencial V, la precision epsilon, la energia minima a considerar minimo y la maxima maximo*/
double refinador(V, epsilon, minimo, maximo)
{
	double superior = maximo;
	double inferior = minimo;
	double distancia = superior - inferior;
	do
	{
		double mitad = (superior + inferior)/2;
		double *psi_superior = integrador1(V, superior);
		double *psi_mitad = integrador1(V, mitad); 
		double paso_superior = psi_superior[n_pasos - 1] - psi_superior[n_pasos - 2];
		double paso_mitad = psi_mitad[n_pasos - 1] - psi_mitad[n_pasos - 2];
		
		if(paso_superior/sqrt(paso_superior*paso_superior) != paso_mitad/sqrt(paso_mitad*paso_mitad))
		{
			inferior = mitad;
		}
		else
		{
			superior = mitad;
		}
		distancia = superior - inferior;
	}
	while(distancia > epsilon)
	return (inferior + superior)/2;
}

/*Calcula el estado base dado un potencial V*/
double *integrador_final(V)
{
	double *encontradas = energias(V, minimo_energia, maximo_energia, delta_energia);
	double E = refinador(V, epsilon_energia, encontradas[0], encontradas[1]);
	return integrador1(V, E);
}

