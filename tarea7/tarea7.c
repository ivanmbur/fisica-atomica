#include <stdio.h>
#include <math.h>
#include <stdlib.h>

double *integrador_final(double *u, double psi_0, double *V, double epsilon_derivada, double minimo_derivada, double maximo_derivada, double epsilon_energia, double minimo_energia, double maximo_energia, double delta);

int main()
{
	/*Se definen las constantes principales*/
	
	/*Numero atómico del átomo de dos electrones*/
	const int Z = 2;

	/*Numero de pasos en la integración*/
	const int N_PASOS = 100000; 

	/*Punto de inicio*/
	const float U_0 = 0.00001;

	/*Punto final*/
	const float U_N = 40;

	/*Numero de veces que se va a iterar el metodo de Hartree*/
	const int N_ITERACIONES = 10;

	/*Minima derivada inicial de las funciones de onda*/
	const double MINIMO_DERIVADA = -5;

	/*Maxima derivada inicial de las funciones de onda*/
	const double MAXIMO_DERIVADA = 0;

	/*Criterio derivadas*/
	const double EPSILON_DERIVADA = 0.01;
	
	/*Minima energia que se considera*/
	const double MINIMO_ENERGIA = -1;

	/*Maxima energia que se considera*/
	const double MAXIMO_ENERGIA = -0.3;

	/*Finura de la busqueda de energias*/
	const double DELTA_ENERGIA = 0.0001;
	
	/*Preción energias*/
	const double EPSILON_ENERGIA = 0.0000000001;

	/*Valor inicial de la funcion de onda*/
	const double PSI_0 = 1;

	/*Se discretiza la posición*/
	double *u = malloc(N_PASOS*sizeof(double));
	int i;
	for(i=0;i<N_PASOS;i++)
	{
		u[i] = pow(10,log10(U_0) + ((log10(U_N)-log10(U_0))*i/(N_PASOS-1)));
	} 

	/*Se inicializa el potencial*/
	double *V = malloc(sizeof(u));
	int j;
	for(j=0;j<N_PASOS;j++)
	{
		V[j] = -2.0*(1.0+((Z-1.0)*exp(-u[j])))/(Z*u[j]);
	}	
	
	double *psi = integrador_final(u, PSI_0, V, EPSILON_DERIVADA, MINIMO_DERIVADA, MAXIMO_DERIVADA, EPSILON_ENERGIA, MINIMO_ENERGIA, MAXIMO_ENERGIA, DELTA_ENERGIA); 
	
	return 0;
}

/*Calcula la segunda derivada de la función de onda dado el punto en el que se quiere calcular u, su valor en el punto anterior psi, su derivada en el punto anterior Dpsi, el potencial en el punto V y la energía de la solucion E*/
double D2psi(double u, double psi, double Dpsi, double V, double E)
{
	return ((psi*(V-E)) - (2*Dpsi/u));
} 

/*Calcula una solucion a la ecuacion de Schrödinger radial dadas la posiciones u, las condiciones iniciales sobre la función psi y su derivada Dpsi, el potencial V, la energía de la solucion E y el numero de valores que se quiere calcular n*/	
double * integrador(double *u, double psi_0, double Dpsi_0, double *V, double E, int n)
{
	double *solucion = malloc(n*sizeof(float));
	double Dsolucion = Dpsi_0;
	solucion[0] = psi_0;
	int k;
	for(k=1;k<n;k++)
	{	
		double delta_u = u[k] - u[k-1];
		double D2solucion = D2psi(u[k-1], solucion[k-1], Dsolucion, V[k-1], E);
		Dsolucion = Dsolucion + (D2solucion*delta_u);
		solucion[k] = solucion[k-1] + (Dsolucion*delta_u);
	}
	return solucion;
}

/*Calcula la primera derivada de la función de onda que la hace continua dadas el potencial V y la energía E, la distancia entre las derivadas a la cual se considera aceptable el resultado epsilon, el minimo valor que se considera minimo y el maximo valor que se considera maximo. El metodo que se utiliza es de biseccion. Vale la pena decir que si epsilon es muy pequeño este metodo no va a converger pues en efecto la solucion tiene un cambio en las derivadas*/
double Dpsi(double *u, double psi_0, double *V, double E, double epsilon, double minimo, double maximo)
{
	double inferior = minimo;
	double superior = maximo;
	double *resultado = integrador(u, psi_0, (superior + inferior)/2, V, E, 3); 
	double distancia = sqrt((((resultado[2]-resultado[1])/(u[2] - u[1]))-((resultado[1]-resultado[0])/(u[1] - u[0])))*(((resultado[2]-resultado[1])/(u[2] - u[1]))-((resultado[1]-resultado[0])/(u[1] - u[0])))); 
	double distancia1;
	double distancia2;
	do
	{
		double mitad = (superior + inferior)/2;
		resultado = integrador(u, psi_0, (superior + mitad)/2, V, E, 3); 
		distancia1 = sqrt((((resultado[2]-resultado[1])/(u[2] - u[1]))-((resultado[1]-resultado[0])/(u[1] - u[0])))*(((resultado[2]-resultado[1])/(u[2] - u[1]))-((resultado[1]-resultado[0])/(u[1] - u[0]))));
		resultado = integrador(u, psi_0, (inferior + mitad)/2, V, E, 3); 
		distancia2 = sqrt((((resultado[2]-resultado[1])/(u[2] - u[1]))-((resultado[1]-resultado[0])/(u[1] - u[0])))*(((resultado[2]-resultado[1])/(u[2] - u[1]))-((resultado[1]-resultado[0])/(u[1] - u[0])))); 
		if(distancia1<distancia2)
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
	while(distancia>epsilon);
	free(resultado);
	return (superior + inferior)/2;
}		

/*Calcula una solucion a la ecuacion de Schrödinger radial dadas el potencial V y la energía de la solucion E*/	
double *integrador1(double *u, double psi_0, double *V, double E, double epsilon, double minimo, double maximo)
{
	return integrador(u, psi_0, Dpsi(u, psi_0, V, E, epsilon, minimo, maximo), V, E, sizeof(u)/sizeof(double)); 
}

/*Encuentra energias dado el potencial V, la energia minima a considerar minimo, la maxima maximo y la finura con la cual se van a buscar energias delta*/
double *energias(double *u, double psi_0, double *V, double epsilon_derivada, double minimo_derivada, double maximo_derivada, double minimo_energia, double maximo_energia, double delta)
{
	int n = sizeof(u)/sizeof(double);
	double *encontradas = malloc(50*sizeof(double));
	int posicion = 0;
	double pasado = minimo_energia;
	double *psi_pasado = integrador1(u, psi_0, V, pasado, epsilon_derivada, minimo_derivada, maximo_derivada);
	double paso_pasado = psi_pasado[n - 1] - psi_pasado[n - 2];
	double E;
	double paso_presente;
	double *psi_presente;
	for(E=pasado+delta;E<maximo_energia;E = E+delta)
	{	
		psi_presente = integrador1(u, psi_0, V, E, epsilon_derivada, minimo_derivada, maximo_derivada);
		paso_presente = psi_presente[n - 1] - psi_presente[n - 2];
		if(paso_presente/sqrt(paso_presente*paso_presente) != paso_pasado/sqrt(paso_pasado*paso_pasado))
		{
			encontradas[posicion] = pasado;
			encontradas[posicion + 1] = E;
			posicion = posicion + 2;
		}
		pasado = E;
		paso_pasado = paso_presente;
	}
	free(psi_pasado);
	free(psi_presente);
	return encontradas;
}
		

/*Refina la energia dado el potencial V, la precision epsilon, la energia minima a considerar minimo y la maxima maximo*/
double refinador(double *u, double psi_0, double *V, double epsilon_derivada, double minimo_derivada, double maximo_derivada, double epsilon_energia, double minimo_energia, double maximo_energia)
{
	int n = sizeof(u)/sizeof(double);
	double superior = maximo_energia;
	double inferior = minimo_energia;
	double distancia = superior - inferior;
	double *psi_superior;
	double *psi_mitad;
	do
	{
		double mitad = (superior + inferior)/2;
		psi_superior = integrador1(u, psi_0, V, superior, epsilon_derivada, minimo_derivada, maximo_derivada);
		psi_mitad = integrador1(u, psi_0, V, mitad, epsilon_derivada, minimo_derivada, maximo_derivada); 
		double paso_superior = psi_superior[n - 1] - psi_superior[n - 2];
		double paso_mitad = psi_mitad[n - 1] - psi_mitad[n - 2];
		
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
	while(distancia > epsilon_energia);
	free(psi_superior);
	free(psi_mitad);
	return (inferior + superior)/2;
}

/*Calcula el estado base dado un potencial V*/
double *integrador_final(double *u, double psi_0, double *V, double epsilon_derivada, double minimo_derivada, double maximo_derivada, double epsilon_energia, double minimo_energia, double maximo_energia, double delta)
{
	double *encontradas = energias(u, psi_0, V, epsilon_derivada, minimo_derivada, maximo_derivada, minimo_energia, maximo_energia, delta);
	double E = refinador(u, psi_0, V, epsilon_derivada, minimo_derivada, maximo_derivada, epsilon_energia, encontradas[0], encontradas[1]);
	return integrador1(u, psi_0, V, E, epsilon_derivada, minimo_derivada, maximo_derivada);
}

