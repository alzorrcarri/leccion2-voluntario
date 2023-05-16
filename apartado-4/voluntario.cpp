/* SIMULACIÓN DEL MODELO DE HOPFIELD DE RED NEURONAL USANDO EL MODELO DE ISING CON UN SOLO PATRÓN */
/* EJERCICIO VOLUNTARIO 2 CUARTO APARTADO DE LA LECCIÓN 2 DE FÍSICA COMPUTACIONAL */

#include <iostream>
#include <cmath>
#include <fstream>
#include <stdio.h>
#include <cstdlib>
#include "gsl_rng.h"

#define FILAS 20
#define COLUMNAS 20
#define mu 30

using namespace std;

gsl_rng *tau;

int main (void)
{
    /* Definimos inicialmente algunas variables necesarias. */
    int i,j,k,l,m,q,r,pasos,memoria;
    double h,p,temp,aux,omega,s;

    temp=0.0001;
    pasos=300*FILAS*COLUMNAS;
    s=0.01; // Sesgo en la generación de patrones aleatorios.

    double a[mu],solapamiento[mu];
    short int matriz[FILAS][COLUMNAS];
    short int patron[FILAS][COLUMNAS][mu];
    double umbral[FILAS][COLUMNAS];
    ofstream memoria_fich;

    memoria_fich.open("/Users/alzorrcarri/Documents/cphys/leccion2-voluntario/apartado-4/memoria.dat");

    extern gsl_rng *tau; //Puntero al estado del número aleatorio
    int semilla=123456; //Semilla del generador de números aleatorios

    tau=gsl_rng_alloc(gsl_rng_taus); //Inicializamos el puntero
    gsl_rng_set(tau,semilla); //Inicializamos la semilla

    /* Generamos los patrones aleatorios y los almacenamos en las matrices 'patron' */
    for ( m = 0; m <= mu-1; m++)
    {
        for ( i = 0; i <= FILAS-1; i++)
        {
            for ( j = 0; j <= COLUMNAS-1; j++)
            {
                patron[i][j][m]=0;
                aux=gsl_rng_uniform(tau);
                if (aux<=s)
                {
                    patron[i][j][m]=1;
                }
            }
                
        }
    }

    /* Calculamos el vector a para cada patrón */
    for ( m = 0; m < mu; m++)
    {
        a[m]=0.0;
        for ( i = 0; i <= FILAS-1; i++)
        {
            for ( j = 0; j <= COLUMNAS-1; j++)
            {
                a[m]=a[m]+patron[i][j][m];
            }
        }
        a[m]=a[m]/(1.0*FILAS*COLUMNAS);
    }

    /* Calculamos la matriz umbral de disparo (umbral) */
    for ( i = 0; i <= FILAS-1; i++)
    {
        for ( j = 0; j <= COLUMNAS-1; j++)
        {
            umbral[i][j]=0.0;
            for ( k = 0; k <= FILAS-1; k++)
            {
                for ( l = 0; l <= COLUMNAS-1; l++)
                {
                    if (!(k==i && l==j))
                    {
                        omega=0.0;
                        for ( m = 0; m < mu; m++)
                        {
                            omega=omega+(patron[i][j][m]-a[m])*(patron[k][l][m]-a[m]);
                        }
                        umbral[i][j]=umbral[i][j]+0.5*omega/(1.0*FILAS*COLUMNAS);
                    }
                    
                }
                
            }
            
        }
        
    }
    
    
    /* En cada iteración consideramos un patrón más. */
    for ( q = 1; q <= mu; q++)
    {
        /* Inicializamos valores a cero */
        memoria=0;

        /* Generamos una configuración aleatoria para la matriz. */
        for ( i = 0; i <= FILAS-1; i++)
        {
            for ( j = 0; j <= COLUMNAS-1; j++)
            {
                matriz[i][j]=gsl_rng_uniform_int(tau,2);
            }
            
        }
         
    /* Hacemos los pasos Montecarlo, en los que se hace el intento de cambio. */
        
        for ( r = 1; r <= pasos+1; r++)
        {
            i=gsl_rng_uniform_int(tau,FILAS);
            j=gsl_rng_uniform_int(tau,COLUMNAS);
            h=0.0;
            for ( k = 0; k <= FILAS-1; k++)
            {
                for ( l = 0; l <= COLUMNAS-1; l++)
                {
                    if (!(k==i && l==j))
                    {
                        omega=0.0;
                        for ( m = 0; m < q; m++)
                        {
                            omega=omega+(patron[i][j][m]-a[m])*(patron[k][l][m]-a[m]);
                        }
                        omega=omega/(1.0*FILAS*COLUMNAS);
                        h=h+omega*matriz[k][l];
                    }
                    
                }
                
            }
            h=(1-2*matriz[i][j])*(umbral[i][j]-h);
            

            p=exp(-h/temp);
            if (p >= 1.0)
            {
                p=1.0;
            }
            
            aux=gsl_rng_uniform(tau);
            if (aux < p)
            {
                matriz[i][j]=1-matriz[i][j]; 
            }
        }

        /* Solapamiento final de la matriz con cada patrón. Escribimos la memoria para cada número de patrones almacenados.*/
        
        for ( m = 0; m < q; m++)
        {
            solapamiento[m]=0.0;
            for ( i = 0; i <= FILAS-1; i++)
            {
                for ( j = 0; j <= COLUMNAS-1; j++)
                {
                    solapamiento[m]=solapamiento[m]+(patron[i][j][m]-a[m])*(matriz[i][j]-a[m]);
                }
            }
            solapamiento[m]=solapamiento[m]/(1.0*FILAS*COLUMNAS*a[m]*(1-a[m]));
            if (solapamiento[m]>=0.75 || solapamiento[m]<=-0.75)
            {
                memoria++;
            }
        }
        memoria_fich << q << '\t' << memoria << '\t' << memoria/(1.0*FILAS*COLUMNAS) << '\n';
    }
    memoria_fich.close();
    return 0;
}