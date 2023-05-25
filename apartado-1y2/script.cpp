/* SIMULACIÓN DEL MODELO DE HOPFIELD DE RED NEURONAL USANDO EL MODELO DE ISING CON UN SOLO PATRÓN */
/* EJERCICIO VOLUNTARIO 2 PRIMER Y SEGUNDO APARTADO DE LA LECCIÓN 2 DE FÍSICA COMPTACIONAL */

#include <iostream>
#include <cmath>
#include <fstream>
#include <stdio.h>
#include <cstdlib>
#include "gsl_rng.h"

#define FILAS 100   // 50
#define COLUMNAS 82   // 47

using namespace std;

gsl_rng *tau;

int main (void)
{
    /* Definimos inicialmente algunas variables necesarias. */
    int i,j,k,l,m,n,pasos;
    double a,h,p,temp,aux,deformacion,solapamiento;
    short int matriz[FILAS][COLUMNAS];
    short int patron[FILAS][COLUMNAS];
    double umbral[FILAS][COLUMNAS];
    char digito;
    ofstream matriz_fich;
    ofstream solap_fich;
    ifstream patron_fich;
    

    temp=0.1;

    pasos=200*FILAS*COLUMNAS;

    cout << "¿Configuración inicial aleatoria (0) o patrón deformado(1)?:";
    cin >> n;

    matriz_fich.open("/Users/alzorrcarri/Documents/cphys/leccion2-voluntario/apartado-1y2/datos_hopfield.dat");
    patron_fich.open("/Users/alzorrcarri/Documents/cphys/leccion2-voluntario/apartado-1y2/piplup.txt");
    solap_fich.open("/Users/alzorrcarri/Documents/cphys/leccion2-voluntario/apartado-1y2/solapamiento.dat");


    extern gsl_rng *tau; //Puntero al estado del número aleatorio
    int semilla=87412853; //Semilla del generador de números aleatorios

    tau=gsl_rng_alloc(gsl_rng_taus); //Inicializamos el puntero
    gsl_rng_set(tau,semilla); //Inicializamos la semilla

    /* Copiamos el dibujo original en una matriz (patron) */
    for ( i = 0; i <= FILAS-1; i++)
        {
            for ( j = 0; j <= COLUMNAS-1; j++)
            {
                patron_fich >> digito;
                digito=digito-'0';
                patron[i][j]=digito;
            }
        }
    
    /* Generamos una configuración aleatoria o no según el valor de la variable aleatorio. */
    if (n==0)
    {
        for ( i = 0; i <= FILAS-1; i++)
        {
            for ( j = 0; j <= COLUMNAS-1; j++)
            {
                matriz[i][j]=gsl_rng_uniform_int(tau,2);
            }
        
        }
    } else if (n==1)
    {
        cout << "Porcentaje de deformación inicial:";
        cin >> deformacion;

        for ( i = 0; i <= FILAS-1; i++)
        {
            for ( j = 0; j <= COLUMNAS-1; j++)
            {
                matriz[i][j]=patron[i][j];
            }
        }

        for ( k = 1; k <= (FILAS*COLUMNAS*deformacion/100.0) ; k++)
        {
            i=gsl_rng_uniform_int(tau,FILAS);
            j=gsl_rng_uniform_int(tau,COLUMNAS);
            matriz[i][j]=1-matriz[i][j];
        }
    }
    
    /* Calculamos el parámetro a */
    a=0.0;
    for ( i = 0; i <= FILAS-1; i++)
    {
        for ( j = 0; j <= COLUMNAS-1; j++)
        {
            a=a+patron[i][j];
        }
        
    }
    a=a/(1.0*FILAS*COLUMNAS);

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
                    if (!(i==k && j==l))
                    {
                        umbral[i][j]=umbral[i][j]+0.5*(patron[i][j]-a)*(patron[k][l]-a)/(1.0*FILAS*COLUMNAS);
                    }
                    
                }
                
            }
            
        }
        
    }
    
    /* Copiamos la matriz inicial en el fichero */

    for ( i = 0; i <= FILAS-1; i++)
    {
        for ( j = 0; j <= COLUMNAS-2; j++)
        {
            matriz_fich << matriz[i][j] << ',';
            //cout << matriz[i][j];
        }
        matriz_fich << matriz[i][COLUMNAS-1] << '\n';
        //cout << matriz[i][COLUMNAS-1] << '\n';
        
    }
    matriz_fich << '\n';
    
/* Hacemos los pasos Montecarlo, en los que se hace el intento de cambio. */
    n=1;
    for ( m = 1; m <= pasos+1; m++)
    {
        i=gsl_rng_uniform_int(tau,FILAS);
        j=gsl_rng_uniform_int(tau,COLUMNAS);
        h=0.0;
        for ( k = 0; k <= FILAS-1; k++)
        {
            for ( l = 0; l <= COLUMNAS-1; l++)
            {
                if (!(i==k && j==l))
                {
                    h=h+(patron[i][j]-a)*(patron[k][l]-a)*matriz[k][l]/(1.0*FILAS*COLUMNAS);
                }
                
            }
            
        }
        h=(1-2*matriz[i][j])*(umbral[i][j]-h);

        p=exp(-h/temp);
        if (p > 1.0)
        {
            p=1.0;
        }
        
        aux=gsl_rng_uniform(tau);
        if (aux < p)
        {
            matriz[i][j]=1-matriz[i][j]; 
        }


        /* Escribimos la matriz en el fichero solo en cada paso Montecarlo y el solapamiento */
        if (m%(FILAS*COLUMNAS) == 1)
        {
            solapamiento=0.0;
            for ( i = 0; i <= FILAS-1; i++)
            {
                for ( j = 0; j <= COLUMNAS-2; j++)
                {
                    matriz_fich << matriz[i][j] << ',';
                    //cout << matriz[i][j] << ',';
                    solapamiento=solapamiento+(patron[i][j]-a)*(matriz[i][j]-a);
                }
                matriz_fich << matriz[i][COLUMNAS-1] << '\n';
                solapamiento=solapamiento+(patron[i][COLUMNAS-1]-a)*(matriz[i][COLUMNAS-1]-a);
                //cout << matriz[i][82] << endl; 
            }
            matriz_fich << '\n';
            solapamiento=solapamiento/(FILAS*COLUMNAS*a*(1-a));
            solap_fich << n << '\t' << solapamiento << endl;
            n++;
            
        }
        
    }
    
    matriz_fich.close();
    patron_fich.close();
    solap_fich.close();
    return 0;
}
