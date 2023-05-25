/* SIMULACIÓN DEL MODELO DE HOPFIELD DE RED NEURONAL USANDO EL MODELO DE ISING CON UN SOLO PATRÓN */
/* EJERCICIO VOLUNTARIO 2 TERCER APARTADO DE LA LECCIÓN 2 DE FÍSICA COMPUTACIONAL */

#include <iostream>
#include <cmath>
#include <fstream>
#include <stdio.h>
#include <cstdlib>
#include "gsl_rng.h"

#define FILAS 19
#define COLUMNAS 19

using namespace std;

gsl_rng *tau;

int main (void)
{
    /* Definimos inicialmente algunas variables necesarias. */
    int i,j,k,l,n,m,q,pasos;
    double h,p,temp,aux,deformacion,omega;
    double a[3],solapamiento[3];
    short int matriz[FILAS][COLUMNAS];
    short int patron1[FILAS][COLUMNAS];
    short int patron2[FILAS][COLUMNAS];
    short int patron3[FILAS][COLUMNAS];
    double umbral[FILAS][COLUMNAS];
    char digito;
    ofstream matriz_fich;
    ofstream solap_fich;
    ifstream patron1_fich;
    ifstream patron2_fich;
    ifstream patron3_fich;
    

    temp=0.05;
    pasos=200*FILAS*COLUMNAS;

    

    cout << "¿Configuración inicial aleatoria (0) o patrón deformado(1)?:";
    cin >> n;

    matriz_fich.open("/Users/alzorrcarri/Documents/cphys/leccion2-voluntario/apartado-3/datos_hopfield.dat");
    patron1_fich.open("/Users/alzorrcarri/Documents/cphys/leccion2-voluntario/apartado-3/1.txt");
    patron2_fich.open("/Users/alzorrcarri/Documents/cphys/leccion2-voluntario/apartado-3/2.txt");
    patron3_fich.open("/Users/alzorrcarri/Documents/cphys/leccion2-voluntario/apartado-3/3.txt");
    solap_fich.open("/Users/alzorrcarri/Documents/cphys/leccion2-voluntario/apartado-3/solapamiento.dat");


    extern gsl_rng *tau; //Puntero al estado del número aleatorio
    int semilla=648365902; //Semilla del generador de números aleatorios

    tau=gsl_rng_alloc(gsl_rng_taus); //Inicializamos el puntero
    gsl_rng_set(tau,semilla); //Inicializamos la semilla

    /* Copiamos el dibujo original en una matriz (patron) */
    for ( i = 0; i <= FILAS-1; i++)
        {
            for ( j = 0; j <= COLUMNAS-1; j++)
            {
                patron1_fich >> digito;
                digito=digito-'0';
                patron1[i][j]=digito;

                patron2_fich >> digito;
                digito=digito-'0';
                patron2[i][j]=digito;

                patron3_fich >> digito;
                digito=digito-'0';
                patron3[i][j]=digito;
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
        cout << "¿Patrón a deformar? (1,2 o 3)?:";
        cin >> q;
        cout << "Porcentaje de deformación inicial:";
        cin >> deformacion;

        if (q==1)
        {
            for ( i = 0; i <= FILAS-1; i++)
            {
                for ( j = 0; j <= COLUMNAS-1; j++)
                {
                   matriz[i][j]=patron1[i][j];
                }
            }
        }
        if (q==2)
        {
            for ( i = 0; i <= FILAS-1; i++)
            {
                for ( j = 0; j <= COLUMNAS-1; j++)
                {
                   matriz[i][j]=patron2[i][j];
                }
            }
        }
        if (q==3)
        {
            for ( i = 0; i <= FILAS-1; i++)
            {
                for ( j = 0; j <= COLUMNAS-1; j++)
                {
                   matriz[i][j]=patron3[i][j];
                }
            }
        }

        for ( k = 1; k <= (FILAS*COLUMNAS*deformacion/100.0) ; k++)
        {
            i=gsl_rng_uniform_int(tau,FILAS);
            j=gsl_rng_uniform_int(tau,COLUMNAS);
            matriz[i][j]=1-matriz[i][j];
        }
        
        
    }
    
    /* Calculamos el vector a para cada patrón */
    a[0]=a[1]=a[2]=0.0;
    for ( i = 0; i <= FILAS-1; i++)
    {
        for ( j = 0; j <= COLUMNAS-1; j++)
        {
            a[0]=a[0]+patron1[i][j];
        }
        
    }
    a[0]=a[0]/(FILAS*COLUMNAS);

    for ( i = 0; i <= FILAS-1; i++)
    {
        for ( j = 0; j <= COLUMNAS-1; j++)
        {
            a[1]=a[1]+patron2[i][j];
        }
        
    }
    a[1]=a[1]/(FILAS*COLUMNAS);

    for ( i = 0; i <= FILAS-1; i++)
    {
        for ( j = 0; j <= COLUMNAS-1; j++)
        {
            a[2]=a[2]+patron3[i][j];
        }
        
    }
    a[2]=a[2]/(FILAS*COLUMNAS);

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
                        omega=(patron1[i][j]-a[0])*(patron1[k][l]-a[0])+(patron2[i][j]-a[1])*(patron2[k][l]-a[1])+(patron3[i][j]-a[2])*(patron3[k][l]-a[2]);
                        omega=omega/(1.0*FILAS*COLUMNAS);
                        umbral[i][j]=umbral[i][j]+0.5*omega;
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
                if (!(k==i && l==j))
                {
                    omega=(patron1[i][j]-a[0])*(patron1[k][l]-a[0])+(patron2[i][j]-a[1])*(patron2[k][l]-a[1])+(patron3[i][j]-a[2])*(patron3[k][l]-a[2]);
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


        /* Escribimos la matriz en el fichero solo en cada paso Montecarlo y el solapamiento */
        if (m%(FILAS*COLUMNAS) == 1)
        {
            solapamiento[0]=0.0;
            solapamiento[1]=0.0;
            solapamiento[2]=0.0;
            for ( i = 0; i <= FILAS-1; i++)
            {
                for ( j = 0; j <= COLUMNAS-2; j++)
                {
                    matriz_fich << matriz[i][j] << ',';
                    solapamiento[0]=solapamiento[0]+(patron1[i][j]-a[0])*(matriz[i][j]-a[0]);
                    solapamiento[1]=solapamiento[1]+(patron2[i][j]-a[1])*(matriz[i][j]-a[1]);
                    solapamiento[2]=solapamiento[2]+(patron3[i][j]-a[2])*(matriz[i][j]-a[2]);
                }
                matriz_fich << matriz[i][COLUMNAS-1] << '\n';
                solapamiento[0]=solapamiento[0]+(patron1[i][COLUMNAS-1]-a[0])*(matriz[i][COLUMNAS-1]-a[0]);
                solapamiento[1]=solapamiento[1]+(patron2[i][COLUMNAS-1]-a[1])*(matriz[i][COLUMNAS-1]-a[1]);
                solapamiento[2]=solapamiento[2]+(patron3[i][COLUMNAS-1]-a[2])*(matriz[i][COLUMNAS-1]-a[2]);
            }
            matriz_fich << '\n';
            solapamiento[0]=solapamiento[0]/(FILAS*COLUMNAS*a[0]*(1-a[0]));
            solapamiento[1]=solapamiento[1]/(FILAS*COLUMNAS*a[1]*(1-a[1]));
            solapamiento[2]=solapamiento[2]/(FILAS*COLUMNAS*a[2]*(1-a[2]));
            solap_fich << n << '\t' << solapamiento[0] << '\t' << solapamiento[1] << '\t' << solapamiento[2] << endl;
            n++;
        }
        
    }
    
    matriz_fich.close();
    patron1_fich.close();
    patron2_fich.close();
    patron3_fich.close();
    solap_fich.close();
    return 0;
}
