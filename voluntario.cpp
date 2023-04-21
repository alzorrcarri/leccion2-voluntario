/* SIMULACIÓN DEL MODELO DE HOPFIELD DE RED NEURONAL USANDO EL MODELO DE ISING */
/* EJERCICIO VOLUNTARIO 2 DE LA LECCIÓN 2 DE FÍSICA COMPTACIONAL */

#include <iostream>
#include <cmath>
#include <fstream>
#include <stdio.h>
#include "gsl_rng.h"

using namespace std;

gsl_rng *tau;

int main (void)
{
    /* Definimos inicialmente algunas variables necesarias.*/
    int i,j,k,l,m;
    int aleatorio;
    double a,h,p,temp,aux,pasos;
    double deformacion,solapamiento;
    int matriz[100][82];
    int patron[100][82];
    double omega[100][82][100][82];
    double umbral[100][82];
    char caracter;
    ofstream matriz_fich;
    ofstream solap_fich;
    ifstream patron_fich;
    
    /* Pedimos los valores iniciales de la temperatura del sistema y el número de pasos Montecarlo.*/
    
    //cout << "Introducir la temperatura de la red neuronal en Kelvin:";
    //cin >> temp;
    temp=0.0001;

    //cout << "Introducir el número de pasos Montecarlo:";
    //cin >> pasos;
    pasos=200*100*82;

    cout << "¿Configuración inicial aleatoria (0) o patrón deformado(1)?";
    cin >> aleatorio;

    matriz_fich.open("/Users/alzorrcarri/Documents/cphys/leccion2-voluntario/datos_hopfield.dat");
    patron_fich.open("/Users/alzorrcarri/Documents/cphys/leccion2-voluntario/piplup.txt");
    solap_fich.open("/Users/alzorrcarri/Documents/cphys/leccion2-voluntario/solapamiento.dat");


    extern gsl_rng *tau; //Puntero al estado del número aleatorio
    int semilla=87412853; //Semilla del generador de números aleatorios

    tau=gsl_rng_alloc(gsl_rng_taus); //Inicializamos el puntero
    gsl_rng_set(tau,semilla); //Inicializamos la semilla

    /* Copiamos el dibujo original en una matriz (patron) */
    for ( i = 0; i <= 99; i++)
        {
            for ( j = 0; j <= 81; j++)
            {
                patron_fich >> caracter;
                caracter=caracter-'0';
                patron[i][j]=caracter;
            }
        }
    
    /* Generamos una configuración aleatoria o no según el valor de la variable aleatorio. */
    if (aleatorio==0)
    {
        for ( i = 0; i <= 99; i++)
        {
            for ( j = 0; j <= 81; j++)
            {
                matriz[i][j]=1.0*gsl_rng_uniform_int(tau,2);
            }
        
        }
    } else if (aleatorio==1)
    {
        cout << "Porcentaje de deformación inicial:";
        cin >> deformacion;

        for ( i = 0; i <= 99; i++)
        {
            for ( j = 0; j <= 81; j++)
            {
                matriz[i][j]=patron[i][j];
            }
        }

        for ( k = 1; k <= (100*82*deformacion/100.0) ; k++)
        {
            i=gsl_rng_uniform_int(tau,100);
            j=gsl_rng_uniform_int(tau,82);
            matriz[i][j]=1-matriz[i][j];
        }
    }
    
    /* Calculamos el parámetro a */
    a=0.0;
    for ( i = 0; i <= 99; i++)
    {
        for ( j = 0; j <= 81; j++)
        {
            a=a+patron[i][j];
        }
        
    }
    a=a/(100.0*82);

    /* Calculamos la matriz que describe la interacción entre neuronas (omega) */

    for ( i = 0; i <= 99; i++)
    {
        for ( j = 0; i <= 81; j++)
        {
            for ( k = 0; k <= 99; k++)
            {
                for ( l = 0; l <= 81; l++)
                {
                    if (i == k && j == l)
                    {
                        omega[i][j][k][l]=0.0;
                    }
                    else
                    {
                        omega[i][j][k][l]=(patron[i][j]-a)*(patron[k][l]-a)/(100.0*82);
                    }
                    
                }
                
            }
            
        }
        
    }

    /* Calculamos la matriz umbral de disparo (umbral) */
    for ( i = 0; i <= 99; i++)
    {
        for ( j = 0; j <= 81; j++)
        {
            umbral[i][j]=0.0;
            for ( k = 0; k <= 99; k++)
            {
                for ( l = 0; l <= 81; l++)
                {
                    umbral[i][j]=umbral[i][j]+0.5*omega[i][j][k][l];
                }
                
            }
            
        }
        
    }
    
    /* Copiamos la matriz inicial en el fichero */

    for ( i = 0; i <= 99; i++)
    {
        for ( j = 0; j <= 80; j++)
        {
            matriz_fich << matriz[i][j];
            //cout << matriz[i][j];
        }
        matriz_fich << matriz[i][81] << '\n';
        //cout << matriz[i][81] << '\n';
        
    }
    matriz_fich << '\n';
    
/* Hacemos los pasos Montecarlo, en los que se hace el intento de cambio de los espines */
    for ( m = 1; m <= pasos+1; m++)
    {
        i=gsl_rng_uniform_int(tau,100);
        j=gsl_rng_uniform_int(tau,82);
        h=0.0;
        for ( k = 0; k <= 99; k++)
        {
            for ( l = 0; l <= 81; l++)
            {
                h=h+omega[i][j][k][l]*matriz[k][l];
            }
            
        }
        h=(1-2*matriz[i][j])*(umbral[i][j]-omega[i][j][i][j]*(0.5+matriz[i][j])-h);
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
        if (m%(100*82) == 1)
        {
            solapamiento=0.0;
            for ( i = 0; i <= 99; i++)
            {
                for ( j = 0; j <= 80; j++)
                {
                    matriz_fich << matriz[i][j] << ',';
                    //cout << matriz[i][j] << ',';
                    solapamiento=solapamiento+(patron[i][j]-a)*(matriz[i][j]-a);
                }
                matriz_fich << matriz[i][81] << endl;
                solapamiento=solapamiento+(patron[i][81]-a)*(matriz[i][81]-a);
                //cout << matriz[i][82] << endl; 
            }
            matriz_fich << endl;
            solapamiento=solapamiento/(100.0*82*a*(1-a));
            solap_fich << solapamiento << endl;
            
        }
        
    }
    
    matriz_fich.close();
    patron_fich.close();
    solap_fich.close();
    return 0;
}