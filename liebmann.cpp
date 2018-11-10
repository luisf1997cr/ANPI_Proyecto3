/*
 * liebmann.cpp
 *
 *  
 */
#include <stdio.h>
#include "iostream"
#include <iomanip>
#include <cmath>
#include <math.h>
#include <armadillo>
#include <fstream>
#include <omp.h>
#include <ctime>
using namespace std;
using namespace arma;

/*! Parámetros: 
 * I: fila en la que se encuentra un nodo
 * J: Columna en la que se encuentra un nodo
 * fila: fila en la que se escriben datos del vector de soluciones
 * tam:tamaño de la placa
 * arriba,abajo,izquierda,derecha: temperaturas en cada uno de los lados de la placa  */
/*!
 * La función de este método es aplicar la ecuación de liebmann al nodo de la placa ubicado
 * en i , j y para encontrar una ecuacion para dicho nodo donde los datos de esta ecuacion se
 * almacenan en una matriz de tamaño tamxtam
 */

mat temperaturaPunto(int I, int J, int fila, int tam, colvec *colum, double arriba, double abajo, double izquierda,
                     double derecha)
{
    int i = I;
    int j = J;
    mat mat(tam, tam);
    mat.zeros();
    mat.at(i, j) = -4;

    if (i - 1 == 0)
    {

        colum->at(fila, 0) = colum->at(fila, 0) + izquierda;
    }
    else
    {
        mat.at(i - 1, j) = 1;
    }
    if (i + 1 == (tam))
    {
        colum->at(fila, 0) = colum->at(fila, 0) + derecha;
    }
    else
    {
        mat.at(i + 1, j) = 1;
    }
    if (j - 1 == 0)
    {
        colum->at(fila, 0) = colum->at(fila, 0) + abajo;
    }

    else
    {
        mat.at(i, j - 1) = 1;
    }
    if (j + 1 == tam)
    {
        colum->at(fila, 0) = colum->at(fila, 0) + arriba;
    }
    else
    {
        mat.at(i, j + 1) = 1;
    }
    return mat;
}

/*!
 * Parámetros: parcial: matriz de datos de una ecuación calculada en temperaturaPunto()
 * final: matriz que representa un sistema de ecuaiones
 * tam:tamaño de la placa
 * fila: fila de la placa a la cual corresponde el nodo de la matriz parcial
 */
/*!
 * este método recibe un matriz que representa las constanttes de una ecuacion
 * y les agrega una matriz que representa un sistema de ecuaciones
 */
void matrizToFilaToMatriz(mat parcial, mat *final, int tam, int fila)
{

    int z = 0;
    for (int i = 1; i < tam; i++)
    {
        for (int j = 1; j < tam; j++)
        {
            final->at(fila, z) = parcial.at(i, j);
            z++;
        }
    }
}

/*!
 * Parámetros: tam: tamaño de la placa matriz: matriz que contiene las constantes de un sistema de ecuaciones
 *
 */
/*!
 * Resuelve el sistema de ecuaciones que representa matriz usando el método de
 * Gauss-Seidel y retorna un vector con las soluciones obtenidas
 */
mat liebmannSolver(int tam, mat matriz)
{

    double error = 0.001, y; //error= máximo error permitido
    int numEcua, flag = 0;
    numEcua = ((tam - 1) * (tam - 1));
    mat matrizFinal(numEcua, numEcua + 1); // matriz de ecuacion
    matrizFinal = matriz;
    mat x(1, numEcua); //matriz con la ecuaciones
    x.zeros();

    int i = 0;
    int j = 0;
    for (i = 0; i < numEcua; i++)
    { //matriz diagonalmente dominante
        for (int k = i + 1; k < numEcua; k++)
        {
            if (abs(matrizFinal.at(i, i)) < abs(matrizFinal.at(k, i)))
            {
#pragma omp parallel for shared(matrizFinal) num_threads(2)
                for (j = 0; j <= numEcua; j++)
                {

                    int temp = matrizFinal.at(i, j);
                    matrizFinal.at(i, j) = matrizFinal.at(k, j);
                    matrizFinal.at(k, j) = temp;
                }
            }
        }
    }

    do //cálculo de temperaturas (gauss seidel)
    {
#pragma omp parallel for shared(matrizFinal, x, y) num_threads(2)
        for (int i = 0; i < numEcua; i++)
        {
            y = x.at(0, i);
            x.at(0, i) = matrizFinal.at(i, numEcua);
            for (int j = 0; j < numEcua; j++)
            {
                if (j != i)
                    x.at(0, i) = x.at(0, i) - matrizFinal.at(i, j) * x.at(j);
            }
            x.at(0, i) = x.at(0, i) / matrizFinal.at(i, i);
            if (abs(x.at(0, i) - y) <= error)
            { //Comparación con el valor anterior
                flag++;
            }
        }

    } while (flag < numEcua);

    /*cout<<"\n vector solución:\n";
    for (int i=0;i<numEcua;i++){
        temp=x.at(0,i);
        if((temp<=0 && temp>-0.1)||(temp>0 && temp<0.1)){
                    x.at(0,i)=0;
                }
        //cout<<"x"<<i<<" = "<<x.at(0,i)<<endl;
            }*/
    return x;
}

mat liebmannSolverSinMejora(int tam, mat matriz)
{

    double error = 0.0001, y; //error= máximo error permitido
    int numEcua, flag = 0;
    numEcua = ((tam - 1) * (tam - 1));
    mat matrizFinal(numEcua, numEcua + 1); // matriz de ecuacion
    matrizFinal = matriz;
    mat x(1, numEcua); //matriz con la ecuaciones
    x.zeros();
    for (int i = 0; i < numEcua; i++) //matriz diagonalmente dominante

        for (int k = i + 1; k < numEcua; k++)
            if (abs(matrizFinal.at(i, i)) < abs(matrizFinal.at(k, i)))

                for (int j = 0; j <= numEcua; j++)
                {
                    int temp = matrizFinal.at(i, j);
                    matrizFinal.at(i, j) = matrizFinal.at(k, j);
                    matrizFinal.at(k, j) = temp;
                }

    do //cálculo de temperaturas
    {

        for (int i = 0; i < numEcua; i++)
        {
            y = x.at(0, i);
            x.at(0, i) = matrizFinal.at(i, numEcua);
            for (int j = 0; j < numEcua; j++)
            {
                if (j != i)
                    x.at(0, i) = x.at(0, i) - matrizFinal.at(i, j) * x.at(j);
            }
            x.at(0, i) = x.at(0, i) / matrizFinal.at(i, i);
            if (abs(x.at(0, i) - y) <= error)
            {
                flag++;
            } //Comparación con el valor anterior
        }

    } while (flag < numEcua);

    return x;
}

/*!
 * Parámetros: arriba,abajo,izquierda,derecha: temperaturas de los lados de la placa
 * tam tamaño de la placa
 */
/*!
 * encuentra las ecuaciones de temperatura para todos los nodos de la placa y
 * completa el sistema de ecuaciones en una matriz
 */
/**mat temperatura(int arriba,int abajo, int izquierda,int derecha,int tam){
    mat matriz((tam-1)*(tam-1),(tam-1)*(tam-1));
    matriz.zeros();
    colvec column((tam-1)*(tam-1),1);
    column.zeros();

    int z=0;
    int i,j;
    
    for(i=1;i<=tam-1;i++){ 
            for(j=1;j<=tam-1;j++){
                
                matrizToFilaToMatriz(temperaturaPunto(j,i,z,tam,&column,arriba,abajo,izquierda,derecha),&matriz,tam,z);
                
                z++;

            }
           
        }
  

     matriz.insert_cols(((tam-1)*(tam-1)),-1*column);
     




    return matriz;
}

---------No me funciona tengo que revisar bien como analizo la matriz--------------------------------------------
**/
/*!
 * Parámetros: resp: vector con las temperaturas obtenidas para la placa
 * tam: tamaño de la placa
 */
/*!
 * convierte de un vector a una matriz con la posición correcta de cada temperatura en la placa
 */
mat rowToMatrix(mat resp, int tam)
{
    mat matriz(tam - 1, tam - 1);
    int z = 0;
    int i, j;
#pragma omp for collapse(2)
    for (j = 0; j < tam - 1; j++)
    {
        for (i = tam - 2; i >= 0; i--)
        {

            matriz.at(i, j) = resp.at(0, z);
            z++;
        }
    }
    return matriz;
}
mat rowToMatrixSinMejora(mat resp, int tam)
{
    mat matriz(tam - 1, tam - 1);
    int z = 0;
    int i, j;

    for (j = 0; j < tam - 1; j++)
    {
        for (i = tam - 2; i >= 0; i--)
        {

            matriz.at(i, j) = resp.at(0, z);
            z++;
        }
    }
    return matriz;
}
/*!
 * Parámetros: arriba,abajo,izquierda,derecha: temperaturas de los lados de la placa
 * temps: matriz de temperaturas
 */
/*!
 * calcula un vector para los nodos de la placa
 */
mat vectoresFlujo(mat temps, int tam, int arriba, int abajo, int izquierda, int derecha)
{
    mat matriz(1, (tam - 1) * (tam - 1) * 2);
    double qx = 0, qy = 0, k = 0.49, valxA = 0, valxB = 0, valyA = 0, valyB = 0;
    int z = 0, dX = 10, dY = 10;

    for (int i = 0; i < tam - 1; i++)
    {
        for (int j = 0; j < tam - 1; j++)
        {

            if (j - 1 < 0)
            {
                valxB = izquierda;
                valxA = temps.at(i, j + 1);
            }
            else if (j + 1 >= (tam - 1))
            {
                valxA = derecha;
                valxB = temps.at(i, j - 1);
            }
            else
            {

                valxA = temps.at(i, j + 1);

                valxB = temps.at(i, j - 1);
            }
            qx = -k * ((valxA - valxB) / (2 * dX));

            if (i + 1 >= (tam - 1))
            {
                valyB = abajo;
                valyA = temps.at(i - 1, j);
            }
            else if (i - 1 < 0)
            {
                valyA = arriba;
                valyB = temps.at(i + 1, j);
            }
            else
            {
                valyA = temps.at(i - 1, j);
                valyB = temps.at(i + 1, j);
            }

            qy = -k * ((valyA - valyB) / (2 * dY));

            double magni = sqrt((qx * qx) + (qy * qy));
            double unitX = 0.5 * (qx / magni);
            double unitY = 0.5 * (qy / magni);

            /*if((unitX<0 && unitX>-0.1)||(unitX>0 && unitX<0.1)){
                    unitX=0;
                }
                if((unitY<0 && unitY>-0.1)||(unitY>0 && unitY<0.1)){
                    unitY=0;
                }*/
            matriz.at(0, z) = unitX;
            matriz.at(0, z + 1) = unitY;

            z = z + 2;
        }
    }

    return matriz;
}
