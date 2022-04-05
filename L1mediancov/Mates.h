#ifndef _MATES_H

#define _MATES_H
 
/* Proyectar ortogonalmente el vector columna Vector sobre la base 
vectorial Matriz, y guardar el vector proyección resultado 
en VectorResultado.  */

void Proyectar(double * const Vector,double * const Matriz,
        double * const VectorResultado,
        int Dimension,int NumVectoresBase);

/* Proyectar ortogonalmente el vector columna Vector sobre la base 
vectorial Matriz, y guardar el vector proyección resultado 
en VectorResultado y la expresión del vector proyección en coordenadas
de la base en VectorResultadoEnBase.  */

void ProyectarExtra(double * const Vector,double * const Matriz,
        double * const VectorResultado,
        double * const VectorResultadoEnBase,
        int Dimension,int NumVectoresBase);
        
/* Hallar el vector diferencia entre dos vectores. */
void Diferencia(double * const VectorEntrada1,double * const VectorEntrada2,
    double * const VectorResultado,int Dimension);
    
/* Hallar el cuadrado de la norma euclídea de un vector */
void NormaCuadrado(double * const Vector,double * const Resultado,int Dimension);

/* Hallar los autovalores y autovectores de una matriz real y simétrica, 
clasificados según el orden descendente de autovalores */
void AutoValVec(double * const Matriz,double * const AutoValores,
    double * const AutoVectores,int Dimension);

/* Ortonormalizar una base vectorial */
void Ortonormalizar(double *BaseEntrada,double *BaseOrto,
    int Dimension,int NumVectoresBase);

/* Calcula la inversa de la matriz cuadrada A, utilizando la descomposición LU.
El resultado es algo más exacto que el que produce la función InversaNorma() */
void Inversa(double *A,double *InversaA,int Dimension);

/* Calcula la inversa de la matriz cuadrada A, y también calcula la norma
2 (o norma espectral, o "norma" a secas) de la matriz A, y también la norma de la inversa
de A. Utiliza para todo ello SVD. Podemos hacer que no calcule una norma o que no 
calcule la inversa poniendo el parámetro correspondiente a NULL */
void InversaNorma(double *A,double *InversaA,double *NormaA,
    double *NormaInversaA,int Dimension);

/* Producto de un escalar por una matriz. Soporta que Matriz==Result */
void ProductoEscalarMatriz(double Escalar,double *Matriz,double *Result,
    int NumFilas,int NumColumnas);
        
/* Suma de matrices. Soporta que alguno de los sumandos coincida con el resultado */
void SumaMatricial(double *A,double *B,double *Result,int NumFilas,int NumColumnas);

/* Diferencia de matrices */
void DiferenciaMatricial(double *A,double *B,double *Result,int NumFilas,int NumColumnas);

/* Producto de matrices */
void ProductoMatricial(double *A,double *B,double *Result,int NumFilasA,
    int NumColumnasA,int NumColumnasB);  

/* Halla la diagonal del producto de dos matrices, es decir,
 Result = diag ( A * B ), donde Result es un vector. Para que 
 pueda realizarse esto es preciso que el número de filas de A
 coincida con el número de columnas de B. */
void DiagonalProductoMatricial(double *A,double *B,double *Result,
    int NumFilasA,int NumColumnasA);
    
/* Matriz traspuesta */
void Traspuesta(double *A,double *TraspuestaA,int NumFilasA,int NumColumnasA);  

/* Sumar una matriz diagonal a una matriz cuadrada. Si Result==NULL,
se hace sobre la propia matriz A */
void SumarMatrizDiagonal(double *A,double *MatrizDiagonal,double *Result,int Dimension);

/* Sumar una constante a los elementos de la diagonal de una matriz cuadrada. Si Result==NULL,
se hace sobre la propia matriz A */
void SumarDiagonalConstante(double *A,double Valor,double *Result,int Dimension);

/* Extraer la diagonal principal de una matriz cuadrada */
void ExtraerDiagonal(double *A,double *DiagonalA,int Dimension);


#endif

