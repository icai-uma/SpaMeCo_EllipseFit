#include "Mates.h"
#include <string.h>
#include <float.h>
#include <stddef.h>
#include "mex.h"
#include "eig.h"
#include "svd.h"
#include "lu.h"
#include "nrutil.h"

 
/* Proyectar ortogonalmente el vector columna Vector sobre la base 
vectorial Matriz, y guardar el vector proyección resultado 
en VectorResultado. */

void Proyectar(double * const Vector,double * const Matriz,
        double * const VectorResultado,int Dimension,int NumVectoresBase)
{
    register double *Limite;
    register double *MiComponente;
    register double *MiElemento;
    int ndx;
    register double Resultado;
    
    memset(VectorResultado,0,sizeof(double)*Dimension);    
    MiElemento=Matriz;
    for(ndx=0;ndx<NumVectoresBase;ndx++)
    {
        /* Hallar el producto escalar del vector de entrada y
        este vector base */
        MiComponente=Vector;
        Limite=MiElemento+Dimension;
        Resultado=0.0;
        while (MiElemento<Limite)
        {
            Resultado+=(*MiElemento)*(*MiComponente);
            MiElemento++;
            MiComponente++;
        }    
        /* Hallar la contribución al vector proyección de este
        vector base*/
        MiComponente=VectorResultado;
        MiElemento-=Dimension;
        while (MiElemento<Limite)
        {
            (*MiComponente)+=Resultado*(*MiElemento);
            MiComponente++;
            MiElemento++;
        }
    }          
}
  
    
/* Proyectar ortogonalmente el vector columna Vector sobre la base 
vectorial Matriz, y guardar el vector proyección resultado 
en VectorResultado y la expresión del vector proyección en coordenadas
de la base en VectorResultadoEnBase. */

void ProyectarExtra(double * const Vector,double * const Matriz,
        double * const VectorResultado,
        double * const VectorResultadoEnBase,
        int Dimension,int NumVectoresBase)
{
    register double *Limite;
    register double *MiComponente;
    register double *MiElemento;
    register double *MiResultadoEnBase;
    int ndx;
    register double Resultado;
    

    memset(VectorResultado,0,sizeof(double)*Dimension);
    MiElemento=Matriz;
    MiResultadoEnBase=VectorResultadoEnBase;
    for(ndx=0;ndx<NumVectoresBase;ndx++)
    {
        /* Hallar el producto escalar del vector de entrada y
        este vector base */
        MiComponente=Vector;
        Limite=MiElemento+Dimension;
        Resultado=0.0;
        while (MiElemento<Limite)
        {
            Resultado+=(*MiElemento)*(*MiComponente);
            MiElemento++;
            MiComponente++;
        }    
        (*MiResultadoEnBase)=Resultado;
        MiResultadoEnBase++;
        /* Hallar la contribución al vector proyección de este
        vector base*/
        MiComponente=VectorResultado;
        MiElemento-=Dimension;
        while (MiElemento<Limite)
        {
            (*MiComponente)+=Resultado*(*MiElemento);
            MiComponente++;
            MiElemento++;
        }
    }           
}  

/* Hallar el vector diferencia entre dos vectores. */
void Diferencia(double * const VectorEntrada1,double * const VectorEntrada2,
    double * const VectorResultado,int Dimension)
{
    register double *MiComponenteEntrada1;
    register double *MiComponenteEntrada2;
    register double *MiComponenteResultado;
    register int ndx;
    
    MiComponenteEntrada1=VectorEntrada1;
    MiComponenteEntrada2=VectorEntrada2;
    MiComponenteResultado=VectorResultado;
    for (ndx=0;ndx<Dimension;ndx++)
    {
        (*MiComponenteResultado)=(*MiComponenteEntrada1)-
                (*MiComponenteEntrada2);
        MiComponenteEntrada1++;
        MiComponenteEntrada2++;
        MiComponenteResultado++;                
    }      
}    
    
/* Hallar el cuadrado de la norma euclídea de un vector */
void NormaCuadrado(double * const Vector,double * const Resultado,int Dimension)
{
    register double *MiComponente;
    register int ndx;
    register double MiResultado;    
    
    MiComponente=Vector;
    MiResultado=0.0;
    for (ndx=0;ndx<Dimension;ndx++)
    {
        MiResultado+=(*MiComponente)*(*MiComponente);
        MiComponente++;                
    }    
    (*Resultado)=MiResultado;  
}
    
/* Hallar los autovalores y autovectores de una matriz real y simétrica, 
clasificados según el orden descendente de autovalores */
void AutoValVec(double * const Matriz,double * const AutoValores,
    double * const AutoVectores,int Dimension)
{
    int i,j;
    double *ptrElemento; /* Para recorrer Matriz */
    double **a; /* Matriz de entrada */
    double **v; /* Matriz de autovectores */
    double *e; /* Diagonal secundaria */
    int nrot; /* Número de iteraciones para el método de Jacobi */
    
    
    /* Reservar memoria para las matrices auxiliares comunes a ambos
    métodos */
    a=dmatrix(1,Dimension,1,Dimension);   
      
    /* Copiar los datos de entrada en a */
    ptrElemento=Matriz;
    for (j=1;j<=Dimension;j++)
    {
        for (i=1;i<=Dimension;i++)
        {
            a[i][j]=(*ptrElemento);
            ptrElemento++;
        }
    }        
    
 
    /* El método de Jacobi se recomienda para matrices de tamaño menor 
    o igual que 10 en los Numerical Recipes */
    if (Dimension<=10)
    {
        /* Método de Jacobi */
        /* Reservar memoria específica de este método */
        v=dmatrix(1,Dimension,1,Dimension);
        /* Poner a cero la matriz v */
        memset(v[1]+1,0,sizeof(double)*Dimension*Dimension);
        /* Poner a cero el vector AutoValores */
        memset(AutoValores,0,sizeof(double)*Dimension);        
         
        /* Ejecutar el método de Jacobi */
        jacobi (a, Dimension, AutoValores-1, v, &nrot);
          
        /* Ordenar los resultados */
        eigsrt (AutoValores-1, v, Dimension); 
        /* Copiar los autovectores de salida en AutoVectores */
        ptrElemento=AutoVectores;
        for (j=1;j<=Dimension;j++)
        {
            for (i=1;i<=Dimension;i++)
            {
                (*ptrElemento)=v[i][j];
                ptrElemento++;                
            }
        }       
        /* Liberar memoria específica de este método */
        free_dmatrix(v,1,Dimension,1,Dimension);
    } 
    else
    {
        /* Método Householder + QL */
        /* Reservar memoria específica de este método */
        e=dvector(1,Dimension);
        /* Poner a cero el vector e */
        memset(e+1,0,sizeof(double)*Dimension);
        /* Poner a cero el vector AutoValores */
        memset(AutoValores,0,sizeof(double)*Dimension);
        
        /* Note que tred2 destruye lo que hubiera en a2copia */
        tred2(a, Dimension, AutoValores-1, e);   
        tqli(AutoValores-1, e, Dimension, a);
        eigsrt (AutoValores-1, a, Dimension);
        /* Copiar los autovectores de salida en AutoVectores */
        ptrElemento=AutoVectores;
        for (j=1;j<=Dimension;j++)
        {
            for (i=1;i<=Dimension;i++)
            {
                (*ptrElemento)=a[i][j];
                ptrElemento++;                
            }
        }   
        /* Liberar memoria específica de este método */
        free_dvector(e,1,Dimension);   
    }    
    
    /* Liberar memoria común a ambos métodos */
    free_dmatrix(a,1,Dimension,1,Dimension);         
}

/* Ortonormalizar una base vectorial */
void Ortonormalizar(double *BaseEntrada,double *BaseOrto,
    int Dimension,int NumVectoresBase)
{
    int i,j;
  
  double **a;
  double **v;
  double *w; 
  double *ptrElemento;

  a=dmatrix(1,Dimension,1,NumVectoresBase);
  v=dmatrix(1,NumVectoresBase,1,NumVectoresBase);
  w=dvector(1,NumVectoresBase);

  /* Poner a cero la matriz v */
  memset(v[1]+1,0,sizeof(double)*NumVectoresBase*NumVectoresBase);
  /* Poner a cero el vector w */
  memset(w+1,0,sizeof(double)*NumVectoresBase);
  
  
  /* Copiar los datos de entrada en a */
  ptrElemento=BaseEntrada;
  for (j=1;j<=NumVectoresBase;j++)
  {
      for (i=1;i<=Dimension;i++)
      {
          a[i][j]=(*ptrElemento);
          ptrElemento++;
      }
  }      
  
  
  /* Cálculo del SVD */ 
  svdcmp(a,Dimension,NumVectoresBase,w,v);

  /* Copiar los vectores de salida en BaseOrto */
  ptrElemento=BaseOrto;
  for (j=1;j<=NumVectoresBase;j++)
  {
      for (i=1;i<=Dimension;i++)
      {
          (*ptrElemento)=a[i][j];          
          ptrElemento++;                
      }
  }

  /* Liberar memoria reservada */  
  free_dmatrix(a,1,Dimension,1,NumVectoresBase);
  free_dmatrix(v,1,NumVectoresBase,1,NumVectoresBase);
  free_dvector(w,1,NumVectoresBase);

}    

/* Calcula la inversa de la matriz cuadrada A, utilizando la descomposición LU.
El resultado es algo más exacto que el que produce la función InversaNorma() */
void Inversa(double *A,double *InversaA,int Dimension)
{
    double **a;
    register double *ptrElemento;
    
    int i,j,*indx;
    double d;
    double *col;

    /* Reservar memoria */
    a=dmatrix(1,Dimension,1,Dimension);    
    
    /* Copiar los datos de entrada en a */
    ptrElemento=A;
    for (j=1;j<=Dimension;j++)
    {
        for (i=1;i<=Dimension;i++)
        {
            a[i][j]=(*ptrElemento);
            ptrElemento++;
        }
    }   
    
    col=dvector(1,Dimension);
    indx=ivector(1,Dimension);
    
    ludcmp(a,Dimension,indx,&d); /*Decompose the matrix just once.*/
    
    for(j=1;j<=Dimension;j++) 
    {   
        /*Find inverse by columns.*/
        for(i=1;i<=Dimension;i++) 
            col[i]=0.0;
        col[j]=1.0;
        lubksb(a,Dimension,indx,col);
        memcpy(InversaA+Dimension*(j-1),col+1,Dimension*sizeof(double));
    }
    
    free_dvector(col,1,Dimension);
    free_ivector(indx,1,Dimension);
    free_dmatrix(a,1,Dimension,1,Dimension);
    
}    

/* Calcula la inversa de la matriz cuadrada A, y también calcula la norma
2 (o norma espectral, o "norma" a secas) de la matriz A, y también la norma de la inversa
de A. Utiliza para todo ello SVD. Podemos hacer que no calcule una norma o que no 
calcule la inversa poniendo el parámetro correspondiente a NULL */
void InversaNorma(double *A,double *InversaA,double *NormaA,double *NormaInversaA,int Dimension)
{
  double **a;
  double **v;
  double *w; 
  register double *ptrElemento;
  register double Suma,Auxiliar;
  register int i,j,k;

  a=dmatrix(1,Dimension,1,Dimension);
  v=dmatrix(1,Dimension,1,Dimension);
  w=dvector(1,Dimension);

  /* Poner a cero la matriz v */
  memset(v[1]+1,0,sizeof(double)*Dimension*Dimension);
  /* Poner a cero el vector w */
  memset(w+1,0,sizeof(double)*Dimension);
  
  
  /* Copiar los datos de entrada en a */
  ptrElemento=A;
  for (j=1;j<=Dimension;j++)
  {
      for (i=1;i<=Dimension;i++)
      {
          a[i][j]=(*ptrElemento);
          ptrElemento++;
      }
  }      
  
  
  /* Cálculo del SVD */ 
  svdcmp(a,Dimension,Dimension,w,v);

  /* La norma de A se calcula como el mayor de los valores singulares
  (ver help norm en MATLAB) */
  if (NormaA!=NULL)
  {
      (*NormaA)=w[1];
      for(i=2;i<=Dimension;i++)
      {
          if (w[i]>(*NormaA))
          {
              (*NormaA)=w[i];
          }    
      } 
  }    
  
  /* La norma de la inversa de A se calcula como el inverso del menor
  de los valores singulares de A (ver Numerical Recipes, apartado SVD) */
  if (NormaInversaA!=NULL)
  {
      Auxiliar=w[1];
      for(i=2;i<=Dimension;i++)
      {
          if (w[i]<Auxiliar)
          {
              Auxiliar=w[i];
          }    
      } 
      (*NormaInversaA)=1.0/Auxiliar;
  } 
  
  /* Dada la descomposición SVD, A=U*W*V', la matriz inversa viene dada
  por inv(A)=V*inv(W)*U', donde inv(W) se calcula invirtiendo cada uno de
  los elementos del propio W. Esto se debe a que tanto V como U son matrices
  ortogonales, y por tanto inv(V)=V', inv(U)=U'. Obtenido de 
  http://kwon3d.com/theory/jkinem/svd.html */
  if (InversaA!=NULL)
  {
      /* Calcular V*inv(W) sobre el propio v */
      for (j=1;j<=Dimension;j++)
      {
          for (i=1;i<=Dimension;i++)
          {
              v[i][j]/=w[j];
          }
      }   
      
      /* Calcular InversaA=(V*inv(W))*U' */
      ptrElemento=InversaA;
      for (j=1;j<=Dimension;j++)
      {
          for (i=1;i<=Dimension;i++)
          {
              Suma=0.0;
              for(k=1;k<=Dimension;k++)
              {
                  /* Nótese que en a está U, y no U' */
                  Suma+=v[i][k]*a[j][k];
              }    
              (*ptrElemento)=Suma;
              ptrElemento++;
          }
      }  
  }     
  
  /* Liberar memoria reservada */  
  free_dmatrix(a,1,Dimension,1,Dimension);
  free_dmatrix(v,1,Dimension,1,Dimension);
  free_dvector(w,1,Dimension);
}

/* Producto de un escalar por una matriz. Soporta que Matriz==Result */
void ProductoEscalarMatriz(double Escalar,double *Matriz,double *Result,
    int NumFilas,int NumColumnas)
{
    register double Factor;
    register double *ptr;
    register double *ptrres;
    register int ndx;
    register int NumElementos;
    
    ptrres=Result;
    ptr=Matriz;
    Factor=Escalar;
    NumElementos=NumFilas*NumColumnas;
    for(ndx=0;ndx<NumElementos;ndx++)
    {
        (*ptrres)=Factor*(*ptr);
        ptrres++;
        ptr++;
    }    
    
}    
/* Suma de matrices. Soporta que alguno de los sumandos coincida con el resultado */
void SumaMatricial(double *A,double *B,double *Result,int NumFilas,int NumColumnas)
{
    register double *ptra;
    register double *ptrb;
    register double *ptrres;
    register int ndx;
    register int NumElementos;
    
    ptra=A;
    ptrb=B;
    ptrres=Result;
    NumElementos=NumFilas*NumColumnas;
    for(ndx=0;ndx<NumElementos;ndx++)
    {
        (*ptrres)=(*ptra)+(*ptrb);
        ptrres++;
        ptra++;
        ptrb++;
    }    
}

/* Diferencia de matrices */
void DiferenciaMatricial(double *A,double *B,double *Result,int NumFilas,int NumColumnas)
{
    register double *ptra;
    register double *ptrb;
    register double *ptrres;
    register int ndx;
    register int NumElementos;
    
    ptra=A;
    ptrb=B;
    ptrres=Result;
    NumElementos=NumFilas*NumColumnas;
    for(ndx=0;ndx<NumElementos;ndx++)
    {
        (*ptrres)=(*ptra)-(*ptrb);
        ptrres++;
        ptra++;
        ptrb++;
    }    
}

/* Producto de matrices */
void ProductoMatricial(double *A,double *B,double *Result,int NumFilasA,
    int NumColumnasA,int NumColumnasB)
{
    register double *ptra;
    register double *ptrb;
    register double *ptrres;
    register int i;
    register int j;
    register int k;
    register double Suma;
    
    ptrres=Result;
    for(j=0;j<NumColumnasB;j++)
    {
        for(i=0;i<NumFilasA;i++)
        {
            Suma=0.0;
            ptrb=B+NumColumnasA*j;
            ptra=A+i;
            for(k=0;k<NumColumnasA;k++)
            {
                Suma+=(*ptra)*(*ptrb);
                ptra+=NumFilasA;
                ptrb++;
            }    
            (*ptrres)=Suma;
            ptrres++;
        }
    }            
}   

/* Halla la diagonal del producto de dos matrices, es decir,
 Result = diag ( A * B ), donde Result es un vector. Para que 
 pueda realizarse esto es preciso que el número de filas de A
 coincida con el número de columnas de B. */
void DiagonalProductoMatricial(double *A,double *B,double *Result,
    int NumFilasA,int NumColumnasA)
{
    register double *ptra;
    register double *ptrb;
    register double *ptrres;
    register int i;
    register int k;
    register double Suma;
    
    ptrres=Result;
    for(i=0;i<NumFilasA;i++)
    {
        Suma=0.0;
        ptrb=B+NumColumnasA*i;
        ptra=A+i;
        for(k=0;k<NumColumnasA;k++)
        {
            Suma+=(*ptra)*(*ptrb);
            ptra+=NumFilasA;
            ptrb++;
        }    
        (*ptrres)=Suma;
        ptrres++;
    }
         
}   

/* Matriz traspuesta */
void Traspuesta(double *A,double *TraspuestaA,int NumFilasA,int NumColumnasA)
{
    register int NdxFila;
    register int NdxCol;
    register double *ptrA;
    
    ptrA=A;
    for(NdxCol=0;NdxCol<NumColumnasA;NdxCol++)
    {
        for(NdxFila=0;NdxFila<NumFilasA;NdxFila++)
        {
            (*(TraspuestaA+NdxFila*NumColumnasA+NdxCol))=(*ptrA);
            ptrA++;
        }
    }        
}    

/* Sumar una matriz diagonal a una matriz cuadrada. Si Result==NULL,
se hace sobre la propia matriz A */
void SumarMatrizDiagonal(double *A,double *MatrizDiagonal,double *Result,int Dimension)
{
    register int NdxElemento;
    register double *ptrDiagonal;
    register double *ptrResult;
    
    /* Copiar la matriz A en el resultado, en su caso*/
    if (Result!=NULL)
    {
        memcpy(Result,A,sizeof(double)*Dimension*Dimension);
    }
    else
    {
        Result=A;
    }        
    
    /* Sumar la matriz diagonal al resultado*/
    ptrDiagonal=MatrizDiagonal;
    ptrResult=Result;
    for(NdxElemento=0;NdxElemento<Dimension;NdxElemento++)
    {
        (*ptrResult)+=(*ptrDiagonal);
        ptrResult+=(Dimension+1);
        ptrDiagonal++;
    }  
}   

/* Sumar una constante a los elementos de la diagonal de una matriz cuadrada. Si Result==NULL,
se hace sobre la propia matriz A */
void SumarDiagonalConstante(double *A,double Valor,double *Result,int Dimension)
{
    register int NdxElemento;
    register double *ptrResult;
    
   /* Copiar la matriz A en el resultado, en su caso*/
    if (Result!=NULL)
    {
        memcpy(Result,A,sizeof(double)*Dimension*Dimension);
    }
    else
    {
        Result=A;
    }
    
    /* Sumar la constante a la diagonal del resultado*/
    ptrResult=Result;
    for(NdxElemento=0;NdxElemento<Dimension;NdxElemento++)
    {
        (*ptrResult)+=Valor;
        ptrResult+=(Dimension+1);
    }  
}    

/* Extraer la diagonal principal de una matriz cuadrada */
void ExtraerDiagonal(double *A,double *DiagonalA,int Dimension)
{
    register int NdxElemento;
    register double *ptrResult;
    register double *ptrA;
    
    ptrResult=DiagonalA;
    ptrA=A;
    for(NdxElemento=0;NdxElemento<Dimension;NdxElemento++)
    {
        (*ptrResult)=(*ptrA);
        ptrA+=(Dimension+1);
        ptrResult++;
    }  
}

