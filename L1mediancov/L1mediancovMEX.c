#include "mex.h"
#include "Mates.h"
#include "Depurar.h"
#include <math.h>
#include <float.h>
#include <stdlib.h>

#ifndef log2
double log2(double x)
{
	return(1.442695040888963*log(x));
}
#endif

/* Para compilar esta función, poner en el inductor de MATLAB:
>> mex L1mediancovMEX.c depurar.c nrutil.c eig.c lu.c svd.c Mates.c


function [mX,covX]=L1mediancovMEX(X,tol,maxstep);

% L1MEDIAN calculates the multivariate L1 median and an estimation of the covariance matrix
% I/O: [mX,covX]=L1median(X,tol);
%
% X is the data matrix CADA MUESTRA ES UN VECTOR COLUMNA
% tol is the convergence criterium; the iterative process stops when ||m_k - m_{k+1}|| < tol.
% maxstep is the maximum number of steps of the iterative process
%
% Ref: Hossjer and Croux (1995) "Generalizing Univariate Signed Rank Statistics for Testing
% and Estimating a Multivariate Location Parameter", Non-parametric Statistics, 4, 293-308.
% Translated from the Gauss code of Hossjer and Croux (1995) in matlab by Sabine Verboven, Antwerp University. 
% http://www.econ.kuleuven.be/public/NDBAE06/programs/pca/L1median.txt
% (accessed 29th December 2008)

% Estimacion de la matriz de covarianzas basada en la página 1666 de:
% "Multivariate Spatial U-Quantiles: a Bahadur-Kiefer Representation, a Theil-Sen Estimator for Multiple
% Regression, and a Robust Dispersion Estimator", Weihua Zhou, Robert
% Serfling: Journal of Statistical Planning and Inference 138 (2008) 1660 – 1678

*/

/* Calcular la mediana componente a componente */
void Mediana(int DimensionEspacio,int NumMuestras,double *Datos,double *Salida);

/* Calcular la mediana L1 */
void L1median(int DimensionEspacio,int NumMuestras,double Tolerancia,double MaxStep,double *Datos,double *Salida);

/* Calcular la estimación de la matriz de covarianzas basada en la mediana L1 */
void L1cov(int DimensionEspacio,int NumMuestras,double Tolerancia,double MaxStep,double *Datos,double *Salida);


void mexFunction(int nlhs, mxArray* plhs[],
                 int nrhs, const mxArray* prhs[])
{   
    double *ptrMuestras,*ptrL1median,*ptrL1cov;
    double Tolerancia;
    int DimensionEspacio,NumMuestras,MaxStep;
    const int *DimMuestras;            

    /* Obtener variables de trabajo */
    ptrMuestras=mxGetPr(prhs[0]);
    if (nrhs>1)
    {
         Tolerancia=*mxGetPr(prhs[1]);
         if (nrhs>2)
         {
              MaxStep=*mxGetPr(prhs[2]);
         }
         else
         {
              MaxStep=200;
         }
    }
    else
    {
         Tolerancia=1.0E-8;
         MaxStep=200;
    }
    DimMuestras=mxGetDimensions(prhs[0]);
    DimensionEspacio=DimMuestras[0];
    NumMuestras=DimMuestras[1];
    
    /* Crear matrices de salida */
    plhs[0]=mxCreateDoubleMatrix(DimensionEspacio,1,mxREAL);
    ptrL1median=mxGetPr(plhs[0]);
    if (nlhs>1)
    {
         plhs[1]=mxCreateDoubleMatrix(DimensionEspacio,DimensionEspacio,mxREAL);
         ptrL1cov=mxGetPr(plhs[1]);
    }
    
    /* Realizar cálculos */
    L1median(DimensionEspacio,NumMuestras,Tolerancia,MaxStep,ptrMuestras,ptrL1median);
    if (nlhs>1)
    {
         L1cov(DimensionEspacio,NumMuestras,Tolerancia,MaxStep,ptrMuestras,ptrL1cov);
    }
}    

/* Comparar dos valores double para qsort() */
int CompararDoubles(const void * a, const void * b)
{
    double MiA,MiB;
    
    MiA=*(double*)a;
    MiB=*(double*)b;
    if (MiA>MiB)
    {
         return 1;
    }
    else
    {
        if (MiA<MiB)
        {
             return -1;
        }
        else
        {
            return 0;
        }
    }
}


/* Calcular la mediana componente a componente */
void Mediana(int DimensionEspacio,int NumMuestras,double *Datos,double *Salida)
{
     int NdxDim,NdxMuestra;
     double *ptrLista;
     
     ptrLista=mxMalloc(NumMuestras*sizeof(double));
     
     /* Procesar cada componente por separado */
     for(NdxDim=0;NdxDim<DimensionEspacio;NdxDim++)
     {
          /* Colocar los datos en una lista y ordenarlos */
          for(NdxMuestra=0;NdxMuestra<NumMuestras;NdxMuestra++)
          {
               ptrLista[NdxMuestra]=Datos[NdxDim+NdxMuestra*DimensionEspacio];
          }
          qsort(ptrLista,NumMuestras,sizeof(double),CompararDoubles);
          
          /* Averiguar cuál es la mediana */
          if ((NumMuestras%2)==0)
          {
               Salida[NdxDim]=0.5*(ptrLista[NumMuestras/2]+ptrLista[(NumMuestras/2)-1]);
          }
          else
          {
               Salida[NdxDim]=ptrLista[(NumMuestras-1)/2];
          }
     }
     
     mxFree(ptrLista);
     
}

/* Calcular la función objetivo mrobj */
double MRObj(int DimensionEspacio,int NumMuestras,double *Datos,double *MiMed)
{
     int NdxDim,NdxMuestra;
     double Salida,MiNorma;
     double *ptrDif;
     
     ptrDif=mxMalloc(DimensionEspacio*1*sizeof(double));
     
     Salida=0.0;
     for(NdxMuestra=0;NdxMuestra<NumMuestras;NdxMuestra++)
     {
         DiferenciaMatricial(Datos+NdxMuestra*DimensionEspacio,
              MiMed,ptrDif,DimensionEspacio,1);
         NormaCuadrado(ptrDif,&MiNorma,DimensionEspacio);
         Salida+=sqrt(MiNorma);
     }

     mxFree(ptrDif);                                                               
     return Salida;
}

/* Calcular la mediana L1 */
void L1median(int DimensionEspacio,int NumMuestras,double Tolerancia,double MaxStep,double *Datos,double *Salida)
{
     int NdxStep,NdxDim,NdxMuestra,MyStep,MaxHalf;
     double *ptrMold,*ptrXext,*ptrDiferencia,*ptrW,*ptrDelta;
     double SumaW,MiNorma;
     
     
     /* Reservar memoria auxiliar */
     ptrMold=mxMalloc(DimensionEspacio*1*sizeof(double));
     ptrDelta=mxMalloc(DimensionEspacio*1*sizeof(double));
     ptrDiferencia=mxMalloc(DimensionEspacio*NumMuestras*sizeof(double));
     ptrXext=mxMalloc((DimensionEspacio+1)*NumMuestras*sizeof(double));
     ptrW=mxMalloc(NumMuestras*1*sizeof(double));
     
     /* Inicializar Salida==m con la mediana hallada componente a componente */
     Mediana(DimensionEspacio,NumMuestras,Datos,Salida);
     
     /* Algoritmo iterativo */
     NdxStep=1;
     while(NdxStep<=MaxStep)
     {
          /* mold=m; */
          memcpy(ptrMold,Salida,sizeof(double)*DimensionEspacio);
          
          /* Ordenar las muestras según la norma del error cuadrático respecto a la mediana */
          for(NdxMuestra=0;NdxMuestra<NumMuestras;NdxMuestra++)
          {
               /* Colocar el dato (vector columna) */
               memcpy(ptrXext+NdxMuestra*(DimensionEspacio+1)+1,
                  Datos+NdxMuestra*DimensionEspacio,
                  sizeof(double)*DimensionEspacio);
               /* Colocar la norma del error como primera componente del vector columna */
               DiferenciaMatricial(Datos+NdxMuestra*DimensionEspacio,Salida,
                    ptrDiferencia,DimensionEspacio,1);
               NormaCuadrado(ptrDiferencia,
                  ptrXext+NdxMuestra*(DimensionEspacio+1),DimensionEspacio);
               ptrXext[NdxMuestra*(DimensionEspacio+1)]=sqrt(ptrXext[NdxMuestra*(DimensionEspacio+1)]);
          }
          qsort(ptrXext,NumMuestras,(DimensionEspacio+1)*sizeof(double),CompararDoubles);
          
          /* Ahora ptrXext+NdxMuestra*(DimensionEspacio+1) es la norma del error de la muestra NdxMuestra,
          y la muestra propiamente dicha está en ptrXext+NdxMuestra*(DimensionEspacio+1)+1.*/
          
          /* Ver cuáles de esas normas son no nulas, y hallamos sus inversas. Para las
          que sean nulas, ponemos cero */
          SumaW=0.0;
          for(NdxMuestra=0;NdxMuestra<NumMuestras;NdxMuestra++)
          {
               if(ptrXext[NdxMuestra*(DimensionEspacio+1)]>0.0)
               {
                    ptrW[NdxMuestra]=1.0/ptrXext[NdxMuestra*(DimensionEspacio+1)];
                    SumaW+=ptrW[NdxMuestra];
               }
               else
               {
                   ptrW[NdxMuestra]=0.0;
               }
          }

          if (SumaW<=0.0)
          {
               mexWarnMsgTxt("L1mediancovMEX: Datos identicos?");
               break;
          }
          
          /* Hallar (X-repmat(m,n,1)).*repmat(w,1,p) sobre ptrDiferencia */
          for(NdxMuestra=0;NdxMuestra<NumMuestras;NdxMuestra++)
          {
               DiferenciaMatricial(ptrXext+NdxMuestra*(DimensionEspacio+1)+1,Salida,
                    ptrDiferencia+NdxMuestra*DimensionEspacio,DimensionEspacio,1);
               ProductoEscalarMatriz(ptrW[NdxMuestra],ptrDiferencia+NdxMuestra*DimensionEspacio,
                   ptrDiferencia+NdxMuestra*DimensionEspacio,DimensionEspacio,1);
          }
          
          /* Hallar nd==MiNorma==norme(delta), donde delta es de dimension 1xDimensionEspacio */
          MiNorma=0.0;
          for(NdxDim=0;NdxDim<DimensionEspacio;NdxDim++)
          {
               ptrDelta[NdxDim]=0.0;
               for(NdxMuestra=0;NdxMuestra<NumMuestras;NdxMuestra++)
               {
                    ptrDelta[NdxDim]+=ptrDiferencia[NdxDim+NdxMuestra*DimensionEspacio];
               }
               ptrDelta[NdxDim]/=SumaW;
               MiNorma+=ptrDelta[NdxDim]*ptrDelta[NdxDim];
          }
          MiNorma=sqrt(MiNorma);
          
          /* Ver si nd<tol */
          if (MiNorma<Tolerancia)
          {
               MaxHalf=0.0;
          }
          else
          {
              MaxHalf=log2(MiNorma/Tolerancia);
          }
          
          /* Nueva estimación: m=mold+delta */
          SumaMatricial(ptrMold,ptrDelta,Salida,DimensionEspacio,1);
          
          /* Iteracion interior */
          MyStep=0;
          while ( (MRObj(DimensionEspacio,NumMuestras,Datos,Salida)>=
                      MRObj(DimensionEspacio,NumMuestras,Datos,ptrMold)) && (MyStep<=MaxHalf) )
          {
                MyStep=MyStep+1;
                /* m=mold+delta./(2^nstep); */
                ProductoEscalarMatriz(1.0/pow(2.0,(double)MyStep),ptrDelta,ptrDiferencia,
                    DimensionEspacio,1);
                SumaMatricial(ptrMold,ptrDiferencia,Salida,DimensionEspacio,1);
          }
          
          if (MyStep>MaxHalf)
          {
              memcpy(Salida,ptrMold,sizeof(double)*DimensionEspacio);
              break;
          }
          
          NdxStep++;
     }
     
     /* Ver si no se ha conseguido la convergencia */
     if (NdxStep>MaxStep)
     {          
          mexWarnMsgTxt("L1mediancovMEX: Iteracion fallida");
     }
     
     /* Liberar memoria auxiliar */
     mxFree(ptrMold);
     mxFree(ptrXext);
     mxFree(ptrDiferencia);
     mxFree(ptrW);
     mxFree(ptrDelta);
}

/* Calcular la estimación de la matriz de covarianzas basada en la mediana L1 */
void L1cov(int DimensionEspacio,int NumMuestras,double Tolerancia,double MaxStep,double *Datos,double *Salida)
{
     int NdxMuestra1,NdxMuestra2,NumProductos,NdxProducto;
     double *ptrDiferencia,*ptrProductos,*ptrVecino1;
     
     NumProductos=(NumMuestras*(NumMuestras-1))/2;
     ptrDiferencia=mxMalloc(DimensionEspacio*sizeof(double));
     ptrProductos=mxMalloc(DimensionEspacio*DimensionEspacio*NumProductos*sizeof(double));
     
     NdxProducto=0;
     for(NdxMuestra1=0;NdxMuestra1<(NumMuestras-1);NdxMuestra1++)
     {
          /* MiVecino1=MuestrasVecinas(:,NdxVecino1); */
          ptrVecino1=Datos+NdxMuestra1*DimensionEspacio;
          for(NdxMuestra2=NdxMuestra1+1;NdxMuestra2<NumMuestras;NdxMuestra2++)
          {
               /* MiVecino2=MuestrasVecinas(:,NdxVecino2);
                    Diferencia=MiVecino1-MiVecino2; */ 
               DiferenciaMatricial(ptrVecino1,
                    Datos+NdxMuestra2*DimensionEspacio,ptrDiferencia,DimensionEspacio,1);
               ProductoMatricial(ptrDiferencia,ptrDiferencia,
                   ptrProductos+NdxProducto*DimensionEspacio*DimensionEspacio,
                   DimensionEspacio,1,DimensionEspacio);
               NdxProducto++;
          }
     }
     ProductoEscalarMatriz(0.5,ptrProductos,ptrProductos,DimensionEspacio*DimensionEspacio,NumProductos);
     
     /* Hallar la mediana L1 de los productos */
     L1median(DimensionEspacio*DimensionEspacio,NumProductos,Tolerancia,MaxStep,
          ptrProductos,Salida);
     
     mxFree(ptrDiferencia);
     mxFree(ptrProductos);
}


