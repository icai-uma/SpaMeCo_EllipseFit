#include "Depurar.h"
#include "mex.h"

void ImprimirValores(double *Valores,int NumValores)
{
    int ndx;
    
    for(ndx=0;ndx<NumValores;ndx++)
        mexPrintf("%lf\n",Valores[ndx]);
    mexPrintf("\n\n");
}
  
    
void ImprimirMatriz(double *Matriz,int NumFilas,int NumCols)
{
    int NdxFila,NdxCol;
    
    for(NdxFila=0;NdxFila<NumFilas;NdxFila++)
    {
        for(NdxCol=0;NdxCol<NumCols;NdxCol++)
        {
            mexPrintf("%lf\t",Matriz[NdxCol*NumFilas+NdxFila]);
        }
        mexPrintf("\n");
    }        
    mexPrintf("\n\n");
}  
