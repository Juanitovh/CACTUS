
#include "f_g_compression.h"
#include "vector_funcions.h"
#include <math.h>

double f_compress_point(double * point , double *mass_gravity,  double **bounding_box)
{

    double res = 0;
    double dista = norm(point, mass_gravity,3);
    //double my_e = exp(-dista);
    double my_e = dista*dista;

    return my_e;
}

double f_compress_box(double * point , double *mass_gravity,  double **bounding_box)
{
    double res = 0;
    double aux;  
    for(int i = 0  ; i < 3 ; i++  )
    {
        aux =   point[i] + point[3] - bounding_box[1][i];
        if ( aux >0)
        {
            res += (aux * aux);
            continue;
        }


        aux = bounding_box[0][i] - point[i] - bounding_box[0][i] ;

        if ( aux >0)
            res += (aux * aux);
    }
    //double my_e = exp(-dista);

    return res;
}


double F_CostFunc_compression(vector<Strand> strands , double *mass_gravity, double **bounding_box ,double (*cost_function)(double* , double* , double **))
{
    double res = 0;
    double  * point;
    for(int i = 0 ; i < strands.size() ; i++)
    {
        Strand  &strand =strands[i];
        if (strand.n_control ==1)
            for(int  j  = 0 ; j < strand.n_control ; j++)
            {
                copy_pointer(strand.params , j*4        ,&point);
                res+= cost_function(point,  mass_gravity , bounding_box);
            }
    }
    return res;
}

///////////////

void  g_compress_point(double * point , double *g_point  ,  double *mass_gravity , double **bounding_box )
{
    double res = 0;
    double dista = norm(point, mass_gravity,3);
    //double my_e = exp(-dista);
    for(int i = 0 ; i < 3 ; i++)
        g_point[i] +=  2*dista; //* (point[i] - mass_gravity[i]);
                                //g_point[i] += - my_e* (point[i] - mass_gravity[i]);

}

void  g_compress_box(double * point , double *g_point  ,  double *mass_gravity , double **bounding_box )
{
    double res = 0;
    //double my_e = exp(-dista);
    //g_point[i] += - my_e* (point[i] - mass_gravity[i]);
    //
    double aux;  
    for(int i = 0  ; i < 3 ; i++  )
    {
        aux =   point[i]+ point[3] - bounding_box[1][i];
        if ( aux >0)
        {
            g_point[i] +=  2*aux ;
            continue;
        }

        aux = point[i] - point[3] -  bounding_box[0][i]   ;

        if ( aux < 0)
            g_point[i] +=  2*aux ;
    }
    //double my_e = exp(-dista);


}

void g_CostFunc_compression(vector<Strand> strands , double *mass_gravity , double **bounding_box, void (*cost_function)(double* , double* , double * , double **))
{
    double res = 0;
    double  * point;
    double  * g_point;
    for(int i = 0 ; i < strands.size() ; i++)
    {
        Strand  &strand =strands[i];
        if (strand.n_control ==1)
            for(int  j  = 0 ; j < strand.n_control ; j++)
            {
                copy_pointer(strand.params , j*4        ,&point);
                copy_pointer(strand.g_params , j*4        ,&g_point);
                cost_function(point, g_point,  mass_gravity , bounding_box);
            }
    }

}
