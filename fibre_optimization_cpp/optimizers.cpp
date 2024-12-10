/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   optimizers.cpp
 * Author: jppl
 * 
 * Created on April 9, 2021, 11:14 AM
 */

#include <omp.h>
#include "omp_config.h"
#include "optimizers.h"



optimizer::optimizer()
{

}


optimizer::optimizer(int n , double alfa)
{

}

void optimizer::compute_step(double * gradient )
{

}
void optimizer::reset_optimizer()
{
  
}

void optimizer::copy_to_pointer_scaled(double * v1 , double alpha  )
{
#pragma omp parallel for num_threads(GLOBAL_NUM_THREADS)
	for(int i = 0 ; i < gt.size() ; i++)
	{
		//cout<<v1[i]<<" -*- "<<alfa*v2[i] <<" , ";
		v1[i] = alpha*gt[i];
	}
	//cout<<endl<<endl;
}

void optimizer::set_gradient(double * gradient )
{
#pragma omp parallel for num_threads(GLOBAL_NUM_THREADS)
    for(int i = 0 ; i < gt.size() ; i++)
        gt[i] = gradient[i];
}

optimizer_Adam::optimizer_Adam( int n_total_params , double alfa1)
{
    gt.assign (n_total_params, 0);
    pt = 0*gt;
    qt = 0*gt;
    eps_v.assign(n_total_params , eps);
    alfa = alfa1;
}

optimizer_Adam::optimizer_Adam()
{

}




void optimizer_Adam::compute_step(double * gradient)
{
    set_gradient(gradient );

    m1t *= m1;
    m2t *= m2;

    pt = m1*pt + (1-m1)*gt;
    qt = m2*qt + (1-m2)*(gt*gt);

    p_hat = (1/(1-m1t))*pt;
    q_hat = (1/(1-m2t))*qt;

    gt =  alfa*p_hat/(sqrt_v(q_hat) + eps_v);
}

void optimizer_Adam::reset_optimizer()
{
    m1 = 0.9;
    m2 = 0.999;

    m1t = 1;
    m2t = 1;
    pt = 0*gt;
    qt = 0*gt;

}





optimizer_Adadelta::optimizer_Adadelta( int n_total_params  , double alfa1 )
{
    gt.assign (n_total_params, 0);
    eps_v.assign(n_total_params , eps);
    GT = gt-gt;
    alfa = alfa1;
}





void optimizer_Adadelta::compute_step(double * gradient )
{

    set_gradient(gradient );

    GT = eta*gt*gt + (1-eta)*GT;

    gt = alfa*gt/(sqrt_v(GT) + eps_v);

}



void optimizer_Adadelta::reset_optimizer()
{
    GT = GT -GT;

}





optimizer_Adagrad::optimizer_Adagrad( int n_total_params , double alfa1)
{
    gt.assign (n_total_params, 0);
    eps_v.assign(n_total_params , eps);
    GT = gt-gt;
    alfa = alfa1;
}





void optimizer_Adagrad::compute_step(double * gradient )
{

    set_gradient(gradient );

    GT =GT + gt*gt ;

    gt = alfa* gt/(sqrt_v(GT) + eps_v);

}



void optimizer_Adagrad::reset_optimizer()
{
    GT = GT -GT;

}



optimizer_fixed_alpha::optimizer_fixed_alpha( int n_total_params , double alfa1)
{
    gt.assign (n_total_params, 0);
    alfa = alfa1;
}





void optimizer_fixed_alpha::compute_step(double * gradient )
{
    set_gradient(gradient );
    gt = alfa* gt;
}



void optimizer_fixed_alpha::reset_optimizer()
{


}



