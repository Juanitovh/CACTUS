/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   optimizers.h
 * Author: jppl
 *
 * Created on April 9, 2021, 11:14 AM
 */

#include"vector_funcions.h"
#ifndef OPTIMIZERS_H
#define OPTIMIZERS_H

class optimizer{
public:

  double alfa;
  vf gt ;
  int n;
  void set_gradient(double * gradient  );
  optimizer();
  virtual ~optimizer(){};

  optimizer(int n_total , double alfa);

  // virtual void initialize( int n_total_params , double alfa1);
  void copy_to_pointer_scaled(double * gradient , double alpha);

  virtual void compute_step(double * gradient );
  virtual void reset_optimizer();
};


class optimizer_Adam : public optimizer
  {
  public:
    vf pt;
    vf qt ;
    vf p_hat ;
    vf q_hat;
    vf wt;


    double m1 = 0.9;
    double m2 = 0.999;

    double m1t = 1;
    double m2t = 1;

    double eps = 1e-8;
    vf eps_v ;

    optimizer_Adam( int n_total_params , double alfa1);
    optimizer_Adam();
    ~optimizer_Adam(){};

    // void initialize( int n_total_params , double alfa1);
    void compute_step(double * gradient );
    void reset_optimizer();

  };



class optimizer_Adadelta : public optimizer
  {
  public:

    vf GT;
    vf eps_v ;
    double alfa_adadelta = 0.1;
    double eta = 0.9;
    double eps = 1e-8;

    optimizer_Adadelta( int n_total_params  , double alfa1 );
    optimizer_Adadelta();
    ~optimizer_Adadelta(){};

    // void initialize( int n_total_params  , double alfa1 );
    void compute_step(double * gradient );

    void reset_optimizer();


  };



class optimizer_Adagrad : public optimizer
  {
  public:

    vf GT;
    vf eps_v ;
    double alfa_adadelta = 0.1;

    double eps = 1e-8;

    optimizer_Adagrad( int n_total_params , double alfa1);
    optimizer_Adagrad();
    ~optimizer_Adagrad(){};

    // void initialize( int n_total_params , double alfa1);
    void compute_step(double * gradient );
    void reset_optimizer();

  };




class optimizer_fixed_alpha :public optimizer
  {
  public:

    double alfa;

    double eps = 1e-8;

    optimizer_fixed_alpha();
    ~optimizer_fixed_alpha(){};
    optimizer_fixed_alpha( int n_total_params , double alfa1);
    void initialize( int n_total_params , double alfa1);

    void compute_step(double * gradient );
    void reset_optimizer();

  };



#endif /* OPTIMIZERS_H */
