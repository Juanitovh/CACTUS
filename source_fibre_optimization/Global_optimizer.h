
/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   strand_class.h
 * Author: jppl
 *
 * Created on February 4, 2021, 12:49 PM
 */



#ifndef GLOBAL_OPTIMIZER_H
#define GLOBAL_OPTIMIZER_H

#include <iostream>
#include <math.h>
#include<vector>
#include "strand_class.h"

#include "optimizers.h"
#include "vector_funcions.h"
#include <map>
#include <unordered_map>
#include <string>

#define pii pair<int,int>
#define v1p vector<pii>
#define v2p vector<v1p>
#define v3p vector<v2p>
#define v4p vector<v3p>

using namespace std; 

class Global_optimizer
  {
  public:
    string str_partial;
    string experiment_name ; //= "optimized_";
    string file;
    string folder;
    string out_folder;
    string profiler_ovlp;
    string profiler_non_ovlp;
    int i_start_optim;
    int max_iterations;

    map<string , double> scale_map;
    int n_total_params;
    int n_global_params ;
    int n_strands;


    int batch_size;
    int batch_id;
    float voxel_size;
    float alpha_multiplier;
    float curve_multiplier;
    double **trim_box;
    double **mass_gravity;
    double bounding_box_2[2][3];
    double bounding_box[2][3];


    double max_rad = 0;

    double *params;
    double *g_params;
    double *g_params_final;
    vector<vector < vector<double> >>  strand_list;
    vector<Strand> strands;

    vector<int> all_permuted_indexes;
    vector<int> batch_indexes;

    double step ; //= max_rad*1.5;
    int dim_grid[3];
    double offset_bb ; // 4*max_rad;

    v4p grid;
    vector<double> val_costf;


    int save_iteration;

    //		optimizer_Adadelta optimizer_ovelap(n_total_params , .02*alpha_multiplier);
    //optimizer_Adam     optimizer_ovelap(n_total_params , .02*alpha_multiplier);

    // optimizer_Adagrad  optimizer_ovelap(n_total_params , .03 *alpha_multiplier);
    // optimizer *optimizer_overlap ;// = new optimizer_Adam();
    optimizer_Adagrad *optimizer_overlap ; //(n_total_params , 1e-15*alpha_multiplier);

    //optimizer_Adagad   optimizer_lenght  (n_total_params , .05);
    optimizer_Adagrad   *optimizer_lenght ; //(n_total_params , .02);
    // optimizer *optimizer_lenght; 

    Global_optimizer();
    Global_optimizer(int argc , char **argv);
    ~Global_optimizer();
    void read_and_initialize_box();
    void declare_batch_variables(int batch_size, int batch_id);


    void delete_strand();
    void initiate(vector<vector<double>> v ,  double * global_params , double *global_g_params );
    void set_fixed_endpoints();
    void update_strands();

    void clean_gradient();

    // double f_cost_non_overlapping( );
    void grad_non_overlapping();

    double f_cost_total( );
    double f_cost_non_overlapping( );
    void grad_cost_total();
    void begin_environment();
    void clear_environment();

    void saving_strands( int save_iteration, bool final_output);
    void step_overlap();
    void step_non_overlap();

    void print_F_and_logger( int show_f_iteration);

    void write_profiler(string a);
    void check_previous_optimization();


  private:
    int get_total_parameters();
    int get_batches_parameters();
    void get_batch_indexes(int batch_index);

  };

#endif /* global optimizer */
