#include <cstdlib>
#include <algorithm>
#include<iostream>
#include <vector>
#include<math.h>
#include <ctime>
#include <string>
#include <fstream>
#include <stdlib.h>

#include"colors.cpp"

#include"vector_funcions.h"
#include"read_strands_txt.h"
#include "strand_class.h"
#include "Grid_space.h"
#include <omp.h>
//this is just

#include "f_g_overlap_GRID.h"
#include "f_g_2CtrlPoints.h"
#include "f_g_direction_strand.h"
#include "f_g_3CtrlPoints.h"
#include "f_g_radii.h"
#include "f_g_compression.h"
#include "optimizers.h"


#include "read_cin.h"

using namespace std;



//#include"catch_kill.cpp"
#include <unistd.h>
#include <cstdlib>
#include <signal.h>


map<string, float> scale_map;

int my_exit = 1;
void signal_callback_handler(int signum) {
  cout <<MAGENTA << "Caught KILL ... waiting 2 kill" <<RESET<< endl;
  my_exit = 0;
}



void declare_scaler(map<string, double> *scale_map)
{
  (*scale_map)["lenght"] = .02;
  (*scale_map)["lenght"] = .2;
  //(*scale_map)["lenght"] = 1;
  // (*scale_map)["lenght"] = .02;
  // (*scale_map)["lenght"] = .2;
  (*scale_map)["rad_smoth"] =.001;
  (*scale_map)["rad_fix"] =.3;
  // (*scale_map)["curve"] = .05; // original
   (*scale_map)["curve"] = .05; // original
 //(*scale_map)["curve"] = .25; // original
  // (*scale_map)["curve"] = .5; // original
  // (*scale_map)["curve"] = .05  * curve_multiplier;
  (*scale_map)["main_angle"] =  .05;
  (*scale_map)["gravity"] = 1;
  (*scale_map)["box"] = 10;

  (*scale_map)["overlap_fs"] = 1e1;
  (*scale_map)["overlap_gs"] = 1e12;


  (*scale_map)["fk_ovlp"]      = 1000000;
  (*scale_map)["fk_ovlp_last"] = 1500000;

  // cout <<endl << BOLDYELLOW<<  "Batch  overlap " <<  scale_map["fk_ovlp"] + scale_map["fk_ovlp_last"] << endl << RESET;

  (*scale_map)["fk_curve"] = 100;
  (*scale_map)["fk_lenght"] = 100;
  (*scale_map)["fk_rad"] = 100;

  //print all keys and values

}


#include "Global_optimizer.h"


// set some values:

int main(int argc, char** argv)
{

  int n_threads = 1;
  // omp_set_num_threads(n_threads);

  if (argc <2  )
  {
    cerr<<"error, parameters needed, ask for help with -h or --help"<<endl;
    return 0;
  }


  string file = argv[1] ; //  my_parameters.initial_conditions; parameters my_parameters ; 
  if (file =="-h" || file =="--help")
  {
    cout<<BOLDBLUE;
    cout<<"my_optimizer [file : txt file] [alpha: double] [#batch: int]"<<endl;

    cout<<BOLDGREEN;
    cout<<"[file]   -> strand file generated with python script"<<endl;
    cout<<"[alpha]  -> step size for gradient descensce, suggested 1"<<endl;
    cout<<"[#batch] -> size of batch to start optimizing. To optimize all strands input -1"<<endl;

    cout<<RESET<<endl;

  }

  cout <<"before global_optimizer"<<endl;
  Global_optimizer voxel_opt(argc, argv);


  cout << "---> "<< voxel_opt.alpha_multiplier << " " << voxel_opt.max_iterations << endl;
  //my_parameters.read_conf_file(config_file);



  map<string, float> scale_map;
  declare_scaler(&voxel_opt.scale_map);



  /// declare initialize box and strands
  voxel_opt.read_and_initialize_box();


int show_f_iteration = 10;


  int total_batches=  ceil ( (float) voxel_opt.n_strands  / (float) voxel_opt.batch_size );
  for (int batch_i = 0 ; batch_i < total_batches ; batch_i ++)
    {
      cout << RED<<"BATCH "<< batch_i <<" / " << total_batches << " &&  batch_size = "<< voxel_opt.batch_size << RESET <<  endl;

      // declare batch function

      voxel_opt.declare_batch_variables(voxel_opt.batch_size , batch_i);

      double last_fk =100000;
      double fk = 1e9;

      int save_iteration = 100;
      int show_f_iteration = 10;



      my_exit = 1;
      signal(SIGINT, signal_callback_handler);

      cout <<"optimizing"<<endl;



      while( voxel_opt.scale_map["fk_ovlp_last"] + voxel_opt.scale_map["fk_ovlp"] > 5000  && my_exit && voxel_opt.i_start_optim < voxel_opt.max_iterations)
 //while( voxel_opt.scale_map["fk_ovlp_last"] + voxel_opt.scale_map["fk_ovlp"] > .25*voxel_opt.batch_size  && my_exit && voxel_opt.i_start_optim < voxel_opt.max_iterations)
      // while( voxel_opt.scale_map["fk_ovlp_last"] + voxel_opt.scale_map["fk_ovlp"] > .02*voxel_opt.batch_size  && my_exit)
        {
          cout<<"("<<flush;
          voxel_opt.i_start_optim++;
          //if (i%8 ==0)
          voxel_opt.begin_environment();

          voxel_opt.step_overlap();
          voxel_opt.step_non_overlap();
          voxel_opt.saving_strands( save_iteration ,0);

          voxel_opt.print_F_and_logger( show_f_iteration );
          voxel_opt.clear_environment();
          cout<<")"<<flush;
          cout <<endl;
        }
      cout <<endl << BOLDYELLOW<<  "Batch  overlap " << voxel_opt.scale_map["fk_ovlp"] + voxel_opt.scale_map["fk_ovlp_last"] << endl << RESET;

 voxel_opt.saving_strands( save_iteration ,show_f_iteration); //

      if (!my_exit )
        break;
    }


  return 1;
}
