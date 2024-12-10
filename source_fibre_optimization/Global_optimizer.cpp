
#include <cstdlib>

#include <algorithm>
#include <ctime>
#include <fstream>
#include <iostream>
#include <math.h>
#include <stdlib.h>
#include <string>
#include <vector>

#include "colors.cpp"

#include "Grid_space.h"
#include "read_strands_txt.h"
#include "strand_class.h"
#include "vector_funcions.h"
#include <omp.h>
#include "omp_config.h"

// this is just
#include "Global_optimizer.h"
#include "Grid_space.h"
#include "f_g_2CtrlPoints.h"
#include "f_g_3CtrlPoints.h"
#include "f_g_compression.h"
#include "f_g_direction_strand.h"
#include "f_g_overlap_GRID.h"
#include "f_g_radii.h"
#include "optimizers.h"

#include "read_cin.h"

#include <ctime>

#include <algorithm>
using namespace std;

//#include"catch_kill.cpp"create_matrix
#include <cstdlib>
#include <iostream>
#include <signal.h>
#include <unistd.h>


#include <iostream>
#include <chrono>
#include <ctime>
#include <iomanip>


int Global_optimizer::get_total_parameters() {
  int total = 0;
  for (int i = 0; i < strand_list.size(); i++)
    total += strand_list[i].size() * 4;
  return total;
}

Global_optimizer::Global_optimizer(int argc, char **argv) {
  file = argv[1];
  cout << "*** ---- >> current file " << file << endl;
  // split file by character .
  vector<string> splits = split_string(file, ".");
  folder = splits[0];

  // make folder if it doesn't exist
  create_folder(folder);

  alpha_multiplier = atof(argv[2]);
  batch_size = atoi(argv[3]);
  max_iterations = atof(argv[4]);

  str_partial = ".partial";
  experiment_name = folder + "/optimized_";
  val_costf.assign(5, -1);
  profiler_ovlp = folder + "/profiler_ovlp.log";
  profiler_non_ovlp = folder + "/profiler_non_ovlp.log";

  check_previous_optimization();
}

void Global_optimizer::check_previous_optimization() {
  vector<string> files_partial = list_files(folder, str_partial);

  if (files_partial.size() == 0) {
    cout << "No previous optimization found, starting at 0" << endl;
    i_start_optim = 0;
  } else {
    int last = files_partial.size() - 1;

    file = folder + "/" + files_partial[last];
    cout << "Previous optimization found, starting at " << file << endl;

    string last_file = files_partial[last];
    vector<string> splits = split_string(last_file, ".");
    last_file = splits[0];
    splits = split_string(last_file, "_");
    last_file = splits[1];

    i_start_optim = atoi(last_file.c_str());
    cout << "Found previous optimization at " << i_start_optim << endl;
  }
}

Global_optimizer::~Global_optimizer() {

  cout << "cleaning" << endl;
  for (int i = 0; i < strands.size(); i++)
    strands[i].delete_strand();

  delete[] params;
  // delete [] params2;
  delete[] g_params_final;
  delete[] g_params;
  delete[] trim_box[0];
  delete[] trim_box[1];
  delete[] trim_box;
}

void Global_optimizer::read_and_initialize_box() {
  // scale_map["curve"] *= curve_multiplier;
  scale_map["overlap_gs"] *= alpha_multiplier;

  cout << "Reading strands" << endl;
  strand_list = read_generic_cylinder_list(file, &voxel_size);
  cout << " Readed " << strand_list.size() << " -strands  " << n_global_params
       << " - parameters " << n_global_params / (float)(strand_list.size() * 4)
       << "  -- average CP in strands" << endl;
  batch_size = strand_list.size();

  // optimizer_Adagrad  optimizer_ovelap(n_total_params , .03
  // *alpha_multiplier);

  trim_box = new double *[2];

  trim_box[0] = new double[3];
  trim_box[1] = new double[3];

  for (int i = 0; i < 3; i++) {
    trim_box[0][i] = -voxel_size / 2;
    trim_box[1][i] = voxel_size / 2;
  }

  n_global_params = get_total_parameters();
  n_total_params = n_global_params;
  cout << " ### " << n_global_params << " ### " << endl;

  optimizer_overlap =
      new optimizer_Adagrad(n_total_params, .03 * alpha_multiplier);
  // optimizer_lenght = new optimizer_Adagrad  (n_total_params , .02);
  optimizer_lenght = new optimizer_Adagrad(n_total_params, .005);

  n_strands = strand_list.size();

  all_permuted_indexes.clear();

  for (int i = 0; i < n_strands; i++)
    all_permuted_indexes.push_back(i);

  random_shuffle(all_permuted_indexes.begin(), all_permuted_indexes.end());
}

void Global_optimizer::get_batch_indexes(int batch_index) {
  // batch_inde
  this->batch_id = batch_index;

  int start_batch = n_strands / batch_size;
  start_batch = batch_size * batch_id;

  int end_batch = min(batch_size, n_strands - start_batch);
  batch_indexes.clear();
  for (int i = 0; i < end_batch; i++) {
    // indexes.push_back(start_batch + i);
    batch_indexes.push_back(start_batch + i);
  }
}

void Global_optimizer::declare_batch_variables(int batch_size, int batch_id) {

  cout << RED << "batch size " << batch_size << RESET << endl;
  int total_batches = ceil((float)n_strands / (float)batch_size);

  get_batch_indexes(batch_id);

  n_total_params = get_batches_parameters();

  params = new double[n_total_params];
  g_params = new double[n_total_params];
  g_params_final = new double[n_total_params];

  int current_index = 0;

  for (int i = 0; i < batch_indexes.size(); i++) {
    Strand auxliar_strand;
    auxliar_strand.initiate(strand_list[batch_indexes[i]],
                            params + current_index, g_params + current_index);
    auxliar_strand.set_fixed_endpoints();
    strands.push_back(auxliar_strand);
    current_index += auxliar_strand.n_control * 4;
  }

  for (int i = 0; i < n_total_params / 4; i++) {

    max_rad = max(max_rad, params[4 * i + 3]);
  }

  // cout << "max_rad " << max_rad <<endl;

  // step = max_rad * 1.5;
  // step = max_rad * 2;
  // step = max_rad * 5;
  step = max_rad * 3;
  step = 4;
  step = 2;
  // step = max_rad * 4;
  // step = max_rad * 10;
  // step = max_rad * 5;
  // step = 2; // max_rad*1.5;:
  cout << "step " << step << endl;
  offset_bb = 4 * max_rad;

  get_bounding_box(params, n_total_params, bounding_box, offset_bb);

  // print bounding box

  for (int i = 0; i < 3; i++) {
    bounding_box_2[0][i] = bounding_box[0][i] - step / 2.0;
    bounding_box_2[1][i] = bounding_box[1][i] - step / 2.0;

    cout << "bounding box " << i << " " << bounding_box[0][i] << " "
         << bounding_box[1][i] << endl;
  }

  double mass_gravity[3];
  for (int i = 0; i < 3; i++) {
    dim_grid[i] = ceil((bounding_box[1][i] - bounding_box[0][i]) / step);
    mass_gravity[i] = (bounding_box[1][i] - bounding_box[0][i]) / 2.0;
    //								cout
    //<<bounding_box[0][i]
    // cout << bounding_box[1][i] << dim_grid[i] << " parts" << endl;
    cout << dim_grid[i] << " parts" << endl;
  }

  grid = create_grid(dim_grid);
  // cout << "grid created" <<endl;
  // cout <<grid.size()<<endl;
  // cout <<grid[0].size()<<endl;
  // cout <<grid[0][0].size()<<endl;
}
void Global_optimizer::begin_environment() {

  cout << "[ " << flush;
    auto wall_start = std::chrono::high_resolution_clock::now();
    std::clock_t cpu_start = std::clock();   

  if (i_start_optim % 2 == 1)
    fill_grid(grid, strands, bounding_box, dim_grid, step);
  else
    fill_grid(grid, strands, bounding_box_2, dim_grid, step);



     std::clock_t cpu_end = std::clock();
    auto wall_end = std::chrono::high_resolution_clock::now();
   double cpu_elapsed_secs = double(cpu_end - cpu_start) / CLOCKS_PER_SEC;
    std::chrono::duration<double> wall_elapsed_secs = wall_end - wall_start;
    cout << std::fixed << std::setprecision(1);


  cout << cpu_elapsed_secs << " , " << wall_elapsed_secs.count()  << " ]" << flush;
}

void Global_optimizer::clear_environment() {

  //cout << "&" << flush;
  clear_grid(grid, dim_grid);

}

void Global_optimizer::step_overlap() {

  if (i_start_optim % 2 == 0) {
    optimizer_overlap->reset_optimizer();
    //cout << BLUE << " reset opt " << RESET;
  }

  // cout<<"%% filled"<<endl;

  cout << "<" << flush;

    auto wall_start = std::chrono::high_resolution_clock::now();
    std::clock_t cpu_start = std::clock();   



  grad_cost_total();

     std::clock_t cpu_end = std::clock();
    auto wall_end = std::chrono::high_resolution_clock::now();
   double cpu_elapsed_secs = double(cpu_end - cpu_start) / CLOCKS_PER_SEC;
    std::chrono::duration<double> wall_elapsed_secs = wall_end - wall_start;
    cout << std::fixed << std::setprecision(1);

  cout << cpu_elapsed_secs << " , " << wall_elapsed_secs.count()  << " >" << flush;

  optimizer_overlap->compute_step(g_params_final);

  // for(int kk = 0 ; kk < n_total_params  ; kk++)
  //		g_params_final[kk] = optimizer_ovelap.gt[kk];

  resta_alfav(params, optimizer_overlap->gt, 1, n_total_params);
}

void Global_optimizer::step_non_overlap() {

 //if (i_start_optim % 4 == 1)
  if (i_start_optim % 3 == 1)
  // if (1)
  {
    
    cout << "{" << flush;

    auto wall_start = std::chrono::high_resolution_clock::now();
    std::clock_t cpu_start = std::clock();   

    optimizer_lenght->reset_optimizer();

    int aux = 0;

    double last_fk_local = 11;
    double relative_diff_local = 10;
    // double fk_local = f_cost_non_overlapping() ;

    f_cost_non_overlapping();
    while (aux < 4) {
      if (aux % 2 == 1)
        optimizer_lenght->reset_optimizer();

      cout << "-" << flush;
      // relative_diff_local = abs(fk_local - last_fk_local)/fk_local;
      aux++;
      grad_non_overlapping();

      optimizer_lenght->compute_step(g_params_final);
      resta_alfav(params, optimizer_lenght->gt, 1, n_total_params);
      // optimizer_lenght->copy_to_pointer_scaled(g_params_final , 1);
      // resta_alfav(params, optimizer_lenght.gt, 1, n_total_params);

      // last_fk_local  = fk_local;
    }
    f_cost_non_overlapping();

     std::clock_t cpu_end = std::clock();
    auto wall_end = std::chrono::high_resolution_clock::now();
   double cpu_elapsed_secs = double(cpu_end - cpu_start) / CLOCKS_PER_SEC;
    std::chrono::duration<double> wall_elapsed_secs = wall_end - wall_start;
    cout << std::fixed << std::setprecision(1);

  cout << cpu_elapsed_secs << " , " << wall_elapsed_secs.count()  << " }" << flush;
    //cout << "}" << flush;
  }
}

void Global_optimizer::print_F_and_logger(int show_f_iteration) {

  if (i_start_optim % show_f_iteration == 0) {
    // double fk_local = f_cost_non_overlapping() ;

    // last_fk = fk;
    cout << BOLDMAGENTA<<" "  << i_start_optim << "th-step " << RESET;
    cout << "¿" << flush;
    //start counting time

    auto wall_start = std::chrono::high_resolution_clock::now();
    std::clock_t cpu_start = std::clock();   


    double fk = f_cost_total();

     std::clock_t cpu_end = std::clock();
    auto wall_end = std::chrono::high_resolution_clock::now();
   double cpu_elapsed_secs = double(cpu_end - cpu_start) / CLOCKS_PER_SEC;
    std::chrono::duration<double> wall_elapsed_secs = wall_end - wall_start;
    cout << std::fixed << std::setprecision(1);

  //  cout << "?" << flush;

  cout << cpu_elapsed_secs << " , " << wall_elapsed_secs.count()  << " ?" << flush;
    scale_map["fk_ovlp_last"] = (double)scale_map["fk_ovlp"];
    scale_map["fk_ovlp"] = fk;

    // write into the loger, or something
    //  write_profiler(profiler_ovlp);
  }
}

void Global_optimizer::saving_strands(int save_iteration, bool final_output) {

  if (i_start_optim % save_iteration == 0 && final_output == false) {
    string old_string = to_string(i_start_optim);
    string new_name =
        string(5 - old_string.length(), '0') + old_string + str_partial;

    // update_strands();
    write_generic_cylinder_list(strands, experiment_name + new_name,
                                voxel_size);
  }
  if (final_output) {

    // update_strands();
    write_generic_cylinder_list(strands, experiment_name + "final.txt",
                                voxel_size);
  }
}

// void Global_optimizer::update_strands() {
//
//   cout << "updating strands!" << endl;
//   for (int index_i = 0; index_i < batch_indexes.size(); index_i++) {
//     for (int i = 0; i < strands[index_i].n_control; i++) {
//       for (int j = 0; j < 4; j++) {
//         strand_list[batch_indexes[index_i]][i][j] =
//             strands[index_i].params[i * 4 + j];
//       }
//     }
//   }
// }

int Global_optimizer::get_batches_parameters() {
  int total = 0;
  for (int i = 0; i < batch_indexes.size(); i++)
    total += strand_list[batch_indexes[i]].size() * 4;
  return total;
}

void Global_optimizer::write_profiler(string file) {

  string filename(file);
  ofstream file_out;

  file_out.open(filename, ios_base::app);
  if (val_costf[0] == -1) {
    cout << "############################################## emergency?!"
         << endl;
    return;
  }

  int i = 0;
  while (val_costf[i] != -1) {
    file_out << val_costf[i] << " ";
    val_costf[i] = -1;
    i++;
  }
  file_out << endl;
}

void Global_optimizer::grad_cost_total() {
  // fill_n(g_params, n_total_params, 0);
  // fill_n(g_params_final, n_total_params, 0);

 parallel_fill_n(g_params, n_total_params, 0);
 parallel_fill_n(g_params_final, n_total_params, 0);

  double flag_g;
  bool compress = 0;

  g_overlap_all_strands_grid(grid, strands, g_params, n_total_params, compress);

  // flag_g = scale_map["overlap_gs"];
  resta_alfav(g_params_final, g_params, -1 * scale_map["overlap_gs"],
              n_total_params);
}

double Global_optimizer::f_cost_total() {
  double res2 =
      F_AllStrands_2CtrlPoints(strands, F_CostFunc_Lenght_2CntrlPoints_hooke) *
      scale_map["lenght"];
  double res3 = F_AllStrands_2CtrlPoints(strands, F_cost_angle_strand) *
                scale_map["main_angle"];
  // double res3=0 ;// F_AllStrands_2CtrlPoints(strands ,
  // F_CostFunc_Rads_2CntrlPoints) * scale_rad_smoth;
  double res4 = f_cost_rads(strands) * scale_map["rad_fix"];
  double res5 =
      F_AllStrands_3CtrlPoints(strands, F_CostFunc_Curve_3CntrlPoints) *
      scale_map["curve"];

  double res6 = 0; // F_CostFunc_compression(strands ,  mass_gravity , trim_box
                   // ,   f_compress_point scale_map["gravity"];
  double res7 = 0; // F_CostFunc_compression(strands ,  mass_gravity, trim_box
                   // ,f_compress_box ) * scale_map["box"];

  double res1 =
      f_overlap_all_strands_grid(grid, strands) * scale_map["overlap_fs"];
  val_costf[0] = res1;

  cout << endl
       << BOLDGREEN << " Ovlp: " << res1 << "  Len: " << res2
       << "  Angº: " << res3 << "  FxR: " << res4 << "  Curve: " << res5
       << "  Comp: " << res6 << " Box: " << res7
       << "  ===  " << res1 + res2 + res3 + res4 + res5 + res6 + res7 << " "
       << RESET << endl;

  write_profiler(profiler_ovlp);

  val_costf[0] = res2;
  val_costf[1] = res5;
  write_profiler(profiler_non_ovlp);
  return res1; //+ res2+res3 +res4+ res5;
}

double Global_optimizer::f_cost_non_overlapping() {
  double res2 =
      F_AllStrands_2CtrlPoints(strands, F_CostFunc_Lenght_2CntrlPoints_hooke) *
      scale_map["lenght"];
  double res3 = 0; //  F_AllStrands_2CtrlPoints(strands , F_cost_angle_strand) *
                   //  scale_map["main_angle"];
  // double res3=0 ;// F_AllStrands_2CtrlPoints(strands ,
  // F_CostFunc_Rads_2CntrlPoints) * scale_rad_smoth;
  double res4 = 0; // f_cost_rads(strands) * scale_map["rad_fix"];
  double res5 =
      F_AllStrands_3CtrlPoints(strands, F_CostFunc_Curve_3CntrlPoints) *
      scale_map["curve"];

  double res6 = 0; // F_CostFunc_compression(strands ,  mass_gravity , trim_box
                   // ,   f_compress_point scale_map["gravity"];
  double res7 = 0; // F_CostFunc_compression(strands ,  mass_gravity, trim_box
                   // ,f_compress_box ) * scale_map["box"];

  double res1 = 0; // f_overlap_all_strands_grid( grid ,  strands)
                   // *scale_map["overlap_fs"];
  // cout <<endl << BOLDGREEN   << " Ovlp: "  << res1 << "  Len: " << res2   <<
  // "  Angº: " << res3  << "  FxR: "<< res4 <<"  Curve: "<< res5  << "  Comp: "
  // << res6 <<" Box: " << res7   << "  ===  "<< res1 + res2  + res3  + res4  +
  // res5 + res6+ res7<< " "<< RESET<<endl;
  val_costf[0] = res2;
  val_costf[1] = res5;
  write_profiler(profiler_non_ovlp);

  return res1; //+ res2+res3 +res4+ res5;
}

void Global_optimizer::grad_non_overlapping() {
  // fill_n(g_params, n_total_params, 0);
  // fill_n(g_params_final, n_total_params, 0);

  parallel_fill_n(g_params, n_total_params, 0);
  parallel_fill_n(g_params_final, n_total_params, 0);

  G_AllStrands_2CtrlPoints(strands, G_CostFunc_Lenght_2CntrlPoints_hooke);
  resta_alfav(g_params_final, g_params, -1 * scale_map["lenght"],
              n_total_params);
  // resta_alfav(g_params_final, g_params, 0, n_total_params);

  // fill_n(g_params, n_total_params, 0);

  parallel_fill_n(g_params, n_total_params, 0);
  
  G_AllStrands_3CtrlPoints(strands, G_CostFunc_Curve_3CntrlPoints);
  resta_alfav(g_params_final, g_params, -1 * scale_map["curve"],
              n_total_params);

  // fill_n(g_params , n_total_params , 0);
  // G_AllStrands_3CtrlPoints( strands , G_CostFunc_angle_strand );
  // resta_alfav(g_params_final, g_params, -scale_map["1main_angle"],
  // n_total_params);

  /*
       fill_n(g_params , n_total_params , 0);
       G_AllStrands_2CtrlPoints( strands , G_CostFunc_Rads_2CntrlPoints );
       resta_alfav(g_params_final, g_params, -1*scale_map["rad_smoth"],
     n_total_params);
       */

  // fill_n(g_params , n_total_params , 0);
  // grad_cost_rads(strands);
  // resta_alfav(g_params_final, g_params,
  // -scale_map["rad_fix"],n_total_params);

  // fill_n(g_params , n_total_params , 0);
  // g_CostFunc_compression( strands , mass_gravity , trim_box   ,
  // g_compress_point  ); resta_alfav(g_params_final, g_params,
  // -scale_map["1gravity"], n_total_params);

  // fill_n(g_params , n_total_params , 0);
  // g_CostFunc_compression( strands , mass_gravity  , trim_box  ,
  // g_compress_box  ); resta_alfav(g_params_final, g_params, -scale_map["box"],
  // n_total_params);
}
