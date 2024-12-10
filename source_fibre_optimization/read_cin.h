#include <iostream>
#include <fstream>
#include <map>
#include <sstream>
#include <string>
#include <vector>

using namespace std;

vector<string> split_string(string s, string del );
vector<string> list_files(string path , string extension);


int ** create_matrix_int(int sizeY , int sizeX);

void clear_matrix(int ** matrix , int sizeY , int sizeX);

int str_dist(string s, string t);


void create_folder( string folder) ; 


class parameters 
{
    public:
        string path;
        string initial_conditions;
        string path_out;

        double scale_rad_smoth ;
        double scale_rad_fix;
        double scale_lenght ;
        double scale_curve  ;
        double scale_main_angle ;

        double scale_overlap_fs ;
        double scale_overlap_gs;

        parameters();

        void read_conf_file(string file);

};





// int main00(int argc , char** argv)
// {
//     if (argc <2 )
//     {
//         cout <<"Missing config file"<<endl;
//         return 0;
//     }

//     parameters my_parameters ; 

//     string config_file = argv[1];
//     cout << " -- " << config_file <<endl;

//     my_parameters.read_conf_file(config_file);

//     cout << my_parameters.path<<endl;
//     return 0;
// }
