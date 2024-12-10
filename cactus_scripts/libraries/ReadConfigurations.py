from itertools import product
from libraries import ascii_art
import os

#
def process_string(value , list_format=False):
    if value.isdigit():
        value = int(value)
        if list_format:
            value = [value]
    elif '.' in value:
        value = float(value)
        if list_format:
            value = [value]
    else:
        value = str(value)

    return value

#dictionary to object
class MyConfig:
    def __init__(self, d=None , list_format = False):
        if d is not None:
            for key, value in d.items():
                # if value is a string of int or float, convert items to that type
                # print(key ,value)
                if type(value) == str:
                    value = process_string(value, list_format = list_format)
                elif type(value) == list:
                    for i in range(len(value)):
                       value[i] = process_string(value[i] , list_format = False) 
                else: 
                    print("Error pocessing this type of value")
                    print(f" It is type {type(value)}")


                setattr(self, key, value)

            if 'number_of_repetitions' in d:
                if list_format:
                    n = self.number_of_repetitions[0]
                else:
                    n = self.number_of_repetitions
                self.number_of_repetitions = list(range(n))
            if 'gamma_theta' in d and 'gamma_kappa' in d:
                # check if gamma_theta is number int or float
                if type(self.gamma_theta) == int or type(self.gamma_theta) == float:
                    self.gamma_theta = [self.gamma_theta]
                    self.gamma_kappa = [self.gamma_kappa]

                if len(self.gamma_theta) != len(self.gamma_kappa):
                    raise ValueError("gamma_theta and gamma_kappa must have the same shape")
                gamma_pairs = zip(self.gamma_theta, self.gamma_kappa)
                self.gamma_pairs = list(gamma_pairs)


    #create print function
    def __str__(self):
        return str(self.__dict__)


def write_my_config(my_obj):
    #iterate over elements of the class 

    # outfile = f"{prefix}_{str(i).zfill(5)}_{i}.txt"
    config_file = my_obj.outfile.split(".")[0] + ".lconf"
    my_dir = dir(my_obj)

    parameters_name , _ =  combination_parameters_iterator(my_obj , return_names = True)


    f = open(config_file , 'w')

    #write items and values of my_obj
    for key ,val in my_obj.__dict__.items():
        a = key
        b = getattr(my_obj, key)
        f.write(f'{a}\t {b}\n')


    """ #write items and values of my_dir """
    """ a = "dispersion"  """
    """ b = my_obj.dispersion """
    """ f.write(f'{a}\t {b}\n') """
    """  """
    """  """
    """ a = "icvf"  """
    """ b = my_obj.icvf """
    """ f.write(f'{a}\t {b}\n') """
    """"""
    """ a = "n_bundles"  """
    """ b = my_obj.n_bundles """
    """ f.write(f'{a}\t {b}\n') """
    """"""
    """ a = "outfile"  """
    """ b = my_obj.outfile """
    """ f.write(f'{a}\t {b}\n') """
    """"""
    """"""
    """ a = "angle_fibers"  """
    """ b = my_obj.angle_fibers """
    """ f.write(f'{a}\t {b}\n') """
    """"""
    """  """
    """ a = "gamma_theta"  """
    """ b = my_obj.gamma_theta """
    """ f.write(f'{a}\t {b}\n') """
    """"""
    """ a = "gamma_kappa"  """
    """ b = my_obj.gamma_kappa """
    """ f.write(f'{a}\t {b}\n') """
    """"""
    """ a = "lenght_side"  """
    """ b = my_obj.lenght_side """
    """ f.write(f'{a}\t {b}\n') """
    
    f.close()


## 

def process_list(random_list):
    print("random_list", random_list)
    if type(random_list) in [float, int]: 
        random_list = [random_list]
        return random_list
    elif type(random_list) == list:
        print("list already")
    else:
        a=1
    return random_list

# read file with two columes: key and value
#ignore lines with # or empty lines
def read_config_file(file_name , list_format = False):
    #check if file exists
    if not os.path.isfile(file_name):
        ascii_art.print_error_message(f"Error!!! File {file_name} does not exist")
        exit()
    d = {}
    with open(file_name) as f:
        for line in f:
            #if line starts with # ignore
            if line.startswith('#') or len(line) == 0:
                continue
            line = [ x for x in line.split() if x != '' ]
            if len(line) == 2:
                key, value = line
            if len(line) > 2:
                key = line[0]
                value = line[1:]
            if len(line) <1:
                continue
                
            # print("line: ", line)
            # print("key " , key)
            # print("value " , value)
            d[key] = value

    return MyConfig(d , list_format = list_format)


def combination_parameters_iterator(cactus_args , return_names = False):

    parameters_to_combine   =  [ 'bias', 'dispersion', 'icvf', "number_of_repetitions", "gamma_pairs"] 
    if return_names:
        return parameters_to_combine , None
    crossing_params = ["overlap", "crossing_angle" , "depth_lenght_bundle"  ] 
    if cactus_args.n_bundles[0] > 1:
        parameters_to_combine = parameters_to_combine + crossing_params

    iterator_combinations  = list(product(*(getattr(cactus_args, key) for key in parameters_to_combine)))

    return parameters_to_combine , iterator_combinations


def get_files_with_extension(extension=".init"):
    """
    Retrieves all files in the current directory with the specified extension.
    """
    return sorted([file for file in os.listdir(".") if file.endswith(extension)])



def count_number_streamlines(file):
    with open(file) as f:
        first_line = f.readline()
        second_line = f.readline()

    n = int(second_line)
    return n 

def read_last_line(file_name):
    """
    Reads the last line of a file and returns the first element split by whitespace.
    It reads the last know cost function value to check if it needs to continue optimizing.
    """
    print(f"Processing file: {file_name}")
    if os.path.isfile(file_name):
        with open(file_name, 'r') as f:
            lines = f.readlines()
        return lines[-1].split()[0]
    else:
        return float('inf')
