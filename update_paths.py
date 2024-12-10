import os

def update_cactus_code_line(file_path):
    """
    Update the 'cactus_code' line in the given Python file with the current directory.
    
    Parameters:
        file_path (str): Path to the CactusPath.py file.
    """
    current_directory = os.getcwd()
    new_cactus_code_line = f'    "cactus_code": str("{current_directory}/cactus_scripts"),\n'

    try:
        # Read the file and modify the 'cactus_code' line
        with open(file_path, 'r') as file:
            lines = file.readlines()

        # Find and update the 'cactus_code' line
        for i, line in enumerate(lines):
            if line.strip().startswith('"cactus_code":'):
                lines[i] = new_cactus_code_line
                break

        # Write the updated content back to the file
        with open(file_path, 'w') as file:
            file.writelines(lines)

        print(f"'cactus_code' updated successfully in {file_path}.")

    except FileNotFoundError:
        print(f"Error: File {file_path} not found.")
    except Exception as e:
        print(f"An error occurred: {e}")

def update_optimizer_path(file_path):
    """
    Update the 'cactus_code' line in the given Python file with the current directory.
    
    Parameters:
        file_path (str): Path to the CactusPath.py file.
    """
    current_directory = os.getcwd()
    new_cactus_code_line = f'    "global_optimizer": str("{current_directory}/cactus_scripts/joint_fibre_optimizer.out)",\n'

    try:
        # Read the file and modify the 'cactus_code' line
        with open(file_path, 'r') as file:
            lines = file.readlines()

        # Find and update the 'cactus_code' line
        for i, line in enumerate(lines):
            if line.strip().startswith('"global_optimizer":'):
                lines[i] = new_cactus_code_line
                break

        # Write the updated content back to the file
        with open(file_path, 'w') as file:
            file.writelines(lines)

        print(f"'optimizer path' updated successfully in {file_path}.")

    except FileNotFoundError:
        print(f"Error: File {file_path} not found.")
    except Exception as e:
        print(f"An error occurred: {e}")

if __name__ == "__main__":

    print("Initializing the CactusPaths.py file...")
# Specify the path to CactusPath.py
    cactus_path_file = "cactus_scripts/libraries/CactusPaths.py"

# Update the 'cactus_code' line
    update_cactus_code_line(cactus_path_file)
    print("CactusPaths.py file initialized successfully.")

    #optimizer_path = "
    update_optimizer_path(cactus_path_file)

