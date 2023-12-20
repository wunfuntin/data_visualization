# Package imports
import os
import pyrosetta
import shutil
import subprocess

from datetime import datetime
import os

def add_directory_to_path(directory):
    """Adds the specified directory to the PATH environment variable."""
    os.environ['PATH'] += os.pathsep + directory

def get_top_pdbs_from_silent(pdb_file_list, silent_file_path, outputDir):


    current_time = datetime.now().strftime('%m%d%y_%H_%M')

    directory_name = (outputDir +
                      f'/pdbs_{current_time}')

    if not os.path.exists(directory_name):
        os.makedirs(directory_name, exist_ok=True)

    # Initialize pyrosetta
    pyrosetta.init()

    file_paths = []
    for pdb in pdb_file_list:

        command = f"silentextractspecific {silent_file_path} {pdb}"


        try:
            subprocess.run(command,
                           shell=True,
                           check=True,
                           capture_output=True,
                           text=True,
                           timeout=30,
                           )

        except subprocess.TimeoutExpired as e:
            return f' Time out error occured as: {e}'

        with open('output.log', 'r') as file:
            file.read()

        output_file_name = f'{pdb}.pdb'
        source_file_path = os.getcwd() + "/" + output_file_name
        destination_file_path = os.path.join(directory_name, output_file_name)
        shutil.move(source_file_path, destination_file_path)

        file_paths.append(os.path.join(directory_name, output_file_name))

    return file_paths
