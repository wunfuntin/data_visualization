# Package imports
import csv
import os
import pyrosetta
import shutil
import subprocess


from datetime import datetime
from pymol import cmd


def add_directory_to_path(directory):
    """Adds the specified directory to the PATH environment variable."""
    os.environ["PATH"] += os.pathsep + directory


def get_top_pdbs_from_silent(pdb_file_list):
    binary_path = "../silent_tools"
    silent_file_path = (
        "/Users/rwalker/Documents/TT_lab/python/dash_prot_char_v1_1/"
        "silent_test_data/gdf8_denovo_IAP.silent"
    )

    add_directory_to_path(binary_path)

    current_time = datetime.now().strftime("%m%d%y_%H_%M")

    directory_name = (
        f"/Users/rwalker/Documents/TT_lab/python/dash_prot_char_v1_1/silent_test_data/extracted_pdbs/"
        f"pdbs_{current_time}"
    )

    if not os.path.exists(directory_name):
        os.makedirs(directory_name, exist_ok=True)

    # Initialize pyrosetta
    pyrosetta.init()

    file_paths = []
    for pdb in pdb_file_list:
        command = f"{binary_path}/silentextractspecific {silent_file_path} {pdb}"

        try:
            print("subprocess being called...")

            subprocess.run(
                command,
                shell=True,
                check=True,
                capture_output=True,
                text=True,
            )

        except Exception as e:
            return f" Time out error occurred as: {e}"

        with open("output.log", "r") as file:
            file.read()

        output_file_name = f"{pdb}.pdb"
        source_file_path = (
            "/Users/rwalker/Documents/TT_lab/python/dash_prot_char_v1_1/silent_tools/"
            + output_file_name
        )
        destination_file_path = os.path.join(directory_name, output_file_name)
        shutil.move(source_file_path, destination_file_path)

        file_paths.append(os.path.join(directory_name, output_file_name))

    return file_paths


def extract_pdb_sequence(pdb_file_path):
    """Extracts the sequence of chain 'A' from the given pdb file."""

    cmd.load(pdb_file_path, "temp_structure")
    fasta_str = cmd.get_fastastr("temp_structure and chain A")
    cmd.delete("temp_structure")

    # Extract the sequence from the fasta string
    sequence = "".join(fasta_str.split("\n")[1:])
    return sequence


def get_pdb_sequences(pdb_directory):
    pdb_directory = pdb_directory
    output_csv = os.path.join(pdb_directory, "sequences.csv")

    with open(output_csv, "w", newline="") as csvfile:
        csvwriter = csv.writer(csvfile)
        csvwriter.writerow(["Filename", "Sequence"])

        for filename in os.listdir(pdb_directory):
            if filename.endswith(".pdb"):
                file_path = os.path.join(pdb_directory, filename)
                sequence = extract_pdb_sequence(file_path)
                csvwriter.writerow([filename, sequence])

    cmd.quit()
