# Package imports
import os

from dash_bio.utils import PdbParser, create_mol3d_style


def get_most_recent_directory(parent_directory):
    parent_directory = parent_directory

    # List all subdirectories in the parent directory
    subdirectories = [
        d
        for d in os.listdir(parent_directory)
        if os.path.isdir(os.path.join(parent_directory, d))
    ]

    # Sort subdirectories by creation time (newest first)
    subdirectories.sort(
        key=lambda d: os.path.getctime(os.path.join(parent_directory, d)), reverse=True
    )

    # Get the most recently created directory
    if subdirectories:
        most_recent_directory = os.path.join(parent_directory, subdirectories[0])
    else:
        most_recent_directory = None

    return most_recent_directory


def modify_style(selected_file):
    if selected_file:
        try:
            parser = PdbParser(selected_file)
            model_data = parser.mol3d_data()
            styles = create_mol3d_style(
                model_data["atoms"],
                visualization_type="cartoon",
                color_element="chain",
            )

            return styles
        except Exception as e:
            print(f"Error styling PDB file: {e}")
        return None


def list_directories(folder_path):
    # print("Folder path received:", folder_path)  # Debug print
    if not folder_path or not os.path.isdir(folder_path):
        return []
    directories = [
        d
        for d in os.listdir(folder_path)
        if os.path.isdir(os.path.join(folder_path, d))
    ]
    # print("Directories found:", directories)  # Debug print
    return directories
