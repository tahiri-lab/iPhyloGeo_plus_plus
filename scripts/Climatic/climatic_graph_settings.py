import os
import sys
import yaml
from yaml.loader import SafeLoader
from main import find_utils


class ClimaticGraphSettings:
    """
    Class that contains the parameters of the program.
    Loads the parameters from a yaml file or from a dictionary.
    """

    PARAMETER_KEYS = {
        "label_color",
        "edge_color",
        "reticulation_color",
        "layout",
        "proportional_edge_lengths",
        "label_internal_vertices",
        "use_leaf_names",
        "show_branch_length",
        "view_type",
    }

    label_color = "black"
    edge_color = "blue"
    reticulation_color = "red"
    layout = "horizontal"
    proportional_edge_lengths = False
    label_internal_vertices = False
    use_leaf_names = True
    show_branch_length = False
    view_type = "network"

    @classmethod
    def load_from_file(cls, params_file):
        """
        Method that loads the parameters from a yaml file.

        args:
            params_file (str): the name of the yaml file, which must be in the utils directory (extension included)
        """
        
        params_file_with_path = find_utils(params_file)
        
        with open(params_file_with_path) as f:
            params = yaml.load(f, Loader=SafeLoader)
            cls.validate_and_set_params(params)

    @classmethod
    def update_from_dict(cls, params_content):
        """
        Method that updates the parameters from a dictionary.

        args:
            params_content (dict): the dictionary with the parameters
        """
        cls.validate_and_set_params(params_content)

    @classmethod
    def validate_and_set_params(cls, params_dict):
        """
        Method that validates and sets the parameters.

        args:
            params_dict (dict): the dictionary with the parameters
        """

        for key, value in params_dict.items():
            if key in cls.PARAMETER_KEYS:
                setattr(cls, key, value)
            else:
                raise ValueError(f"Invalid parameter: {key}")

    @classmethod
    def get_params(cls):
        return {key: getattr(cls, key) for key in cls.PARAMETER_KEYS}
