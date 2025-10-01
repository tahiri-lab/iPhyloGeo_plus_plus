import os

import yaml
from aphylogeo.params import Params
from main import find_utils


class MyDumper(yaml.Dumper):
    """
    Custom YAML Dumper to modify the default indentation and list representation behavior.

    Methods:
        increase_indent(flow=False, indentless=False):
            Increase the indentation level in the YAML output.

        represent_list(data):
            Represent Python lists in a flow style in the YAML output.
    """

    def increase_indent(self, flow=False, indentless=False):
        """
        Increase the indentation level in the YAML output.

        Args:
            flow (bool): Indicates whether the current context is a flow style. Defaults to False.
            indentless (bool): Indicates whether to use an indentless format. This argument is ignored. Defaults to False.

        Returns:
            The result from the superclass's increase_indent method with modified behavior.
        """
        return super(MyDumper, self).increase_indent(flow, False)

    def represent_list(self, data):
        """
        Represent Python lists in a flow style in the YAML output.

        Args:
            data (list): The list to represent in the YAML output.

        Returns:
            The YAML representation of the list in a flow style.
        """
        return self.represent_sequence("tag:yaml.org,2002:seq", data, flow_style=True)


yaml.add_representer(list, MyDumper.represent_list, Dumper=MyDumper)
"""
Add a representer for the list type to the YAML dumper.

This ensures that lists are represented in a flow style using the MyDumper class.

Args:
    list (type): The Python list type to represent.
    MyDumper.represent_list (method): The method that defines how to represent lists.
    Dumper (yaml.Dumper): The custom dumper class to use, in this case, MyDumper.
"""


def update_yaml_param(params, file_path, property_name, new_value):
    """
    Updates a specified property within a YAML file with a new value.

    Args:
        params: An object with an update_from_dict method, typically used for updating parameters.
        file_path (str): The path to the YAML file.
        property_name (str): The name of the property to modify (e.g., 'file_name').
        new_value: The new value to set for the property (can be any valid YAML type).

    Raises:
        FileNotFoundError: If the specified YAML file does not exist.
        KeyError: If the specified property name is not found in the YAML file.
    """
    if isinstance(new_value, list):
        new_value = [element.strip() for element in new_value]
    params.update_from_dict({property_name: new_value})

    yaml_file_with_path = find_utils(file_path)
    # 1. Load existing YAML data
    try:
        with open(yaml_file_with_path, "r") as yaml_file:
            data = yaml.safe_load(yaml_file)  # Use safe_load for security
    except FileNotFoundError:
        if isinstance(params, Params):
            data = params.PARAMETER_KEYS
        else:
            data = params.get_params()

    # 2. Update the specified property
    if property_name in data:
        data[property_name] = new_value
    # 3. Write the updated data back to the file
    os.makedirs(os.path.dirname(yaml_file_with_path), exist_ok=True)  # Ensure the file exists

    with open(yaml_file_with_path, "w") as yaml_file:
        yaml.dump(data, yaml_file, default_flow_style=None, Dumper=MyDumper, sort_keys=False)
