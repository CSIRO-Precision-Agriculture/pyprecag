import json
import os
from collections import OrderedDict

# Read in Config File and use as module constant
CONFIG_FILE = os.path.join(os.path.abspath(os.path.dirname(__file__)), 'config.json')

def get_debug_mode():
    """
    Gets the value of debug_mode for the module
    """
    try:
        return debug_mode
    except NameError:
        set_debug_mode(False)
        return debug_mode

def set_debug_mode(debug):
    """
    Sets the value of debug_mode for the module constant
     Args:
        debug (bool): True or False value, determines if pyprecag is in debug mode
    """
    if not isinstance(debug, bool):
        raise TypeError("debug must be a boolean")

    global debug_mode
    debug_mode = debug

def read_config():
    """
    Reads the JSON Configuration return and returns a dictionary of values
    Returns:
        collections.OrderedDict : Dictionary of configuration settings.
    """

    if not os.path.exists(CONFIG_FILE):
        raise IOError("Invalid path: {}".format(CONFIG_FILE))

    with open(CONFIG_FILE) as f:
        try:
            # use an ordered dict to maintain the order of the json file.
            config_dict = json.load(f, object_pairs_hook=OrderedDict)

        except ValueError:
            message = 'The file {} does not appear to be valid JSON'.format(CONFIG_FILE)
            raise ValueError(message)

    return config_dict


def write_config(json_dict):
    """
    Writes A Json Configuration Dictionary to file.
     Args:
        json_dict (Dict): The configuration dictionary to write to file

    """

    with open(CONFIG_FILE, 'w') as f:
        json.dump(json_dict, f, indent=4)


def get_config_key(key):
    """
    Get an existing configuration key from the json configuration file.
    Currently getting a key from a nested dictionary is not supported.
    ie ['geoCSV']['xCoordinate_ColumnName']

    Args:
        key (str): The configuration key to extract value for
    Returns: The value of the Key
    """
    config_dict = read_config()

    try:
        val = config_dict[key]

    except KeyError as e:
        raise KeyError('Could not find config key')

    return val


def set_config_key(key, new_value):
    """
    Apply a value to and existing key in the configuration json file and write to file.

    Currently applying a key to a nested dictionary is not supported.
    ie ['geoCSV']['xCoordinate_ColumnName']

    Args:
        key (str): The key to apply a new value to.
        new_value (str): The new value for the required key
    Returns:
        config_dict: The updated configuration dictionary

    """
    config_dict = read_config()

    try:
        old_value = config_dict[key]

    except KeyError as e:
        raise KeyError('Could not find config key')

    if old_value != new_value:
        config_dict[key] = new_value
        write_config(config_dict)
