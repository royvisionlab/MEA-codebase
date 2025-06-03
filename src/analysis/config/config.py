import platform
import os

def get_save_path():
    if (platform.system() == 'Linux'):
        if ('mmfs1' in os.getcwd()) or ('gscratch' in os.getcwd()): # Hyak
            SAVE_PATH = '/gscratch/retina/data/sorted/'
        else: # Local
            import sys
            sys.path.append(os.path.expanduser('~/Documents/'))
            import custom_config
            if hasattr(custom_config, 'SAVE_PATH'):
                SAVE_PATH = custom_config.SAVE_PATH
            else:
                SAVE_PATH = custom_config.SORT_PATH
    return SAVE_PATH

def get_data_paths():
    # Set the path to the json files according to OS.
    if (platform.system() == 'Windows'):
        JSON_PATH = ''
        SORT_PATH = ''
    elif (platform.system() == 'Linux'):
        if ('mmfs1' in os.getcwd()) or ('gscratch' in os.getcwd()): # Hyak
            JSON_PATH = '/gscratch/retina/data/metadata/json/'
            SORT_PATH = '/gscratch/retina/data/sorted/'
        else: # Local
            import sys
            sys.path.append(os.path.expanduser('~/Documents/'))
            import custom_config
            JSON_PATH = custom_config.JSON_PATH
            SORT_PATH = custom_config.SORT_PATH
    elif (platform.system() == 'Darwin'):
        # Look for a custom path definition file.
        if os.path.exists(os.path.expanduser('~/Documents/custom_config.py')):
            import sys
            sys.path.append(os.path.expanduser('~/Documents/'))
            import custom_config
            JSON_PATH = custom_config.JSON_PATH
            SORT_PATH = custom_config.SORT_PATH
        elif os.path.exists('/Volumes/data/'): # Look for the NAS array mounted to the computer
            JSON_PATH = '/Volumes/data/data/metadata/json/'
            SORT_PATH = '/Volumes/data/data/sorted/'
        else:
            # Warn the user that the data drive is not connected.
            print('Warning!! Cannot find path to data files!')
            JSON_PATH = ''
            SORT_PATH = ''          
    else:
        # Warn the user that the data drive is not connected.
        print('Warning!! Python is not detecting the operating system!')
        JSON_PATH = ''
        SORT_PATH = ''
    return SORT_PATH, JSON_PATH