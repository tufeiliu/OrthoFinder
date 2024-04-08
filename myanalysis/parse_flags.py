from importlib_resources import files, as_file
import fnmatch
from orthofinder import configfiles
try:
    from configparser import ConfigParser
except ImportError:
    from ConfigParser import ConfigParser  # ver. < 3.0
from typing import Optional, Sequence
import json
import pathlib 

def get_flag_name(input_flag: Optional[str] = None) -> Optional[str]:
    
    for file in files(configfiles).iterdir():
        if file.name == 'usage_options.ini':
            flag_options = ConfigParser()
            flag_options.read(file)
            for option in flag_options.sections():
                for key, val in  flag_options[option].items():
                    if input_flag in json.loads(val):
                        return key

def get_example_data_dir(input_flag: Optional[str] = None, 
                        file_format: Sequence[str] = ["fa", "faa", "fasta", "fas", "pep"]) \
                        -> dict[str, pathlib.PosixPath]:

    flag_name = get_flag_name(input_flag)
    if not flag_name:
        raise Exception("Invalue options flag!")

    data_dir_dict = {}
    if flag_name == "fasta":
        for data_dir in files(example_data).iterdir():
            if data_dir.name in file_format:
                if data_dir.name not in data_dir_dict:
                    data_dir_dict[data_dir.name] = data_dir
    else:
        print("Not an example data type")
    
    return data_dir_dict

    
    

    
    