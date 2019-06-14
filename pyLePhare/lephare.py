import sys
import subprocess as sp
import os
import copy

from default import default_pars, default_config

class lephareException(Exception):
    """Base class for exceptions in this module."""
    pass

class lephareError(lephareException):
    """Exception raised for general errors.

    Attributes:
        message -- explanation of the error
    """

    def __init__(self, message):
        self.message = message



class lePhare:

    CONFIG = copy.deepcopy(default_config)
    PARAMETERS = copy.deepcopy(default_pars)

    def __init__(self, config = None, parameters = None):


        if config is not None:
            assert isinstance(config, dict), "config must be a dictionary."
            for key, val in config.items():
                if key not in lePhare.CONFIG.keys():
                    raise lephareError(f"{key} not valid in config file")
                else:
                    lePhare.CONFIG[key] = val

        if parameters is not None:
            assert isinstance(parameters, dict), "parameters must be a dictionary."
            for key, val in parameters.items():
                if key not in lePhare.PARAMETERS.keys():
                    raise lephareError(f"{key} not valid in parameters file")
                else:
                    lePhare.PARAMETERS[key] = val
        return


    def set_working_directory(self,dir):
        os.environ["LEPHAREWORK"] = dir
        return

    def config(self):
        print("="*30)
        print("Current configuration parameters")
        print("="*30)
        for key, val in lePhare.CONFIG.items():
            print(f"{key}:\t{val}")
        return

    def param(self):
        print("="*30)
        print("Current output parameters")
        print("="*30)
        for key, val in lePhare.PARAMETERS.items():
            if val == 1 :
                print(f"{key}")
        return


if __name__ == "__main__" :

    run1 = lePhare(config = {"ZFIX":"YES"})
    # run1.config()
    # run1.param()
    run1.set_working_directory("/Users/bribeiro/git-projects/pyLePhare")
    print(os.environ.get("LEPHAREWORK"))
