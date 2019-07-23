"""
$LEPHAREDIR/source/sedtolib -t S -c nnnnn.para
$LEPHAREDIR/source/sedtolib -t Q -c nnnnn.para
$LEPHAREDIR/source/sedtolib -t G -c nnnnn.para
$LEPHAREDIR/source/filter -c nnnnn.para
$LEPHAREDIR/source/mag_star -c nnnnn.para
$LEPHAREDIR/source/mag_gal -t Q -c nnnnn.para
$LEPHAREDIR/source/mag_gal -t G -c nnnnn.para
$LEPHAREDIR/source/zphota -c nnnnn.para
"""

import sys
import subprocess as sp
import os
import copy

from .default import default_pars, default_config

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



class LePhare:

    CONFIG = copy.deepcopy(default_config)
    PARAMETERS = copy.deepcopy(default_pars)

    def __init__(self, config = None, parameters = None):


        if config is not None:
            assert isinstance(config, dict), "config must be a dictionary."
            for key, val in config.items():
                if key not in LePhare.CONFIG.keys():
                    raise lephareError(f"{key} not valid in config file")
                else:
                    LePhare.CONFIG[key] = val

        if parameters is not None:
            assert isinstance(parameters, dict), "parameters must be a dictionary."
            for key, val in parameters.items():
                if key not in LePhare.PARAMETERS.keys():
                    raise lephareError(f"{key} not valid in parameters file")
                else:
                    LePhare.PARAMETERS[key] = val
        return


    def set_working_directory(self,dir):
        os.environ["LEPHAREWORK"] = dir
        self.workdir = dir
        return

    def run_cmd(self,cmd,debug=False):
        result = sp.run(cmd,shell=True,stdout=sp.PIPE)
        if debug is True:
            print(result.stdout)
        return result

    def sedtolib(self,lib,parfile,debug=False):
        cmd = f"$LEPHAREDIR/source/sedtolib -t {lib} -c {parfile}"
        result = self.run_cmd(cmd,debug=debug)
        return result

    def filter(self,parfile,debug=False):
        cmd = f"$LEPHAREDIR/source/filter -c {parfile}"
        result = self.run_cmd(cmd,debug=debug)
        return result

    def mag(self,lib,parfile,debug=False):
        if lib == "S":
            cmd = f"$LEPHAREDIR/source/mag_star -c {parfile}"
        else:
            cmd = f"$LEPHAREDIR/source/mag_gal -t {lib} -c {parfile}"
        result = self.run_cmd(cmd,debug=debug)
        return result

    def zphota(self,parfile,deubg=False):
        cmd = f"$LEPHAREDIR/source/zphota -c {parfile}"
        result = self.run_cmd(cmd,debug=debug)
        return result

    def prepare_files(self,parfile):

        if not os.path.isdir(workdir):
            os.mkdir(workdir)

        f = open(f"{self.workdir}/{parfile}","w")
        for key, val in LePhare.CONFIG.items():
            f.write(f"{key}\t{val}\n")
        f.close()

        f = open(f"{LePhare.CONFIG['PARA_OUT']}","w")
        for key, val in LePhare.PARAMETERS.items():
            f.write(f"{key}\t{val}\n")
        f.close()

        return

    def run_all(self,parfile,debug=False):
        self.prepare_files(parfile)

        self.sedtolib("S",parfile,debug=debug)
        self.sedtolib("Q",parfile,debug=debug)
        self.sedtolib("G",parfile,debug=debug)

        self.filter(parfile,debug=debug)

        self.mag("S",parfile,debug=debug)
        self.mag("Q",parfile,debug=debug)
        self.mag("G",parfile,debug=debug)

        self.zphota(parfile,debug=debug)
        return

    def run_gal(self,parfile,debug=False):
        self.prepare_files(parfile)
        self.sedtolib("G",parfile,debug=debug)
        self.filter(parfile,debug=debug)
        self.mag("G",parfile,debug=debug)
        zself.phota(parfile,debug=debug)
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