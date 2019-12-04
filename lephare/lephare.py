"""
$LEPHAREDIR/source/sedtolib -t S -c zphot.para
$LEPHAREDIR/source/sedtolib -t Q -c zphot.para
$LEPHAREDIR/source/sedtolib -t G -c zphot.para
$LEPHAREDIR/source/filter -c zphot.para
$LEPHAREDIR/source/mag_star -c zphot.para
$LEPHAREDIR/source/mag_gal -t Q -c zphot.para
$LEPHAREDIR/source/mag_gal -t G -c zphot.para
$LEPHAREDIR/source/zphota -c zphot.para
"""

import sys
import subprocess as sp
import os
import copy
import glob
import tempfile

from astropy import units as u
from astropy.table import Table
import numpy as np
import matplotlib.pyplot as mpl
import matplotlib.ticker as mpt

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
            assert isinstance(parameters, list), "parameters must be a list of strings."
            for key, val in LePhare.PARAMETERS.items():
                if key not in parameters:
                    LePhare.PARAMETERS[key] = 0
                else:
                    LePhare.PARAMETERS[key] = 1

        return


    def set_working_directory(self,dir):
        os.environ["LEPHAREWORK"] = dir
        self.workdir = dir
        return

    def get_filters(self,filters,SearchFolder=None):
        lepharedir = os.environ['LEPHAREDIR']
        if SearchFolder is None:
            SearchFolder = "*"

        filter_files = []
        for fltr in filters:
            searchPattern = f"{lepharedir}/filt/{SearchFolder}/*{fltr}*"
            file = glob.glob(searchPattern)
            if len(file) == 0:
                raise lephareError(f"Filter {fltr} file not found in {lepharedir}/filt/{SearchFolder}.")
            filter_files.append( "/".join(file[0].split("/")[-2:]) )

        return filter_files

    def set_filters(self,filters,SearchFolder=None):
        filter_files = self.get_filters(filters,SearchFolder=SearchFolder)
        LePhare.CONFIG["FILTER_LIST"] = ",".join(filter_files)
        return None

    def run_cmd(self,cmd,debug=False):
        result = sp.run(cmd,shell=True,stdout=sp.PIPE)
        if debug is True:
            print(result.stdout.decode())
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

    def zphota(self,parfile,debug=False):
        cmd = f"$LEPHAREDIR/source/zphota -c {parfile}"
        result = self.run_cmd(cmd,debug=debug)
        return result

    def prepare_files(self,parfile):
        work_folders = ["filt", "lib_bin", "lib_mag"]
        for fldr in work_folders:
            path = f"{self.workdir}/{fldr}"
            if not os.path.isdir(path):
                os.mkdir(path)

        f = open(f"{parfile}","w")
        for key, val in LePhare.CONFIG.items():
            f.write(f"{key}\t{val}\n")
        f.close()

        f = open(f"{LePhare.CONFIG['PARA_OUT']}","w")
        for key, val in LePhare.PARAMETERS.items():
            if val == 1:
                f.write(f"{key}\n")
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
        self.zphota(parfile,debug=debug)
        return

    def create_lephare_input(self,table,bands,tmpdir=None,errsuffix="_err",longcat=False,fname=None):
        if tmpdir is None:
            tmpdir = os.getcwd()

        if fname is None:
            fpointer,fname = tempfile.mkstemp(prefix="lephare_tmp",dir=tmpdir)
        else:
            fpointer = open(fname,"w")

        if not "ID" in table.colnames:
            table["ID"] = np.arange(1,len(table)+1,dtype=np.int)
        elif table["ID"].dtype == int:
            pass
        else:
            raise lephareError("ID column must be of integer type.")

        cols = ["ID"]
        for bnd in bands:
            if not bnd in table.colnames:
                raise lephareError(f"{bnd} not found in input table.")

            if not f"{bnd}{errsuffix}" in table.colnames:
                raise lephareError(f"{bnd}{errsuffix} not found in input table. Consider changing errsuffix {errsuffix} to match table names.")

            if table[bnd].unit is None:
                raise lephareError(f"{bnd} column has no units.")
            elif table[bnd].unit is u.mag:
                pass
            else:
                table[bnd].convert_unit_to("erg/(s cm2 Hz)")
                cols.append(bnd)
                table[f"{bnd}{errsuffix}"].convert_unit_to("erg/(s cm2 Hz)")
                cols.append(f"{bnd}{errsuffix}")

                nanSel = np.isnan(table[bnd])
                table[bnd][nanSel] = -99.
                table[f"{bnd}{errsuffix}"][nanSel] = -99.

        if longcat is True:
            table["context"]=0
            cols.append("context")
            cols.append("zphot")

        table[cols].write(fpointer,format="ascii.commented_header")
        return fname

    def config(self):
        print("="*30)
        print("Current configuration parameters")
        print("="*30)
        for key, val in LePhare.CONFIG.items():
            print(f"{key}:\t{val}")
        return

    def param(self):
        print("="*30)
        print("Current output parameters")
        print("="*30)
        for key, val in LePhare.PARAMETERS.items():
            if val == 1 :
                print(f"{key}")
        return


    def _get_colnames(self,filename):

        f = open(filename,"r")
        txt = f.readlines()
        f.close()

        bands = []
        colnames = []
        for i,line in enumerate(txt):
            if not "#" in line[0]:
                break

            elif "Output" in line:
                kOutput = i+1

            elif "FILTER_FILE" in line:
                bands = txt[i+1].split()[2:]

        outputColumns = " ".join(txt[kOutput:i-1])

        columns = []
        for cln in outputColumns.split():
            if len(cln)>1 and not cln.isdigit():
                columns.append(cln)

        colBands = []
        k=0
        if not "MAG_OBS()" in columns:
            colBands += []
        else:
            for coln in ["MOBS_%s","MOBSERR_%s"]:
                for b in bands:
                    colBands += [coln%(b)]

            for colmag in ["MAG_OBS()","ERR_MAG_OBS()"]:
                k = columns.index(colmag)
                columns.pop(k)

        if not "MAG_MOD()" in columns:
            colBands += []
        else:
            for coln in ["MMOD_%s"]:
                for b in bands:
                    colBands += [coln%(b)]

                k = columns.index("MAG_MOD()")
                columns.pop(k)

        if not "MAG_ABS()" in columns:
            colBands += []
        else:
            for coln in ["MABS_%s","MABSERR_%s"]:
                for b in bands:
                    colBands += [coln%(b)]

            for colmag in ["MAG_ABS()","ERR_MAG_ABS()"]:
                k = columns.index(colmag)
                columns.pop(k)

        FinalColumns = columns[:k] + colBands + columns[k:]
        return FinalColumns


    def read_output(self,fname=None):
        if fname is None:
            fname = LePhare.CONFIG["CAT_OUT"]

        tableLePhare = Table.read(fname,format="ascii")
        colnames = self._get_colnames(fname)
        for i,c in enumerate(tableLePhare.colnames):
            tableLePhare[c].name = colnames[i]
        return tableLePhare

    def cleanup(self, parfile):
        files = [ LePhare.CONFIG["CAT_OUT"], LePhare.CONFIG["CAT_IN"],\
                  LePhare.CONFIG["PARA_OUT"], parfile]
        for fl in files:
            os.remove(fl)
        return


def plot_SED(filename, ax = None, sedOnly=False, **kwargs):
    ### Open .spec file[s] and read the parameters
    fsp=open(filename,'r')

    bid=fsp.readline()    #header1
    line=fsp.readline()
    line=line.split()
    id=line[0]; zspec=line[1]; zphot=float(line[2])
    #z68low=float(line[3]); z68hig=float(line[4])

    bid=fsp.readline()    #header2
    line=fsp.readline()
    line=line.split()
    nfilt=int(line[1])

    bid=fsp.readline()    #header3
    line=fsp.readline()
    line=line.split()
    npdf=int(line[1])

    bid=fsp.readline()
    #header4:  Type Nline Model Library Nband    Zphot Zinf Zsup Chi2  PDF     Extlaw EB-V Lir Age  Mass SFR SSFR
    models_info=[]
    for i in range(6):
        line=fsp.readline()
        model=line.split()
        models_info.append(model)

    # Read observed mag, err_mag, filters' lambda_eff and width, models' mag
    mag=np.zeros(nfilt); em=np.zeros(nfilt);
    lf=np.zeros(nfilt); dlf=np.zeros(nfilt);
    mmod=np.zeros(nfilt); mfir=np.zeros(nfilt); mphys=np.zeros(nfilt)
    for i in range(nfilt):
        line=fsp.readline(); line=line.split()
        mag[i]=float(line[0]); em[i]=float(line[1]);
        lf[i]=float(line[2]); dlf[i]=float(line[3]);
        mmod[i]=float(line[4]); mfir[i]=float(line[5]); mphys[i]=float(line[6])

    #convert mag(AB syst.) in log(flux)
    ibad=np.where((mmod<=0) | (mmod>35))
    mmod=-0.4*(mmod-23.91)  # uJy
    mmod[ibad]=-10.
    ibad=np.where((mphys<=0) | (mphys>35))
    mphys=-0.4*(mphys-23.91)  # uJy
    mphys[ibad]=-10.
    ibad=np.where((mfir<=0) | (mfir>35))
    mfir=-0.4*(mfir-23.91)  # uJy
    mfir[ibad]=-10.

    zpdf=np.zeros([2,npdf])
    for i in range(npdf):
        line=fsp.readline()
        zpdf[:,i]=np.array(line.split())

    # Read spectra [lambda(A), Mag(AB)]
    # convert in log(F_nu) = -0.4*(mab-23.91) [uJy]
    # Loop over the 6 models (GAL-1, GAL-2, GAL-FIR, GAL-STOCH, QSO, STAR)
    lg=[]; mg=[]
    for m in range(6):
        nline=int(models_info[m][1])
        bid=np.zeros([2,nline])
        if nline>0:
            for i in range(nline):
                line=fsp.readline()
                bid[:,i]=np.array(line.split())
                if (bid[1,i]>35):
                    bid[1,i]=-10.
                else:
                    bid[1,i]=-0.4*(bid[1,i]-23.91)
        lg.append(bid[0,:]/10000.)
        mg.append(bid[1,:])

    fsp.close()

    ##############  PLOT  ###############



    ### Main panel
    if ax is None:
        ### Initialise the figure
        fig=mpl.figure(figsize=(10,8))
        ax1=fig.add_axes([.1,.1,.78,.78], xscale='log',\
        xlabel='$\lambda$ [$\mu$m]',ylabel='log(F$_{\\nu}$ [$\mu$Jy]) ')
    else:
        fig = None
        ax1 = ax


    # only the reliable obs mag will be plotted:
    em=em*2.
    dlf=dlf/2.
    mag1=mag[(mag>0.) & (mag<35) & (em>-3)]
    em1=em[(mag>0.) & (mag<35) & (em>-3)]
    lf1=lf[(mag>0.) & (mag<35) & (em>-3)]/10000.
    dlf1=dlf[(mag>0.) & (mag<35) & (em>-3)]/10000.

    ymin=max(mag1+2.); ymax=min(mag1-4.)
    if ymin>60:
        ymin=30

    ic=[(em1>=0.) & (em1<2.)]
    lf2=lf1[ic]
    mag2=-.4*(mag1[ic]-23.91)
    em2=0.4*em1[ic]
    dlf2=dlf1[ic]
    # low S/N bands:
    ic2=[(em1>=2.) & (em1<8.)]
    lf2b=lf1[ic2]
    mag2b=-.4*(mag1[ic2]-23.91)
    em2b=0.4*em1[ic2]
    dlf2b=dlf1[ic2]

    # set the plot aspect
    ax1.xaxis.set_major_formatter(mpt.ScalarFormatter())
    ax1.set_xticks([0,1,0.2,0.3,0,4,0.5,0.7,1,2,3,4,5,7,10])
    ax1.axis([0.1,max(lf1)*1.2,-0.4*(ymin-23.91),-0.4*(ymax-23.91)])

    ### plot SED and print info of best-fit models
    col_lst=["k",'w','r','g','b','m','y','c']  #each one with a different color
    # plt.figtext(0.15,0.96,' Type: (Model Library Nband) z_phot  Chi^2,  Extlaw  EB-V  Lir  Age  logM*  logSFR',size='small')
    # plt.figtext(0.73,0.12,'ID='+id,size='small')
    iml=0
    for im in range(6):
        if int(models_info[im][2])<0:
          continue   #print only models really used
        iml=iml+1  #counter of models used
        ax1.plot(lg[im],mg[im], **kwargs)  #plot the SED
        del models_info[im][6:8] #do not print z_inf and z_sup
        del models_info[im][-1]  #nor sSFR
        info1=('  '.join(['%.3f']*len(models_info[im][5:7]))) % tuple([float(j) for j in models_info[im][5:7]])
        if float(models_info[im][8])>=0.:    #additional information
            info2=('   '.join(['%.2f']*len(models_info[im][8:]))) % tuple([float(j) for j in models_info[im][8:]])
            info2=',  '+info2+'.'
        else:
            info2='.'
        infol=models_info[im][0] + ': ('+' '.join(models_info[im][2:5])+')  ' + info1  + info2
        # plt.figtext(0.15,0.96-0.02*iml,infol,color=col_lst[im],size='x-small') #print the rest

    # plot the obs mag...
    if sedOnly is True:
        pass
    else:
        if "color" in kwargs.keys():
            color=kwargs["color"]
        else:
            color = "k"

    ax1.errorbar(lf2b,mag2b,yerr=em2b,xerr=dlf2b,fmt='o',color="crimson",mfc=color,ms=10,mew=2)
    ax1.errorbar(lf2,mag2,yerr=em2,xerr=dlf2,fmt='o',color="ForestGreen",mfc=color,ms=10,mew=2)
    #... and upper limits
    iu=np.where(em1<0)
    if len(iu[0])>0 :
        lf3=lf1[iu]
        mag3=-0.4*(mag1[iu]-23.91)
        ax1.quiver(lf3,mag3,0,-1,units='height',width=0.004,headwidth=5,color='k',pivot='tip')

    ### 2nd panel (inset) showing PDF(z)
    base=0.9-0.02*iml #starting position for the inset plot
    if base>0.84:
        base=0.84
    # ax2=fig.add_axes([0.13,base-0.20,0.3,0.20],xlabel='z_phot',title='z_spec='+zspec)
    # ax2.yaxis.set_label_position("right")
    # ax2.yaxis.set_ticks_position("right")
    # ax2.plot(zpdf[0,:],zpdf[1,:],color='r')
    #plot also z_phot with error bar
    #ax2.errorbar(zphot,0.5,fmt='ok',xerr=[[zphot-z68low],[z68hig-zphot]],mfc='none')

    # no chose  window/png for the moment TBI
    #plt.show()
    return fig,ax1


if __name__ == "__main__" :

    run1 = LePhare(config = {"ZFIX":"YES"})
    # run1.config()
    # run1.param()
    run1.set_working_directory("/Users/bribeiro/git-projects/pyLePhare")
    print(os.environ.get("LEPHAREWORK"))
