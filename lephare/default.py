
DEFAULT_PARS ="""IDENT
Z_BEST
Z_BEST68_LOW
Z_BEST68_HIGH
CHI_BEST
MOD_BEST
EBV_BEST
AGE_BEST
AGE_INF
AGE_MED
AGE_SUP
LDUST_BEST
LDUST_INF
LDUST_MED
LDUST_SUP
MASS_BEST
MASS_INF
MASS_MED
MASS_SUP
SFR_BEST
SFR_INF
SFR_MED
SFR_SUP
SSFR_BEST
SSFR_INF
SSFR_MED
SSFR_SUP
LUM_NUV_BEST
LUM_R_BEST
LUM_K_BEST
MAG_OBS()
ERR_MAG_OBS()
MAG_ABS()
ERR_MAG_ABS()"""


def output_pars():
    f = open("lephare-output.para")
    txt = f.readlines()
    f.close()

    pars = {}
    for line in txt:
        if "#" not in line[0]:
            words = line.split()
            if len(words)==0:
                continue
            if words[0] in DEFAULT_PARS.split():
                pars[words[0]] = 1
            else:
                pars[words[0]] = 0
    return pars


default_pars = output_pars()


def config_pars():
    f = open("lephare.para")
    txt = f.readlines()
    f.close()

    config = {}
    for line in txt:
        if "#" not in line[0]:
            words = line.split()
            if len(words)==0:
                continue
            config[words[0]] = words[1]

    return config

default_config = config_pars()
