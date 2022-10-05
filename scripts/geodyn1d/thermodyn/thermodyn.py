
import sys
import numpy as np
import os.path
from scipy import interpolate
import scipy.io as sio
from scipy.io import FortranFile
from copy import deepcopy

from .. import core as geodyn1d
from ..tools import *

from scipy.interpolate import RegularGridInterpolator
from scipy.ndimage.filters import gaussian_filter
from scipy.signal import medfilt2d

class perplex_output():
    def __init__(self):
        self.T           = None
        self.P           = None
        self.nt          = None
        self.np          = None
        self.dP          = None
        self.dT          = None
        self.rho         = None
        self.rhoref      = None
        self.rhoresidual = None
        self.rhomelt     = None
        self.alpha       = None
        self.alpharesidual = None
        self.alphamelt     = None
        self.beta          = None
        self.betaresidual  = None
        self.betamelt      = None
        self.melt          = None
        self.deltarho      = None
        self.composition   = None
        self.rho_ref       = None
        self.rho_ref_std   = None
        self.cp            = None
        self.cpresidual    = None
        self.cpmelt        = None

# columns in perplex table
#                Name              Counter       T(K)                 P(bar)               V,J/bar/mol 
#         H,J/mol              Gruneisen_T          Ks,bar               Gs,bar               v0,km/s   
#         vp,km/s              vs,km/s              vp/vs                rho,kg/m3            G,J/mol    
#         cp,J/K/mol           alpha,1/K            beta,1/bar           S,J/K/mol            n,mol     
#         N,g                  Ks_{T},bar/K         Gs_{T},bar/K         Ks_{P}               Gs_P     
#         v0_{T}               vp_{T}               vs_{T}               v0_{P}               vp_P        
#         vs_{P}               cp/cv                vol,%                wt,%                 mol,%


def read_output_from_perple_x(tablename,path=None,fractionation=False):
    """ Read output from perple_x"""
   
    if ( path == None ):
        path     = os.path.abspath(__file__)
        filename = path.partition('thermodyn.py')[0]+'data/'+tablename+'.phm'
    else:
        filename = path+'/'+tablename+'.phm'

    if not (os.path.isfile(filename)):
        sys.exit("ERROR(read_output_from_perple_x): Either this file does not exist or the path is uncorrect {} ".format(filename))

    if ( fractionation == False ):
        num_lines_header = 13
    else:
        num_lines_header = 9
    
    # The file is automatically closed
    with open(filename, 'r') as f:
        iline = 0
        data                = []
        perplex_composition = []
        ndata               = 0
        perplex             = perplex_output()
        for line in f:
            iline = iline + 1
            line  = line.strip()
            line  = line.replace("Missing data", "Missing_data")
            columns  = line.split()
            if (iline<=num_lines_header):
                if ( fractionation == False ):
                    #print('2d grid input file')
                    # Pmin, Pmax
                    if (iline==9):
                        Pmin = float(columns[0])
                    if (iline==10):
                        perplex.dP = float(columns[0])
                    if (iline==11):
                        perplex.np = int(columns[0])
                    # Tmin, Tmax
                    if (iline==5):
                        Tmin = float(columns[0])
                    if (iline==6):
                        perplex.dT = float(columns[0])
                    if (iline==7):
                        perplex.nt = int(columns[0])
                #else:
                #    # nothing to read in the header
                #    print('Fractionation 1d grid input file')

                # header
                if (iline==num_lines_header):
                    columns_name = columns
                    ncol         = len(columns_name) 
            else:
                source   = {}
                source[columns_name[0]] = columns[0]
                if ((source[columns_name[0]]=='system') | (source[columns_name[0]]=='Missing_data')):
                    ndata = ndata + 1
                    if ( ndata > 1):
                        perplex_composition.append(composition)
                    counter     = 0
                    composition = {}
                    # counter
                    composition[columns_name[1]]  = columns[1]
                else:
                    counter = counter + 1
                    composition[columns_name[0]+str(counter)]  = columns[0]   # Name
                    composition[columns_name[13]+str(counter)] = columns[13]  # rho,kg/m3
                    composition[columns_name[15]+str(counter)] = columns[15]  # cp,J/K/mol
                    composition[columns_name[16]+str(counter)] = columns[16]  # alpha
                    composition[columns_name[17]+str(counter)] = columns[17]  # beta
                    composition[columns_name[19]+str(counter)] = columns[19]  # n,mol
                    composition[columns_name[20]+str(counter)] = columns[20]  # N,g
                    composition[columns_name[32]+str(counter)] = columns[32]  # vol,%
                    composition[columns_name[33]+str(counter)] = columns[33]  # wt,%

                    #composition[columns_name[35]+str(counter)] = columns[35]  # SiO2,wt%
                    #composition[columns_name[36]+str(counter)] = columns[36]  # Al2O3,wt%
                    #composition[columns_name[37]+str(counter)] = columns[37]  # FeO,wt%
                    #composition[columns_name[38]+str(counter)] = columns[38]  # MgO,wt% 
                    #composition[columns_name[39]+str(counter)] = columns[39]  # CaO,wt% 
                    #composition[columns_name[40]+str(counter)] = columns[40]  # Na2O,wt% 
                    #composition[columns_name[41]+str(counter)] = columns[41]  # Cr2O3,wt%
                    
                    #MSiO2  = 60.09
                    #MAl2O3 = 101.96
                    #MFeO   = 71.84
                    #MMgO   = 40.3
                    #MCaO   = 56.08
                    #MNa2O  = 64.6
                    #MCr2O3 = 152
                    ## Molar mass
                    #M = N/n

                for i in range(1,ncol):
                    source[columns_name[i]] = float(columns[i])
                data.append(source)
    perplex_composition.append(composition)
    
    # Name rho,kg/m3 vol,% wt,% == $1, $14, $33, $34
    #awk '{ if ($1=="system") { if (ratio>0) {print ref,density/ratio,Melt ;}  ref=$2 ; ratio=0 ; density=0; Melt=""; }
    #       else { if ($1!="Melt(JH)") { density = density + ($3/100)*$2; ratio = ratio + $3/100 ; }else {Melt="Melt"} } }' caca > caca2
    
    # convert the list of dictionaries into array
    # be carful to extract the full system average and not all phases
    perplex.T           = np.asarray([ dict['T(K)']       for dict in data if ((dict['Name']=='system')|(dict['Name']=='Missing_data'))])
    perplex.P           = np.asarray([ dict['P(bar)']     for dict in data if ((dict['Name']=='system')|(dict['Name']=='Missing_data'))])
    perplex.rho         = np.asarray([ dict['rho,kg/m3']  for dict in data if ((dict['Name']=='system')|(dict['Name']=='Missing_data'))])

    if ( fractionation == True ):
            perplex.np = ndata
            Pmin       = np.min(perplex.P)
            Pmax       = np.max(perplex.P)
            perplex.dP = (Pmax-Pmin)/(ndata-1)
            perplex.nt = perplex.np
            perplex.dT = 9999999
            for i in np.arange(0,ndata):
                perplex.T[i] = 1573.15 + 0.0011765*perplex.P[i]
                if ( i > 0 ):
                    dT = np.abs(perplex.T[i] - perplex.T[i-1])
                    if ( dT < perplex.dT ):
                        perplex.dT = dT
            Tmin       = np.min(perplex.T)
            Tmax       = np.max(perplex.T)

    perplex.cp          = np.asarray([ dict['cp,J/K/mol']  for dict in data if ((dict['Name']=='system')|(dict['Name']=='Missing_data'))])
    #print(perplex.cp)
    n                   = np.asarray([ dict['n,mol']  for dict in data if ((dict['Name']=='system')|(dict['Name']=='Missing_data'))])
    N                   = np.asarray([ dict['N,g']  for dict in data if ((dict['Name']=='system')|(dict['Name']=='Missing_data'))])
    M                   = np.divide(N,n)*1.0e-3
    perplex.cp          = np.divide(perplex.cp,M) # heat capacity should be in J/K/Kg
    #print(perplex.cp)
    #print('read raw data...........')
    #print(perplex.cp)
    
    perplex.alpha       = np.asarray([ dict['alpha,1/K']  for dict in data if ((dict['Name']=='system')|(dict['Name']=='Missing_data'))])
    perplex.beta        = np.asarray([ dict['beta,1/bar'] for dict in data if ((dict['Name']=='system')|(dict['Name']=='Missing_data'))])
    perplex.composition = perplex_composition
    
    # extract melt
    #{'Counter': '4',
    # 'Name1': 'Cpx(JH)',
    # 'vol,%1': '11.5145',
    # 'wt,%1': '11.4392',
    # 'Name2': 'O(JH)',
    # 'vol,%2': '29.7342',
    # 'wt,%2': '30.2045',
    # 'Name3': 'Pl(JH)',
    # 'vol,%3': '1.73598',
    # 'wt,%3': '1.43973',
    # 'Name4': 'Opx(JH)',
    # 'vol,%4': '57.0153',
    # 'wt,%4': '56.9165'}
    
    book_liq_phases = {
    'abL':'albite_liquid',    # NaAlSi3O8
    'anL':'anorthite_liquid', # CaAl2Si2O8
    'diL':'diopside_liquid',  # CaMgSi2O6
    'enL':'enstatite_liquid', # Mg2Si2O6
    'faL':'fayalite_liquid',         # Fe2SiO4
    'fliq':'Fe-liquid_(in_KFMASH)',  # K3Fe0:5Al4Si19:5O47
    'foL':'Forsterite_liquid',       # Mg2SiO4
    'h2oL':'H2O_liquid',             # H2O
    'hliq':'H2O_liquid_(in_KFMASH)', # H2O
    'kspL':'K-feldspar_liquid',      # KAlSi3O8
    'mliq':'Mg_liquid_(in_KFMASH)',  # K3Mg0:5Al4Si19:5O47
    'qL':'Silica_liquid',            # SiO2
    'silL':'Sillimanite_liquid'      # Al2SiO5
    }
    
    perplex.melt = []
    perplex.rhomelt       = []
    perplex.rhoresidual   = []
    perplex.alphamelt     = []
    perplex.alpharesidual = []
    perplex.betamelt      = []
    perplex.betaresidual  = []
    perplex.cpresidual    = []
    perplex.cpmelt        = []
    for i in range(0,len(perplex.composition)):
        counter        = int(perplex.composition[i]['Counter'])
        melt           = 0.0e0
        rhomelt        = 0.0e0
        alphamelt      = 0.0e0
        betamelt       = 0.0e0
        cpmelt         = 0.0e0
        name_melt      = ''
        rhoresidual    = 0.0e0
        alpharesidual  = 0.0e0
        betaresidual   = 0.0e0
        cpresidual     = 0.0e0
        nresidual      = 0.0e0
        Nresidual      = 0.0e0
        Mresidual      = 0.0e0
        ratio_residual = 0.0e0
        for j in range(1,counter+1):
            name = perplex.composition[i]['Name'+str(j)]
            test = name.partition('(')
            value=test[0]
            if ( value == 'Melt'):
                if ( melt != 0 ):
                    print('already melt name {} prev {}'.format(name,name_melt))
                else:
                    name_melt = name
                    melt      = float(perplex.composition[i]['vol,%'+str(j)])/100.0e0
                    rhomelt   = float(perplex.composition[i]['rho,kg/m3'+str(j)])
                    alphamelt = float(perplex.composition[i]['alpha,1/K'+str(j)])
                    betamelt  = float(perplex.composition[i]['beta,1/bar'+str(j)])
                    cpmelt    = float(perplex.composition[i]['cp,J/K/mol'+str(j)])
                    n         = float(perplex.composition[i]['n,mol'+str(j)])
                    N         = float(perplex.composition[i]['N,g'+str(j)])
                    M         = (N/n)*1.0e-3
                    cpmelt    = cpmelt/M # heat capacity should be in J/K/Kg
            else:
                rhoresidual    = rhoresidual    + (float(perplex.composition[i]['vol,%'+str(j)])/100.0e0) * float(perplex.composition[i]['rho,kg/m3'+str(j)])
                alpharesidual  = alpharesidual  + (float(perplex.composition[i]['vol,%'+str(j)])/100.0e0) * float(perplex.composition[i]['alpha,1/K'+str(j)])
                betaresidual   = betaresidual   + (float(perplex.composition[i]['vol,%'+str(j)])/100.0e0) * float(perplex.composition[i]['beta,1/bar'+str(j)])
                
                cpresidual   = cpresidual   + (float(perplex.composition[i]['vol,%'+str(j)])/100.0e0) * float(perplex.composition[i]['cp,J/K/mol'+str(j)])
                nresidual   = nresidual   + (float(perplex.composition[i]['vol,%'+str(j)])/100.0e0) * float(perplex.composition[i]['n,mol'+str(j)])
                Nresidual   = Nresidual   + (float(perplex.composition[i]['vol,%'+str(j)])/100.0e0) * float(perplex.composition[i]['N,g'+str(j)])
                ratio_residual = ratio_residual + float(perplex.composition[i]['vol,%'+str(j)])/100.0e0
            if ( value in book_liq_phases ):
                melt = melt + float(perplex.composition[i]['vol,%'+str(j)])/100.0e0
                print('liquid phase {} melt {} additional melt {}'.format(value,melt,float(perplex.composition[i]['vol,%'+str(j)])/100.0e0))
            
        perplex.melt.append(melt)
        if ( ratio_residual > 0.0e0 ):
            #print('RESIDUAL ratio {} rho {} alpha {} beta {}'.format(ratio_residual,rhoresidual,alpharesidual,betaresidual))
            perplex.rhoresidual.append(rhoresidual/ratio_residual)
            perplex.alpharesidual.append(alpharesidual/ratio_residual)
            perplex.betaresidual.append(betaresidual/ratio_residual)
            Mresidual  = (Nresidual/nresidual)*1.0e-3
            cpresidual = cpresidual/Mresidual # heat capacity should be in J/K/Kg
            perplex.cpresidual.append(cpresidual/ratio_residual)

        else:
            # never with T <1600ºC in our grids
            #print('RESIDUAL ratio(<=0.0) {} rho {} alpha {} beta {}'.format(ratio_residual,rhoresidual,alpharesidual,betaresidual))
            perplex.rhoresidual.append(0.0e0)
            perplex.alpharesidual.append(0.0e0)
            perplex.betaresidual.append(0.0e0)
            perplex.cpresidual.append(0.0e0)
            
        perplex.rhomelt.append(rhomelt)
        perplex.alphamelt.append(alphamelt)
        perplex.betamelt.append(betamelt)
        perplex.cpmelt.append(cpmelt)
    
    perplex.melt          = np.asarray(perplex.melt)
    perplex.rhoresidual   = np.asarray(perplex.rhoresidual)
    perplex.rhomelt       = np.asarray(perplex.rhomelt)
    perplex.alphamelt     = np.asarray(perplex.alphamelt)
    perplex.betamelt      = np.asarray(perplex.betamelt)
    perplex.cpmelt        = np.asarray(perplex.cpmelt)
    perplex.alpharesidual = np.asarray(perplex.alpharesidual)
    perplex.betaresidual  = np.asarray(perplex.betaresidual)
    perplex.cpresidual    = np.asarray(perplex.cpresidual)
    
    #print('read raw data cpresidual ...... (read_output_from_perple_x)')
    #print(perplex.cpresidual)
    
    # convert units
    if ( fractionation == False ):
        num_col_P = 3
    else:
        num_col_P = 2

    if (columns_name[num_col_P].partition('(')[2].partition(')')[0]=='bar'):
        #print('P given in bar')
        units_P   = 1.0e5
    else:
        units_P = 1.0e0
        print('CHECK pressure conversion to Pa.')
        print('Current units is {}'.format(columns_name[num_col_P].partition('(')[2].partition(')')[0]))
    perplex.P    = perplex.P * units_P
    perplex.beta = perplex.beta / units_P
    perplex.betamelt     = perplex.betamelt / units_P
    perplex.betaresidual = perplex.betaresidual / units_P
    
    if ( (ndata != len(perplex.T)) |
         (ndata != len(perplex.P)) |
         (ndata != len(perplex.rho)) |
         (ndata != len(perplex.alpha)) |
         (ndata != len(perplex.beta)) |
         (ndata != len(perplex.cp)) |
         (ndata != len(perplex.melt)) |
         (ndata != len(perplex.rhoresidual)) |
         (ndata != len(perplex.rhomelt)) |
         (ndata != len(perplex.alphamelt)) |
         (ndata != len(perplex.alpharesidual)) |
         (ndata != len(perplex.betamelt)) |
         (ndata != len(perplex.betaresidual) |
         (ndata != len(perplex.cpmelt)) |
         (ndata != len(perplex.cpresidual)) )
       ):
        print('ndata = {} \nT     = {}\nP     = {}\nrho   = {}\nalpha = {}\nbeta  = {}\ncp  = {}\nmelt  = {}\nrhoresidual = {}\nrhomelt = {}\nalphamelt {}\nbetamelt {}\ncpmelt {}\nalpharesidual {}\nbetaresidual {}\ncpresidual {}'.format(
                ndata,
                len(perplex.T),
                len(perplex.P),
                len(perplex.rho),
                len(perplex.alpha),
                len(perplex.beta),
                len(perplex.cp),
                len(perplex.melt),
                len(perplex.rhoresidual),
                len(perplex.rhomelt),
                len(perplex.alphamelt),
                len(perplex.betamelt),
                len(perplex.cpmelt),
                len(perplex.alpharesidual),
                len(perplex.betaresidual),
                len(perplex.cpresidual)
                ))
    
    perplex.dP = perplex.dP*units_P
    print('Pmin {} dP {} Tmin {} dT {} '.format(Pmin*units_P,perplex.dP,Tmin,perplex.dT))
    print('Pmin {} Pmax {} Tmin {} Tmax {} '.format(min(perplex.P),max(perplex.P),min(perplex.T),max(perplex.T)))

    Tref  = 273
    Pref  = 0

    factor1        = 1 + np.multiply(perplex.beta,perplex.P-Pref)
    factor2        = 1 - np.multiply(perplex.alpha,perplex.T-Tref)
    factor         = np.multiply(factor1,factor2)
    perplex.rhoref = np.divide(perplex.rho,factor) 
 
    return perplex

def check_melt_from_perplex(perplex):
    
    rhomelt   = perplex.rhomelt.reshape(  int(perplex.np), int(perplex.nt))
    alphamelt = perplex.alphamelt.reshape(int(perplex.np), int(perplex.nt))
    betamelt  = perplex.betamelt.reshape( int(perplex.np), int(perplex.nt))
    cpmelt    = perplex.cpmelt.reshape( int(perplex.np), int(perplex.nt))
    melt      = perplex.melt.reshape( int(perplex.np), int(perplex.nt))
    
    for j in range(0,int(perplex.np)):
        for i in range(int(perplex.nt)-1,-1,-1):
            if (melt[j,i]>=0.0001):
                if (np.isnan(rhomelt[j,i])|np.isnan(betamelt[j,i])|np.isnan(cpmelt[j,i])|np.isnan(alphamelt[j,i])):
                    meltfixed = False
                    for kj in range(j,-1,-1):
                        for k in range(i,int(perplex.nt)):
                            if ((melt[kj,k]>=0.0001)&(not np.isnan(rhomelt[kj,k]))):
                                rhomelt[j,i]   = rhomelt[kj,k]
                                alphamelt[j,i] = alphamelt[kj,k]
                                betamelt[j,i]  = betamelt[kj,k]
                                cpmelt[j,i]    = cpmelt[kj,k]
                                meltfixed      = True
                                break
                        if (meltfixed):
                            break
                    # if (meltfixed==False):
                    #     print('Melt not fixed j {} i {}'.format(j,i))
                    # else:
                    #     print('FIXED i {} j {} i_used {} j_used {} '.format(i,j,k,kj))
            else:
                melt[j,i]      = 0.0e0
                alphamelt[j,i] = float('nan')
                betamelt[j,i]  = float('nan')
                cpmelt[j,i]    = float('nan')
                rhomelt[j,i]   = float('nan')
    
    perplex.rhomelt   = rhomelt.reshape(perplex.np*perplex.nt)
    perplex.alphamelt = alphamelt.reshape(perplex.np*perplex.nt)
    perplex.betamelt  = betamelt.reshape(perplex.np*perplex.nt)
    perplex.cpmelt    = cpmelt.reshape(perplex.np*perplex.nt)
    
    return perplex

def remove_melt_from_perplex(perplex,melt_percent=-1):
    """ Extrapolate high temperature values to remove melt content using sub-solidus values.
    The assumption is that alpha and beta are constant and temperature-independent at high temperature."""
     
    Tref = 273
    Pref = 0
    
    rho   = perplex.rho.reshape(  int(perplex.np), int(perplex.nt))
    rhoresidual = perplex.rhoresidual.reshape(  int(perplex.np), int(perplex.nt))
    rhomelt     = perplex.rhomelt.reshape(  int(perplex.np), int(perplex.nt))
    T     = perplex.T.reshape(    int(perplex.np), int(perplex.nt))
    P     = perplex.P.reshape(    int(perplex.np), int(perplex.nt))
    alpha = perplex.alpha.reshape(int(perplex.np), int(perplex.nt))
    alpharesidual = perplex.alpharesidual.reshape(int(perplex.np), int(perplex.nt))
    alphamelt = perplex.alphamelt.reshape(int(perplex.np), int(perplex.nt))
    beta  = perplex.beta.reshape( int(perplex.np), int(perplex.nt))
    betaresidual  = perplex.betaresidual.reshape( int(perplex.np), int(perplex.nt))
    betamelt  = perplex.betamelt.reshape( int(perplex.np), int(perplex.nt))
    cp         = perplex.cp.reshape( int(perplex.np), int(perplex.nt))
    cpmelt     = perplex.cpmelt.reshape( int(perplex.np), int(perplex.nt))
    cpresidual = perplex.cpresidual.reshape( int(perplex.np), int(perplex.nt))
    melt  = perplex.melt.reshape( int(perplex.np), int(perplex.nt))
    
    # smoothing alpha and beta along the boundaries to avoid vertical discontinuities not suitable
    n_smooth = 3
    rho_smooth           = []
    rhomelt_smooth       = []
    rhoresidual_smooth   = []
    alpha_smooth         = []
    beta_smooth          = []
    cp_smooth            = []
    alphamelt_smooth     = []
    betamelt_smooth      = []
    cpmelt_smooth        = []
    alpharesidual_smooth = []
    betaresidual_smooth  = []
    cpresidual_smooth    = []
    
    i_smooth = 0
    i_int    = 0
    #alpha_beta_values = False
    for j in range(0,int(perplex.np)):
        if (melt_percent<0):
            are_values = False
            for i in range(int(perplex.nt)-1,-1,-1):
                #print('None T {} P {} melt {}'.format(T[j,i],P[j,i],melt[j,i]))
                if ( melt[j,i] > 0.0e0  ):
                    #print('None T {} P {}'.format(T[j,i],P[j,i]))
                    pass
                else:
                    if (i_smooth<n_smooth):
                        alpha_smooth.append(alpha[j,i])
                        beta_smooth.append(beta[j,i])
                        cp_smooth.append(cp[j,i])
                        cpmelt_smooth.append(cpmelt[j,i])
                        cpresidual_smooth.append(cpresidual[j,i])
                        alphamelt_smooth.append(alphamelt[j,i])
                        betamelt_smooth.append(betamelt[j,i])
                        alpharesidual_smooth.append(alpharesidual[j,i])
                        betaresidual_smooth.append(betaresidual[j,i])
                        rho_smooth.append(rho[j,i])
                        rhomelt_smooth.append(rhomelt[j,i])
                        rhoresidual_smooth.append(rhoresidual[j,i])
                        i_smooth = i_smooth + 1
                    else:
                        alpha_smooth[i_int]         = alpha[j,i]
                        beta_smooth[i_int]          = beta[j,i]
                        cp_smooth[i_int]            = cp[j,i]
                        cpmelt_smooth[i_int]        = cpmelt[j,i]
                        cpresidual_smooth[i_int]    = cpresidual[j,i]
                        alphamelt_smooth[i_int]     = alphamelt[j,i]
                        betamelt_smooth[i_int]      = betamelt[j,i]
                        alpharesidual_smooth[i_int] = alpharesidual[j,i]
                        betaresidual_smooth[i_int]  = betaresidual[j,i]
                        rho_smooth[i_int]           = rho[j,i]
                        rhomelt_smooth[i_int]       = rhomelt[j,i]
                        rhoresidual_smooth[i_int]   = rhoresidual[j,i]
                        i_int = i_int + 1
                        if (i_int>=n_smooth):
                            i_int = 0
                 
                    alpha_used = sum(alpha_smooth)/len(alpha_smooth)
                    beta_used  = sum(beta_smooth)/len(beta_smooth)
                    cp_used    = sum(cp_smooth)/len(cp_smooth)
                    rho_ref    = sum(rho_smooth)/len(rho_smooth) / ( (1+beta_used*(P[j,i]-Pref)) * (1-alpha_used*(T[j,i]-Tref)) )
                 
                    alpha_used_melt = sum(alphamelt_smooth)/len(alphamelt_smooth)
                    beta_used_melt  = sum(betamelt_smooth)/len(betamelt_smooth)
                    cp_used_melt    = sum(cpmelt_smooth)/len(cpmelt_smooth)
                    rho_ref_melt    = sum(rhomelt_smooth)/len(rhomelt_smooth) / ( (1+beta_used_melt*(P[j,i]-Pref)) * (1-alpha_used_melt*(T[j,i]-Tref)) )
                 
                    alpha_used_residual = sum(alpharesidual_smooth)/len(alpharesidual_smooth)
                    beta_used_residual  = sum(betaresidual_smooth)/len(betaresidual_smooth)
                    cp_used_residual    = sum(cpresidual_smooth)/len(cpresidual_smooth)
                    rho_ref_residual    = sum(rhoresidual_smooth)/len(rhoresidual_smooth) / ( (1+beta_used_residual*(P[j,i]-Pref)) * (1-alpha_used_residual*(T[j,i]-Tref)) )

                    #if ( not alpha_beta_values):
                    #    # we use low pressure value for alpha and beta - upper-bound estimation of it then
                    #    alpha_used = alpha[j,i]
                    #    beta_used  = beta[j,i]
                    #    alpha_beta_values = True
                    #rho_ref    = rho[j,i] / ( (1+beta_used*(P[j,i]-Pref)) * (1-alpha_used*(T[j,i]-Tref)) )
                    melt_ref   = 0.0e0
                    are_values = True
                    break
            
            if (are_values):
                for i in range(int(perplex.nt)-1,-1,-1):
                     if ( melt[j,i] > 0.0e0 ):
#                         rho[j,i]   = rho_ref*(1+beta_used*(P[j,i]-Pref))*(1-alpha_used*(T[j,i]-Tref))
                         rho[j,i]   = rho_ref*(1+betaresidual[j,i]*(P[j,i]-Pref))*(1-alpharesidual[j,i]*(T[j,i]-Tref))
                         #alpha[j,i] = alpha_used
                         #beta[j,i]  = beta_used
                         # we do not extrapolate alpha and beta but only rho_ref
                         # we keep alpha and beta from residual in order to keep them P,T dependant  
                         alpha[j,i] = alpharesidual[j,i]
                         beta[j,i]  = betaresidual[j,i]
                         cp[j,i]    = cpresidual[j,i]
                        
                         melt[j,i]      = melt_ref 
                         rhomelt[j,i]   = float('nan')
                         alphamelt[j,i] = float('nan')
                         betamelt[j,i]  = float('nan')
                         cpmelt[j,i]    = float('nan')
                         
                     else:
                         melt[j,i]      = melt_ref
                         rhomelt[j,i]   = float('nan')
                         alphamelt[j,i] = float('nan')
                         betamelt[j,i]  = float('nan')
                         cpmelt[j,i]    = float('nan')
                         break
        else:
            for i in range(int(perplex.nt)-1,-1,-1):
                # print('melt[j,i] {}'.format(melt[j,i]))
                if (melt[j,i]>melt_percent/100.0e0):
                    melt[j,i]  = melt_percent/100.0e0
                    rho[j,i]   = rhoresidual[j,i]*(100.0e0-melt_percent)/100.0e0 + rhomelt[j,i]*melt_percent/100.0e0
                    alpha[j,i] = alpharesidual[j,i]*(100.0e0-melt_percent)/100.0e0 + alphamelt[j,i]*melt_percent/100.0e0
                    beta[j,i]  = betaresidual[j,i]*(100.0e0-melt_percent)/100.0e0  + betamelt[j,i]*melt_percent/100.0e0
                    cp[j,i]  = cpresidual[j,i]*(100.0e0-melt_percent)/100.0e0  + cpmelt[j,i]*melt_percent/100.0e0
                if (np.isnan(rho[j,i])):
                        print('NaN melt {} rho {} rhoresidual {} rhomelt {} alpha {} beta {}'.format(
                                melt[j,i],rho[j,i],rhoresidual[j,i], rhomelt[j,i], alpha[j,i], beta[j,i]))
                        quit()
                    
                    
    perplex.rho   = rho.reshape(perplex.np*perplex.nt)
    perplex.T     = T.reshape(perplex.np*perplex.nt)
    perplex.P     = P.reshape(perplex.np*perplex.nt)

    perplex.alpha = alpha.reshape(perplex.np*perplex.nt)
    perplex.beta  = beta.reshape(perplex.np*perplex.nt)
    perplex.cp    = cp.reshape(perplex.np*perplex.nt)

    perplex.melt  = melt.reshape(perplex.np*perplex.nt)
    perplex.melt  = np.zeros_like(perplex.melt)
    
    perplex.rhomelt   = rhomelt.reshape(perplex.np*perplex.nt)
    perplex.alphamelt = alphamelt.reshape(perplex.np*perplex.nt)
    perplex.betamelt  = betamelt.reshape(perplex.np*perplex.nt)
    perplex.cpmelt    = cpmelt.reshape(perplex.np*perplex.nt)
    
    perplex.rhoresidual   = rhoresidual.reshape(perplex.np*perplex.nt)
    perplex.alpharesidual = alpharesidual.reshape(perplex.np*perplex.nt)
    perplex.betaresidual  = betaresidual.reshape(perplex.np*perplex.nt)
    perplex.cpresidual    = cpresidual.reshape(perplex.np*perplex.nt)

    return perplex

def remove_window_no_data_from_perplex(perplex,window=None):
    book = {'rho':          perplex.rho,
            'alpha':        perplex.alpha,
            'beta':         perplex.beta,
            'cp':           perplex.cp,
            'melt':         perplex.melt,
            'rhomelt':      perplex.rhomelt,
            'rhoresidual':  perplex.rhoresidual,
            'alphamelt':    perplex.alphamelt,
            'alpharesidual':perplex.alpharesidual,
            'betamelt':     perplex.betamelt,
            'betaresidual': perplex.betaresidual,
            'cpmelt':       perplex.cpmelt,
            'cpresidual':   perplex.cpresidual}
    
    if ( window is not None):
        P = perplex.P.reshape(int(perplex.np), int(perplex.nt))
        T = perplex.T.reshape(int(perplex.np), int(perplex.nt))
    
    for key in book:
        value = book[key].reshape(int(perplex.np), int(perplex.nt))
        # we remove local data
        if ( window is not None):
            value = np.where(((T<window[0])|(T>window[1])|(P<window[2])|(P>window[3])),value,float('nan'))
        if (key=='rho'):
            perplex.rho   = value.ravel()
        elif (key=='alpha'):
            perplex.alpha = value.ravel()
        elif (key=='beta'):
            perplex.beta  = value.ravel()
        elif (key=='cp'):
            perplex.cp  = value.ravel()
        elif (key=='melt'):
            perplex.melt  = value.ravel()
        elif (key=='rhomelt'):
            perplex.rhomelt  = value.ravel()
        elif (key=='rhoresidual'):
            perplex.rhoresidual  = value.ravel()
        elif (key=='alphamelt'):
            perplex.alphamelt = value.ravel()
        elif (key=='alpharesidual'):
            perplex.alpharesidual = value.ravel()
        elif (key=='betamelt'):
            perplex.betamelt  = value.ravel()
        elif (key=='betaresidual'):
            perplex.betaresidual  = value.ravel()
        elif (key=='cpmelt'):
            perplex.cpmelt  = value.ravel()
        elif (key=='cpresidual'):
            perplex.cpresidual  = value.ravel()
            
    return perplex

def remove_NaN_values_from_perplex(perplex,method='cubic'):
    ndata = perplex.nt*perplex.np
    #dT = (max(perplex.T)-perplex.T[0])/(perplex.nt-1)
    #dT = perplex.dT
    #T  = np.arange(perplex.T[0],max(perplex.T)+dT/1e9,dT)
    T  = np.linspace(min(perplex.T),max(perplex.T),perplex.nt)
    #dP = (max(perplex.P)-perplex.P[0])/(perplex.np-1)
    #dP = perplex.dP
    #P  = np.arange(perplex.P[0],max(perplex.P)+dP/1e9,dP) 
    P  = np.linspace(min(perplex.P),max(perplex.P),perplex.np)
    tt, pp    = np.meshgrid(T, P)

    book = {'rho':  perplex.rho,
            'alpha':perplex.alpha,
            'beta': perplex.beta,
            'cp':   perplex.cp,
            'melt': perplex.melt,
            'rhoresidual': perplex.rhoresidual,
            'rhomelt': perplex.rhomelt,
            'alphamelt':perplex.alphamelt,
            'alpharesidual':perplex.alpharesidual,
            'betamelt': perplex.betamelt,
            'betaresidual': perplex.betaresidual,
            'cpmelt': perplex.cpmelt,
            'cpresidual': perplex.cpresidual}
    for key in book:
        value     = book[key].reshape(int(perplex.np), int(perplex.nt))
        array     = np.ma.masked_invalid(value)
        t1        = tt[~array.mask]
        p1        = pp[~array.mask]
        newarr    = array[~array.mask]
        
        count_nan = np.count_nonzero(np.isnan(value))
        
        #print(key+' count nan = '+str(count_nan))
               
        if (count_nan>0 and count_nan<ndata):
            book[key] = interpolate.griddata(
                                     (t1, p1),
                                     newarr.ravel(),
                                     (tt, pp),
                                     method=method,rescale=True).reshape(ndata)
            
            if (key=='rho'):
                perplex.rho   = book[key]
            elif (key=='alpha'):
                perplex.alpha = book[key]
            elif (key=='beta'):
                perplex.beta  = book[key]
            elif (key=='cp'):
                perplex.cp  = book[key]
            elif (key=='melt'):
                perplex.melt  = book[key]
            elif (key=='rhomelt'):
                perplex.rhomelt  = book[key]
            elif (key=='rhoresidual'):
                perplex.rhoresidual  = book[key]
            elif (key=='alphamelt'):
                perplex.alphamelt = book[key]
            elif (key=='alpharesidual'):
                perplex.alpharesidual = book[key]
            elif (key=='betamelt'):
                perplex.betamelt  = book[key]
            elif (key=='betaresidual'):
                perplex.betaresidual  = book[key]
            elif (key=='cpmelt'):
                perplex.cpmelt  = book[key]
            elif (key=='cpresidual'):
                perplex.cpresidual  = book[key]

    return perplex

def nan_instead_of_outlier(x,value=None):
    if ( value == None ):
        average = np.nanmean(x)
        sigma   = np.nanstd(x)
        #print('mean {} sigma {}'.format(average,sigma))
        #if ((average==0) & (sigma==0)):
        #    print(x)
        #value   = average + sigma
        value   = sigma
        #value   = 0
    for j in range(1,2):
        for i in range(j,len(x)-j):
            if ( ( (x[i]<x[i-j]-value) & (x[i]<x[i+j]-value) ) |
                 ( (x[i]>x[i-j]+value) & (x[i]>x[i+j]+value) ) ):
                #x[i] = float('nan')
                x[i] = 99.0e9
    value = np.where(value<90.0e9,value,float('nan'))
    return x

def remove_lowT_from_perplex(perplex):
    book = {'rho':  perplex.rho,
            'alpha':perplex.alpha,
            'beta': perplex.beta,
            'cp':   perplex.cp,
            'melt': perplex.melt,
            'rhomelt': perplex.rhomelt,
            'rhoresidual': perplex.rhoresidual,
            'alphamelt':perplex.alphamelt,
            'alpharesidual':perplex.alpharesidual,
            'betamelt': perplex.betamelt,
            'betaresidual': perplex.betaresidual,
            'cpmelt': perplex.cpmelt,
            'cpresidual': perplex.cpresidual}
    
    for key in book:
        value = book[key]  #.reshape(int(perplex.np), int(perplex.nt))
        value[perplex.T<220+273] = float('nan')

    return perplex
        

def remove_outliers_from_perplex(perplex,window=None):
    book = {'rho':  perplex.rho,
            'alpha':perplex.alpha,
            'beta': perplex.beta,
            'cp':   perplex.cp,
            'melt': perplex.melt,
            'rhomelt': perplex.rhomelt,
            'rhoresidual': perplex.rhoresidual,
            'alphamelt':perplex.alphamelt,
            'alpharesidual':perplex.alpharesidual,
            'betamelt': perplex.betamelt,
            'betaresidual': perplex.betaresidual,
            'cpmelt': perplex.cpmelt,
            'cpresidual': perplex.cpresidual}
    
    book_melt = {
            'melt': perplex.melt,
            'rhomelt': perplex.rhomelt,
            'alphamelt':perplex.alphamelt,
            'betamelt': perplex.betamelt,
            'cpmelt': perplex.cpmelt,
            'alpharesidual':perplex.alpharesidual,
            'betaresidual': perplex.betaresidual,
            'cpresidual': perplex.cpresidual
             }
    if ( window is not None):
        P = perplex.P.reshape(int(perplex.np), int(perplex.nt))
        T = perplex.T.reshape(int(perplex.np), int(perplex.nt))
    
    for key in book:
        value = book[key].reshape(int(perplex.np), int(perplex.nt))
        melt  = perplex.melt.reshape(int(perplex.np), int(perplex.nt))
        
        # we remove extreme values
        if (key=='rho'):
            maxvalue = 4500e0
            minvalue = 2000e0
            maxvalue_nomelt = 4500e0
            minvalue_nomelt = 2000e0
        elif (key=='rhomelt'):
            maxvalue = 3500e0
            minvalue = 2000e0
            maxvalue_nomelt = 3500e0
            minvalue_nomelt = 2000e0
        elif (key=='rhoresidual'):
            maxvalue = 4500e0
            minvalue = 2000e0
            maxvalue_nomelt = 4500e0
            minvalue_nomelt = 2000e0
        elif (key=='alpha'):
            maxvalue = 7.0e-5
            minvalue = 1.0e-5
            maxvalue_nomelt = 4.6e-5
            minvalue_nomelt = 1.0e-5
        elif (key=='beta'):
            maxvalue = 5.0e-11
            minvalue = 0.3e-11
            maxvalue_nomelt = 1.3e-11
            minvalue_nomelt = 0.3e-11
        elif (key=='cp'):
            maxvalue = 4000
            minvalue = 300
            maxvalue_nomelt = 4000
            minvalue_nomelt = 300
        elif (key=='melt'):
            maxvalue = 1.0e0
            minvalue = 0.0e0
            maxvalue_nomelt = 1.0e0
            minvalue_nomelt = 0.0e0
        elif (key=='alphamelt'):
            maxvalue = 10.0e-5
            minvalue = 1.0e-5
            maxvalue_nomelt = 10.0e-5
            minvalue_nomelt = 1.0e-5
        elif (key=='alpharesidual'):
            maxvalue = 4.8e-5
            minvalue = 1.0e-5
            maxvalue_nomelt = 4.8e-5
            minvalue_nomelt = 1.0e-5
        elif (key=='betamelt'):
            maxvalue = 10.0e-11
            minvalue = 0.3e-11
            maxvalue_nomelt = 10.0e-11
            minvalue_nomelt = 0.3e-11
        elif (key=='betaresidual'):
            maxvalue = 1.3e-11
            minvalue = 0.3e-11
            maxvalue_nomelt = 1.3e-11
            minvalue_nomelt = 0.3e-11
        elif (key=='cpmelt'):
            maxvalue = 4000
            minvalue = 300
            maxvalue_nomelt = 4000
            minvalue_nomelt = 300
        elif (key=='cpresidual'):
            maxvalue = 4000
            minvalue = 300
            maxvalue_nomelt = 4000
            minvalue_nomelt = 300

        #if ( (key=='cp') | (key=='cpresidual') | (key=='cpmelt') ):
        #    print('check values......... key = '+key)
        #    print(value)
        value[np.isnan(value)] = 99.0e9 # the next comparison cannot be done with nan
        value = np.where((value<=maxvalue)&(value>=minvalue),value,float('nan'))
        
        #print(key)
        #print(('nonan {} nan {}'.format((~np.isnan(value)).sum(),(np.isnan(value)).sum())))
        if ( key not in book_melt ):
            value[np.isnan(value)] = 99.0e9 # the next comparison cannot be done with nan
            #value = np.where((melt>0.0e0)|((melt<=0.0e0)&(value<=maxvalue_nomelt)&(value>=minvalue_nomelt)),value,float('nan'))
            value = np.where(((melt>0.0e0)&(value<=maxvalue)&(value>=minvalue))|((melt<=0.0e0)&(value<=maxvalue_nomelt)&(value>=minvalue_nomelt)),value,float('nan'))
        
        # we remove local statistical outlier
        for j in range(0,int(perplex.np)):
            a = [ abs(value[j,i]-value[j,i+1]) for i in range(0,int(perplex.nt)-1) ]
            a.append(0)
            a = nan_instead_of_outlier(a)
            value[j,:] = np.where(~np.isnan(a),value[j,:],float('nan'))
        for i in range(0,int(perplex.nt)):
            a = [ abs(value[j,i]-value[j+1,i]) for j in range(0,int(perplex.np)-1) ]
            a.append(0)
            a = nan_instead_of_outlier(a)
            value[:,i] = np.where(~np.isnan(a),value[:,i],float('nan'))
        
        # we remove local data
        if ( window is not None):
            value = np.where(((T<window[0])|(T>window[1])|(P<window[2])|(P>window[3])),value,float('nan'))
                
        if (key=='rho'):
            perplex.rho   = value.ravel()
        elif (key=='alpha'):
            perplex.alpha = value.ravel()
        elif (key=='beta'):
            perplex.beta  = value.ravel()
        elif (key=='cp'):
            perplex.cp  = value.ravel()
        elif (key=='melt'):
            perplex.melt  = value.ravel()
        elif (key=='rhomelt'):
            perplex.rhomelt  = value.ravel()
        elif (key=='rhoresidual'):
            perplex.rhoresidual  = value.ravel()
        elif (key=='alphamelt'):
            perplex.alphamelt = value.ravel()
        elif (key=='alpharesidual'):
            perplex.alpharesidual = value.ravel()
        elif (key=='betamelt'):
            perplex.betamelt  = value.ravel()
        elif (key=='betaresidual'):
            perplex.betaresidual  = value.ravel()
        elif (key=='cpmelt'):
            perplex.cpmelt  = value.ravel()
        elif (key=='cpresidual'):
            perplex.cpresidual  = value.ravel()
            
    return perplex

def remove_extrapolated_data_from_merge(data,n,m):
    # depending resolutions of merged grids and host grid
    # the new grid can contain extrapolated value

    # from low T
    #checked = np.ones(10, dtype=bool)
    checked = np.full(m,False)
    for i in range(0,n):
        for j in range(0,m):
            if ( np.isnan(data[j,i]) ):
                pass
            else:
                if ( i==0 ):
                    checked[j] = True
                else:
                    if ( checked[j] ):
                        pass
                    else:
                        data[j,i]  = np.float('nan')
                        checked[j] = True
        
    # from high T
    checked = np.full(m,False)
    for i in reversed(range(0,n)):
        for j in range(0,m):
            if ( np.isnan(data[j,i]) ):
                pass
            else:
                if ( i==n-1 ):
                    checked[j] = True
                else:
                    if ( checked[j] ):
                        pass
                    else:
                        data[j,i]  = np.float('nan')
                        checked[j] = True

    # from low P
    checked = np.full(n,False)
    for j in range(0,m):
        for i in range(0,n):
            if ( np.isnan(data[j,i]) ):
                pass
            else:
                if ( j==0 ):
                    checked[i] = True
                else:
                    if ( checked[i] ):
                        pass
                    else:
                        data[j,i]  = np.float('nan')
                        checked[i] = True

    # from high P
    checked = np.full(n,False)
    for j in reversed(range(0,m)):
        for i in range(0,n):
            if ( j==m ):
                checked[i] = True
            else:
                if ( np.isnan(data[j,i]) ):
                    pass
                else:
                    if ( checked[i] ):
                        pass
                    else:
                        data[j,i]  = np.float('nan')
                        checked[i] = True

    return data

def merge_perplex(perplex1,perplex2,tlim=[273,1923], plim=[0,6e9],method='linear',nt_used=None,np_used=None,melt=1):
    # we define dT <= smallest dT (same for dP)
    # and we define a new grid
    if ( ( nt_used is not None ) & ( np_used is not None ) ):
        dT = (tlim[1]-tlim[0])/(nt_used-1)
        dP = (plim[1]-plim[0])/(np_used-1)
    else:
        dT = 9e99
        dP = 9e99
        perplex_list = {perplex1,perplex2}
        for data in perplex_list:
            if (data.dT<dT):
                dT = data.dT
            if (data.dP<dP):
                dP = data.dP
        nt_used = int(np.ceil((tlim[1]-tlim[0])/dT))
        np_used = int(np.ceil((plim[1]-plim[0])/dP))
    
    T              = np.linspace(tlim[0],tlim[1],nt_used)
    P              = np.linspace(plim[0],plim[1],np_used)
    tt, pp         = np.meshgrid(T, P)
    new_perplex    = perplex_output()
    new_perplex.T  = tt.ravel() ; new_perplex.P  = pp.ravel()
    new_perplex.dT = T[1]-T[0]  ; new_perplex.dP = P[1]-P[0]
    new_perplex.nt = nt_used    ; new_perplex.np = np_used
    ndata_new      = new_perplex.nt*new_perplex.np
    
    # we check either there is overlapping or gaps
    # NOT IMPLEMENTED !!!
    
    # count_nan = np.count_nonzero(np.isnan(perplex1.rho))
    # print('perplex1.rho count nan = '+str(count_nan))
    # count_nan = np.count_nonzero(np.isnan(perplex2.rho))
    # print('perplex2.rho count nan = '+str(count_nan))
    
    if (melt==0):
        book1 = {'rho':  perplex1.rho,
                 'alpha':perplex1.alpha,
                 'beta': perplex1.beta,
                 'cp':   perplex1.cp,
                 }
        book2 = {'rho':  perplex2.rho,
                 'alpha':perplex2.alpha,
                 'beta': perplex2.beta,
                 'cp':   perplex2.cp,
                  }
    else:
        book1 = {'rho':  perplex1.rho,
                 'alpha':perplex1.alpha,
                 'beta': perplex1.beta,
                 'cp': perplex1.cp,
                 'melt': perplex1.melt,
                 'rhoresidual': perplex1.rhoresidual,
                 'rhomelt': perplex1.rhomelt,
                 'alphamelt':perplex1.alphamelt,
                 'alpharesidual':perplex1.alpharesidual,
                 'betamelt': perplex1.betamelt,
                 'betaresidual': perplex1.betaresidual,
                 'cpmelt': perplex1.cpmelt,
                 'cpresidual': perplex1.cpresidual}
        book2 = {'rho':  perplex2.rho,
                 'alpha':perplex2.alpha,
                 'beta': perplex2.beta,
                 'cp': perplex2.cp,
                 'melt': perplex2.melt,
                 'rhoresidual': perplex2.rhoresidual,
                 'rhomelt': perplex2.rhomelt,
                 'alphamelt':perplex2.alphamelt,
                 'alpharesidual':perplex2.alpharesidual,
                 'betamelt': perplex2.betamelt,
                 'betaresidual': perplex2.betaresidual,
                 'cpmelt': perplex2.cpmelt,
                 'cpresidual': perplex2.cpresidual}
    
    T1       = np.linspace(min(perplex1.T),max(perplex1.T),perplex1.nt)
    P1       = np.linspace(min(perplex1.P),max(perplex1.P),perplex1.np)
    tt1, pp1 = np.meshgrid(T1, P1)
    
    T2       = np.linspace(min(perplex2.T),max(perplex2.T),perplex2.nt)
    P2       = np.linspace(min(perplex2.P),max(perplex2.P),perplex2.np)
    tt2, pp2 = np.meshgrid(T2, P2)
    
    for key in book1:
        value1    = book1[key].reshape(int(perplex1.np), int(perplex1.nt))
        array1    = np.ma.masked_invalid(value1)
        t1        = tt1[~array1.mask]
        p1        = pp1[~array1.mask]
        newarr1   = array1[~array1.mask]
        
        value2    = book2[key].reshape(int(perplex2.np), int(perplex2.nt))
        array2    = np.ma.masked_invalid(value2)
        t2        = tt2[~array2.mask]
        p2        = pp2[~array2.mask]
        newarr2   = array2[~array2.mask]
        
        if (perplex1 != perplex2):
            T_merge   = np.concatenate([t1,t2])
            P_merge   = np.concatenate([p1,p2])
        else:
            T_merge     = t1
            P_merge     = p1

        if (key=='rho'):
            if (perplex1 != perplex2):
                rho_merge   = np.concatenate([newarr1,newarr2])
            else:
                rho_merge   = newarr1
            new_perplex.rho = interpolate.griddata((T_merge,P_merge),rho_merge,(tt, pp),method=method,rescale=True)
            new_perplex.rho = remove_extrapolated_data_from_merge(new_perplex.rho,new_perplex.nt,new_perplex.np).reshape(ndata_new)
        elif (key=='alpha'):
            if (perplex1 != perplex2):
                alpha_merge   = np.concatenate([newarr1,newarr2])
            else:
                alpha_merge   = newarr1
            new_perplex.alpha = interpolate.griddata((T_merge,P_merge),alpha_merge,(tt, pp),method=method,rescale=True)
            new_perplex.alpha = remove_extrapolated_data_from_merge(new_perplex.alpha,new_perplex.nt,new_perplex.np).reshape(ndata_new)
        elif (key=='beta'):
            if (perplex1 != perplex2):
                beta_merge   = np.concatenate([newarr1,newarr2])
            else:
                beta_merge   = newarr1
            new_perplex.beta = interpolate.griddata((T_merge,P_merge),beta_merge,(tt, pp),method=method,rescale=True)
            new_perplex.beta = remove_extrapolated_data_from_merge(new_perplex.beta,new_perplex.nt,new_perplex.np).reshape(ndata_new)
        elif (key=='cp'):
            if (perplex1 != perplex2):
                cp_merge   = np.concatenate([newarr1,newarr2])
            else:
                cp_merge   = newarr1
            new_perplex.cp = interpolate.griddata((T_merge,P_merge),cp_merge,(tt, pp),method=method,rescale=True)
            new_perplex.cp = remove_extrapolated_data_from_merge(new_perplex.cp,new_perplex.nt,new_perplex.np).reshape(ndata_new)
        elif (key=='melt'):
            if (perplex1 != perplex2):
                melt_merge   = np.concatenate([newarr1,newarr2])
            else:
                melt_merge   = newarr1
            new_perplex.melt = interpolate.griddata((T_merge,P_merge),melt_merge,(tt, pp),method=method,rescale=True)
            new_perplex.melt = remove_extrapolated_data_from_merge(new_perplex.melt,new_perplex.nt,new_perplex.np).reshape(ndata_new)
            new_perplex.melt[np.isnan(new_perplex.melt)] = 0.0e0
            new_perplex.melt[new_perplex.melt<0.0001e0]    = 0.0e0
            #count_nan = np.count_nonzero(np.isnan(new_perplex.melt))
            #print('melt count nan = '+str(count_nan))
        elif (key=='rhomelt'):
            if (perplex1 != perplex2):
                rhomelt_merge   = np.concatenate([newarr1,newarr2])
            else:
                rhomelt_merge   = newarr1
            new_perplex.rhomelt = interpolate.griddata((T_merge,P_merge),rhomelt_merge,(tt, pp),method=method,rescale=True)
            new_perplex.rhomelt = remove_extrapolated_data_from_merge(new_perplex.rhomelt,new_perplex.nt,new_perplex.np).reshape(ndata_new)
            new_perplex.rhomelt[new_perplex.melt<0.0001e0] = float('nan')
        elif (key=='rhoresidual'):
            if (perplex1 != perplex2):
                rhoresidual_merge   = np.concatenate([newarr1,newarr2])
            else:
                rhoresidual_merge   = newarr1
            new_perplex.rhoresidual = interpolate.griddata((T_merge,P_merge),rhoresidual_merge,(tt, pp),method=method,rescale=True)
            new_perplex.rhoresidual = remove_extrapolated_data_from_merge(new_perplex.rhoresidual,new_perplex.nt,new_perplex.np).reshape(ndata_new)
        elif (key=='alphamelt'):
            if (perplex1 != perplex2):
                alphamelt_merge   = np.concatenate([newarr1,newarr2])
            else:
                alphamelt_merge   = newarr1
            new_perplex.alphamelt = interpolate.griddata((T_merge,P_merge),alphamelt_merge,(tt, pp),method=method,rescale=True)
            new_perplex.alphamelt = remove_extrapolated_data_from_merge(new_perplex.alphamelt,new_perplex.nt,new_perplex.np).reshape(ndata_new)
            new_perplex.alphamelt[new_perplex.melt<0.0001e0] = float('nan')
        elif (key=='alpharesidual'):
            if (perplex1 != perplex2):
                alpharesidual_merge   = np.concatenate([newarr1,newarr2])
            else:
                alpharesidual_merge   = newarr1
            new_perplex.alpharesidual = interpolate.griddata((T_merge,P_merge),alpharesidual_merge,(tt, pp),method=method,rescale=True)
            new_perplex.alpharesidual = remove_extrapolated_data_from_merge(new_perplex.alpharesidual,new_perplex.nt,new_perplex.np).reshape(ndata_new)
        elif (key=='betamelt'):
            if (perplex1 != perplex2):
                betamelt_merge   = np.concatenate([newarr1,newarr2])
            else:
                betamelt_merge   = newarr1
            new_perplex.betamelt = interpolate.griddata((T_merge,P_merge),betamelt_merge,(tt, pp),method=method,rescale=True)
            new_perplex.betamelt = remove_extrapolated_data_from_merge(new_perplex.betamelt,new_perplex.nt,new_perplex.np).reshape(ndata_new)
            new_perplex.betamelt[new_perplex.melt<0.0001e0] = float('nan')
        elif (key=='betaresidual'):
            if (perplex1 != perplex2):
                betaresidual_merge   = np.concatenate([newarr1,newarr2])
            else:
                betaresidual_merge   = newarr1
            new_perplex.betaresidual = interpolate.griddata((T_merge,P_merge),betaresidual_merge,(tt, pp),method=method,rescale=True)
            new_perplex.betaresidual = remove_extrapolated_data_from_merge(new_perplex.betaresidual,new_perplex.nt,new_perplex.np).reshape(ndata_new)
            #new_perplex.betaresidual = np.where((new_perplex.melt<0.01e0),new_perplex.beta,new_perplex.betaresidual)
        elif (key=='cpmelt'):
            if (perplex1 != perplex2):
                cpmelt_merge   = np.concatenate([newarr1,newarr2])
            else:
                cpmelt_merge   = newarr1
            #if (cpmelt_merge==[]):
            if (cpmelt_merge.shape[0]==0):
                print('WARNING: check what is wrong cpmelt aray is fully NaN')
                new_perplex.cpmelt = np.zeros_like(new_perplex.melt)
                new_perplex.cpmelt[new_perplex.melt<0.0001e0] = float('nan')
            else:
                new_perplex.cpmelt = interpolate.griddata((T_merge,P_merge),cpmelt_merge,(tt, pp),method=method,rescale=True)
                new_perplex.cpmelt = remove_extrapolated_data_from_merge(new_perplex.cpmelt,new_perplex.nt,new_perplex.np).reshape(ndata_new)
                new_perplex.cpmelt[new_perplex.melt<0.0001e0] = float('nan')
        elif (key=='cpresidual'):
            if (perplex1 != perplex2):
                cpresidual_merge   = np.concatenate([newarr1,newarr2])
            else:
                cpresidual_merge   = newarr1
            new_perplex.cpresidual = interpolate.griddata((T_merge,P_merge),cpresidual_merge,(tt, pp),method=method,rescale=True)
            new_perplex.cpresidual = remove_extrapolated_data_from_merge(new_perplex.cpresidual,new_perplex.nt,new_perplex.np).reshape(ndata_new)
            #new_perplex.cpresidual = np.where((new_perplex.melt<0.01e0),new_perplex.cp,new_perplex.cpresidual)

    
    new_perplex.dT    = dT
    new_perplex.dP    = dP
    new_perplex.nt    = nt_used
    new_perplex.np    = np_used
    
    if (melt==0):
        # new_perplex.alphamelt.fill(0.0e0)
        # new_perplex.betamelt.fill(0.0e0)
        # new_perplex.melt.fill(0.0e0)
        # new_perplex.rhomelt.fill(0.0e0)
        new_perplex.alphamelt = np.zeros_like(new_perplex.rho)
        new_perplex.betamelt  = np.zeros_like(new_perplex.rho)
        new_perplex.cpmelt    = np.zeros_like(new_perplex.rho)
        new_perplex.melt      = np.zeros_like(new_perplex.rho)
        new_perplex.rhomelt   = np.zeros_like(new_perplex.rho)
        
        new_perplex.rhoresidual   = deepcopy(new_perplex.rho)
        new_perplex.alpharesidual = deepcopy(new_perplex.alpha)
        new_perplex.betaresidual  = deepcopy(new_perplex.beta)
        new_perplex.cpresidual    = deepcopy(new_perplex.cp)
        

    return new_perplex

def update_sampling_perplex(perplex,dT,dP,method='linear'):
    """ Update sampling. useful to compare perplex grid. """
    tlim=[np.min(perplex.T),np.max(perplex.T)]
    plim=[np.min(perplex.P),np.max(perplex.P)]

    if ((tlim[1]-tlim[0])%dT != 0):
        dT = (tlim[1]-tlim[0])/np.ceil((tlim[1]-tlim[0])/dT)
    if ((plim[1]-plim[0])%dP != 0):
        dP = (plim[1]-plim[0])/np.ceil((plim[1]-plim[0])/dP)
    
    T              = np.arange(tlim[0],tlim[1]+dT,dT)
    P              = np.arange(plim[0],plim[1]+dP,dP)
    tt, pp         = np.meshgrid(T, P)
    new_perplex    = perplex_output()
    new_perplex.T  = tt.ravel() ; new_perplex.P  = pp.ravel()
    new_perplex.dT = dT         ; new_perplex.dP = dP
    new_perplex.nt = int(len(T))     ; new_perplex.np = int(len(P))
    ndata_new      = new_perplex.nt*new_perplex.np
    
    T_old     = perplex.T
    P_old     = perplex.P
    rho_old   = perplex.rho
    alpha_old = perplex.alpha
    beta_old  = perplex.beta
    cp_old    = perplex.cp
    melt_old  = perplex.melt
    rhomelt_old     = perplex.rhomelt
    rhoresidual_old = perplex.rhoresidual
    alphamelt_old = perplex.alphamelt
    alpharesidual_old = perplex.alpharesidual
    betamelt_old      = perplex.betamelt
    betaresidual_old  = perplex.betaresidual
    cpmelt_old      = perplex.cpmelt
    cpresidual_old  = perplex.cpresidual
    
    new_perplex.rho   = interpolate.griddata((T_old,P_old),rho_old,(tt, pp),method=method,rescale=True).reshape(ndata_new)
    new_perplex.alpha = interpolate.griddata((T_old,P_old),alpha_old,(tt, pp),method=method,rescale=True).reshape(ndata_new)
    new_perplex.beta  = interpolate.griddata((T_old,P_old),beta_old,(tt, pp),method=method,rescale=True).reshape(ndata_new)
    new_perplex.cp    = interpolate.griddata((T_old,P_old),cp_old,(tt, pp),method=method,rescale=True).reshape(ndata_new)
    new_perplex.melt  = interpolate.griddata((T_old,P_old),melt_old,(tt, pp),method=method,rescale=True).reshape(ndata_new)
    new_perplex.rhomelt  = interpolate.griddata((T_old,P_old),rhomelt_old,(tt, pp),method=method,rescale=True).reshape(ndata_new)
    new_perplex.rhoresidual  = interpolate.griddata((T_old,P_old),rhoresidual_old,(tt, pp),method=method,rescale=True).reshape(ndata_new)
    new_perplex.alphamelt = interpolate.griddata((T_old,P_old),alphamelt_old,(tt, pp),method=method,rescale=True).reshape(ndata_new)
    new_perplex.alpharesidual = interpolate.griddata((T_old,P_old),alpharesidual_old,(tt, pp),method=method,rescale=True).reshape(ndata_new)
    new_perplex.betamelt      = interpolate.griddata((T_old,P_old),betamelt_old,(tt, pp),method=method,rescale=True).reshape(ndata_new)
    new_perplex.betaresidual  = interpolate.griddata((T_old,P_old),betaresidual_old,(tt, pp),method=method,rescale=True).reshape(ndata_new)
    new_perplex.cpmelt      = interpolate.griddata((T_old,P_old),cpmelt_old,(tt, pp),method=method,rescale=True).reshape(ndata_new)
    new_perplex.cpresidual  = interpolate.griddata((T_old,P_old),cpresidual_old,(tt, pp),method=method,rescale=True).reshape(ndata_new)
    
    return new_perplex

def diff_perplex(perplex1,perplex2):
    
    if ( (perplex1.nt == perplex2.nt) &
         (perplex1.np == perplex2.np) &
         (perplex1.dT == perplex2.dT) &
         (perplex1.dP == perplex2.dP)):
        new_perplex        = perplex_output()
        new_perplex.T      = deepcopy(perplex1.T)
        new_perplex.P      = deepcopy(perplex1.P)
        new_perplex.rho    = np.substract(perplex1.rho,perplex2.rho)
        new_perplex.alpha  = np.substract(perplex1.alpha,perplex2.alpha)
        new_perplex.beta   = np.substract(perplex1.beta,perplex2.beta)
        new_perplex.cp     = np.substract(perplex1.cp,perplex2.cp)
        new_perplex.rhomelt       = np.substract(perplex1.rhomelt,perplex2.rhomelt)
        new_perplex.rhoresidual   = np.substract(perplex1.rhoresidual,perplex2.rhoresidual)
        new_perplex.alphamelt  = np.substract(perplex1.alphamelt,perplex2.alphamelt)
        new_perplex.alpharesidual  = np.substract(perplex1.alpharesidual,perplex2.alpharesidual)
        new_perplex.betamelt   = np.substract(perplex1.betamelt,perplex2.betamelt)
        new_perplex.betaresidual   = np.substract(perplex1.betaresidual,perplex2.betaresidual)
        new_perplex.cpmelt       = np.substract(perplex1.cpmelt,perplex2.cpmelt)
        new_perplex.cpresidual   = np.substract(perplex1.cpresidual,perplex2.cpresidual)
        
    else:
        print('The two input perplex objects must have the same size.')
        return
    
    return new_perplex

def smooth_perplex(perplex):
    
    #sigma = 3
    #rho   = np.zeros_like(perplex.rho)
    #gaussian_filter(perplex.rho,sigma=sigma,order=0,output=rho,mode='reflect',cval=0.0,truncate=1.0)
    
    rho         = medfilt2d(perplex.rho.reshape(int(perplex.np), int(perplex.nt)),kernel_size=(9,9))
    perplex.rho = rho.reshape(perplex.np*perplex.nt)
    
    rhoresidual         = medfilt2d(perplex.rhoresidual.reshape(int(perplex.np), int(perplex.nt)),kernel_size=(9,9))
    perplex.rhoresidual = rhoresidual.reshape(perplex.np*perplex.nt)
    
    return perplex

def extrapolate_lowT_perplex2(perplex):
    """ Extrapolate low temperature values that perplex cannot
    compute. The assumption is that alpha and beta are constant and
    follow the general trend at low temperature"""

    extrapolated_perplex = perplex_output()
    extrapolated_perplex = deepcopy(perplex)
    Tref = 273
    T0   = 0
    Pref = 0
    trend = (5.802e9 - 1.646e9)/(1181e0-237e0)

    rho   = extrapolated_perplex.rho.reshape(  int(perplex.np), int(perplex.nt))
    T     = extrapolated_perplex.T.reshape(    int(perplex.np), int(perplex.nt))
    tmax  = np.max(T) ;
    P     = extrapolated_perplex.P.reshape(    int(perplex.np), int(perplex.nt))
    pmax  = np.max(P) ;
    alpha = extrapolated_perplex.alpha.reshape(int(perplex.np), int(perplex.nt))
    beta  = extrapolated_perplex.beta.reshape( int(perplex.np), int(perplex.nt))
    cp    = extrapolated_perplex.cp.reshape( int(perplex.np), int(perplex.nt))
    melt  = extrapolated_perplex.melt.reshape( int(perplex.np), int(perplex.nt))
    rhomelt       = extrapolated_perplex.rhomelt.reshape( int(perplex.np), int(perplex.nt))
    rhoresidual   = extrapolated_perplex.rhoresidual.reshape( int(perplex.np), int(perplex.nt))
    alphamelt     = extrapolated_perplex.alphamelt.reshape(int(perplex.np), int(perplex.nt))
    alpharesidual = extrapolated_perplex.alpharesidual.reshape(int(perplex.np), int(perplex.nt))
    betamelt      = extrapolated_perplex.betamelt.reshape( int(perplex.np), int(perplex.nt))
    betaresidual  = extrapolated_perplex.betaresidual.reshape( int(perplex.np), int(perplex.nt))
    cpmelt        = extrapolated_perplex.cpmelt.reshape( int(perplex.np), int(perplex.nt))
    cpresidual    = extrapolated_perplex.cpresidual.reshape( int(perplex.np), int(perplex.nt))
    

    # values used for T = 273 K
    extra_P_at_T0       = []
    extra_rho           = []
    extra_alpha         = []
    extra_beta          = []
    extra_cp            = []
    extra_melt          = []
    extra_rhomelt       = []
    extra_rhoresidual   = []
    extra_alphamelt     = []
    extra_alpharesidual = []
    extra_betamelt      = []
    extra_betaresidual  = []
    extra_cpmelt        = []
    extra_cpresidual    = []
    for j in range(0,int(perplex.np)):
        for i in range(0,int(perplex.nt)):
            if ( np.isnan(rho[j,i]) ):
                pass
            else:
                extra_P_at_T0.append((P[j,i]-Pref)-trend*(T[j,i]-Tref))
                extra_rho.append(rho[j,i])
                extra_alpha.append(alpha[j,i])
                extra_beta.append(beta[j,i])
                extra_cp.append(cp[j,i])
                extra_melt.append(melt[j,i])
                extra_rhomelt.append(rhomelt[j,i])
                extra_rhoresidual.append(rhoresidual[j,i])
                extra_alphamelt.append(alphamelt[j,i])
                extra_alpharesidual.append(alpharesidual[j,i])
                extra_betamelt.append(betamelt[j,i])
                extra_betaresidual.append(betaresidual[j,i])
                extra_cpmelt.append(cpmelt[j,i])
                extra_cpresidual.append(cpresidual[j,i])
                break
   
    max_extra_P_at_T0 = np.max(extra_P_at_T0)
    # we interpolate 1-D values at T0 on merged grid
    for j in range(0,int(perplex.np)):
        if ( P[j,0] <= max_extra_P_at_T0 ):
            rho[j,0]           = pos2value_1d(extra_P_at_T0,extra_rho,P[j,0])
            alpha[j,0]         = pos2value_1d(extra_P_at_T0,extra_alpha,P[j,0])
            beta[j,0]          = pos2value_1d(extra_P_at_T0,extra_beta,P[j,0])
            cp[j,0]            = pos2value_1d(extra_P_at_T0,extra_cp,P[j,0])
            melt[j,0]          = pos2value_1d(extra_P_at_T0,extra_melt,P[j,0])
            rhomelt[j,0]       = pos2value_1d(extra_P_at_T0,extra_rhomelt,P[j,0])
            rhoresidual[j,0]   = pos2value_1d(extra_P_at_T0,extra_rhoresidual,P[j,0])
            alphamelt[j,0]     = pos2value_1d(extra_P_at_T0,extra_alphamelt,P[j,0])
            alpharesidual[j,0] = pos2value_1d(extra_P_at_T0,extra_alpharesidual,P[j,0])
            betamelt[j,0]      = pos2value_1d(extra_P_at_T0,extra_betamelt,P[j,0])
            betaresidual[j,0]  = pos2value_1d(extra_P_at_T0,extra_betaresidual,P[j,0])
            cpmelt[j,0]        = pos2value_1d(extra_P_at_T0,extra_cpmelt,P[j,0])
            cpresidual[j,0]    = pos2value_1d(extra_P_at_T0,extra_cpresidual,P[j,0])
        else:
            break

    # we interpolate at low temperature between T=273 and the first low T value
    extrapolated_perplex = remove_NaN_values_from_perplex(extrapolated_perplex)


    return extrapolated_perplex

def extrapolate_lowT_perplex(perplex):
    """ Extrapolate low temperature values that perplex cannot 
    compute. The assumption is that alpha and beta are constant and temperature-independent 
    at low temperature which is certainly (reasonably) wrong."""
    
    extrapolated_perplex = perplex_output()
    extrapolated_perplex = deepcopy(perplex)
    
    Tref = 273
    Pref = 0
    
    rho   = extrapolated_perplex.rho.reshape(  int(perplex.np), int(perplex.nt))
    T     = extrapolated_perplex.T.reshape(    int(perplex.np), int(perplex.nt))
    tmax  = np.max(T) ;
    P     = extrapolated_perplex.P.reshape(    int(perplex.np), int(perplex.nt))
    pmax  = np.max(P) ;
    alpha = extrapolated_perplex.alpha.reshape(int(perplex.np), int(perplex.nt))
    beta  = extrapolated_perplex.beta.reshape( int(perplex.np), int(perplex.nt))
    cp    = extrapolated_perplex.cp.reshape( int(perplex.np), int(perplex.nt))
    melt  = extrapolated_perplex.melt.reshape( int(perplex.np), int(perplex.nt))
    rhomelt       = extrapolated_perplex.rhomelt.reshape( int(perplex.np), int(perplex.nt))
    rhoresidual   = extrapolated_perplex.rhoresidual.reshape( int(perplex.np), int(perplex.nt))
    alphamelt     = extrapolated_perplex.alphamelt.reshape(int(perplex.np), int(perplex.nt))
    alpharesidual = extrapolated_perplex.alpharesidual.reshape(int(perplex.np), int(perplex.nt))
    betamelt      = extrapolated_perplex.betamelt.reshape( int(perplex.np), int(perplex.nt))
    betaresidual  = extrapolated_perplex.betaresidual.reshape( int(perplex.np), int(perplex.nt))
    cpmelt        = extrapolated_perplex.cpmelt.reshape( int(perplex.np), int(perplex.nt))
    cpresidual    = extrapolated_perplex.cpresidual.reshape( int(perplex.np), int(perplex.nt))
    
    # smoothing alpha and beta along the boundaries to avoid vertical discontinuities not suitable
    n_smooth = 2
    rho_smooth           = []
    rhomelt_smooth       = []
    rhoresidual_smooth   = []
    alpha_smooth         = []
    beta_smooth          = []
    cp_smooth            = []
    alphamelt_smooth     = []
    betamelt_smooth      = []
    alpharesidual_smooth = []
    betaresidual_smooth  = []
    cpmelt_smooth        = []
    cpresidual_smooth    = []
    
    # -----------------------------------------------------------------------
    # find alpha and beta variations on the lowT boundary
    alpha_lowT = []
    alpha_lowT_prev = []
    beta_lowT  = []
    beta_lowT_prev  = []
    iref = []
    for j in range(0,int(perplex.np)):
        for i in range(0,int(perplex.nt)):
            if ( np.isnan(rho[j,i]) ):
                #print('None T {} P {}'.format(T[j,i],P[j,i]))
                pass
            else:
                alpha_lowT.append(alpha[j,i])
                alpha_lowT_prev.append(alpha[j,i+1])            
                  
                beta_lowT.append(beta[j,i])
                beta_lowT_prev.append(beta[j,i+1])
                  
                iref.append(i)
                jmax = j
                break
    # -----------------------------------------------------------------------
    
    
    i_smooth = 0
    i_int    = 0
    for j in range(0,int(perplex.np)):
        are_values = False
        for i in range(0,int(perplex.nt)):
              if ( np.isnan(rho[j,i]) ):
                  #print('None T {} P {}'.format(T[j,i],P[j,i]))
                  pass
              else:
                  if (i_smooth<n_smooth):
                      alpha_smooth.append(alpha[j,i])
                      beta_smooth.append(beta[j,i])
                      cp_smooth.append(cp[j,i])
                      alphamelt_smooth.append(alphamelt[j,i])
                      betamelt_smooth.append(betamelt[j,i])
                      alpharesidual_smooth.append(alpharesidual[j,i])
                      betaresidual_smooth.append(betaresidual[j,i])
                      cpmelt_smooth.append(cpmelt[j,i])
                      cpresidual_smooth.append(cpresidual[j,i])
                      rho_smooth.append(rho[j,i])
                      rhomelt_smooth.append(rhomelt[j,i])
                      rhoresidual_smooth.append(rhoresidual[j,i])
                      i_smooth = i_smooth + 1
                  else:
                      alpha_smooth[i_int]         = alpha[j,i]
                      beta_smooth[i_int]          = beta[j,i]
                      cp_smooth[i_int]            = cp[j,i]
                      alphamelt_smooth[i_int]     = alphamelt[j,i]
                      betamelt_smooth[i_int]      = betamelt[j,i]
                      alpharesidual_smooth[i_int] = alpharesidual[j,i]
                      betaresidual_smooth[i_int]  = betaresidual[j,i]
                      cpmelt_smooth[i_int]        = cpmelt[j,i]
                      cpresidual_smooth[i_int]    = cpresidual[j,i]
                      rho_smooth[i_int]           = rho[j,i]
                      rhomelt_smooth[i_int]       = rhomelt[j,i]
                      rhoresidual_smooth[i_int]   = rhoresidual[j,i]
                      i_int = i_int + 1
                      if (i_int>=n_smooth):
                          i_int = 0
                    
                  alpha_used = sum(alpha_smooth)/len(alpha_smooth)
                  beta_used  = sum(beta_smooth)/len(beta_smooth)
                  cp_used    = sum(cp_smooth)/len(cp_smooth)
                  rho_ref    = sum(rho_smooth)/len(rho_smooth) / ( (1+beta_used*(P[j,i]-Pref)) * (1-alpha_used*(T[j,i]-Tref)) )
                  
                  alpha_used_melt = sum(alphamelt_smooth)/len(alphamelt_smooth)
                  beta_used_melt  = sum(betamelt_smooth)/len(betamelt_smooth)
                  cp_used_melt    = sum(cpmelt_smooth)/len(cpmelt_smooth)
                  rho_ref_melt    = sum(rhomelt_smooth)/len(rhomelt_smooth) / ( (1+beta_used_melt*(P[j,i]-Pref)) * (1-alpha_used_melt*(T[j,i]-Tref)) )
                 
                  alpha_used_residual = sum(alpharesidual_smooth)/len(alpharesidual_smooth)
                  beta_used_residual  = sum(betaresidual_smooth)/len(betaresidual_smooth)
                  cp_used_residual    = sum(cpresidual_smooth)/len(cpresidual_smooth)
                  rho_ref_residual    = sum(rhoresidual_smooth)/len(rhoresidual_smooth) / ( (1+beta_used_residual*(P[j,i]-Pref)) * (1-alpha_used_residual*(T[j,i]-Tref)) )

                  are_values  = True
                  break
        if (are_values):
            for i in range(0,int(perplex.nt)):
                  if ( np.isnan(rho[j,i]) ):
                      melt[j,i]  = 0.0e0
                     
                      #alpha[j,i] = alpha_used
                      #beta[j,i]  = beta_used
                      for k in range(1,1000):
                          if ((j+k>jmax)|(i+k>=perplex.nt)):
                              alpha[j,i] = alpha_lowT[jmax]
                              beta[j,i]  = beta_lowT[jmax]
                              break
                          else:
                              if ( np.isnan(rho[j+k,i+k]) ):
                                  pass
                              else:
                                  alpha[j,i] = alpha[j+k,i+k]
                                  beta[j,i]  = beta[j+k,i+k]
                                  break

                      rho[j,i]   = rho_ref *(1+beta[j,i]*(P[j,i]-Pref))*(1-alpha[j,i]*(T[j,i]-Tref))
                     
                      cp[j,i]    = cp_used
                    
                      rhomelt[j,i]   = rho_ref_melt*(1+beta_used_melt*(P[j,i]-Pref))*(1-alpha_used_melt*(T[j,i]-Tref))
                      alphamelt[j,i] = alpha_used_melt
                      betamelt[j,i]  = beta_used_melt
                      cpmelt[j,i]    = cp_used_melt
                     
                      rhoresidual[j,i]   = rho_ref_residual*(1+beta_used_residual*(P[j,i]-Pref))*(1-alpha_used_residual*(T[j,i]-Tref))
                      alpharesidual[j,i] = alpha_used_residual
                      betaresidual[j,i]  = beta_used_residual
                      cpresidual[j,i]    = cp_used_residual
                  else:
                      break

    return extrapolated_perplex

def compute_rho_ref_average_perplex(perplex,tlim=[573,1573], plim=[0.5e9,1e9],domainBased=True,lithos=None,inside_lithosphere=True,printOut=True):
    """ Compute a reference density based on average at depth, i.e. without
        correcting from alpha and beta,  with its standard deviation.
        1- domain based estimation
        We assume that shallow rocks can be exhumed without any changes.
        So we use those default P,T conditions to compute 
        an average reference density at surface.
        2- T,P along a given geotherm based"""
        
    Tref  = 273
    Pref  = 0
    rho   = deepcopy(perplex.rho)
    T     = deepcopy(perplex.T)
    P     = deepcopy(perplex.P)
    #alpha = deepcopy(perplex.alpha)
    #beta  = deepcopy(perplex.beta)
    melt  = deepcopy(perplex.melt)
  
    # no correction 
    rho_ref = rho
 
    if (domainBased):
        # domain selection
        sort_rho_ref  = rho_ref[ (melt==0) &
                                 (P<plim[1]) & (P>plim[0]) &
                                 (T<tlim[1]) & (T>tlim[0]) &
                                 (~np.isnan(rho_ref)) ]
    else:
        if ( not lithos ):
            print('Lithosphere lith125 used as default to compute rho_ref')
            lithos = geodyn1d.lithosphere('lith125')
        else:
            if ( isinstance(lithos, str) ):
                lithos = geodyn1d.lithosphere(lithos)
        # selection based on geotherm and lithospheric structure
        # if ( not lithos ):
        #     print('Lithosphere based rho_ref requires a lithopshere class')
        #     return 0
        # else:
            # we define the interpolation function
            pmax = np.max(perplex.P)
            tmax = np.max(perplex.T)
            T = np.linspace(np.min(perplex.T),
                            tmax,
                            perplex.nt)
            P = np.linspace(np.min(perplex.P),
                            pmax,
                            perplex.np)
            f = RegularGridInterpolator((np.array(P)/pmax, np.array(T)/tmax), rho_ref.reshape(perplex.np,perplex.nt))
            
            # we check for geotherm and pressure profile
            if ( not lithos.geotherm):
                lithos.get_steady_state_geotherm(printOut=printOut)
            if ( not lithos.P_rho_DepthProfile):
                lithos.get_pressure_density_profile(printOut=printOut)
            
            sort_rho_ref = []
            Tmoho = pos2value_1d(lithos.P_rho_DepthProfile.Z,lithos.geotherm.T,sum(lithos.thicknesses[0:3]))
            Tlab  = pos2value_1d(lithos.P_rho_DepthProfile.Z,lithos.geotherm.T,sum(lithos.thicknesses[0:5]))
            if (Tlab-Tmoho<=1.0e-4):
                Tlab  = 1569.15
            print('tmoho {} tlab {}'.format(Tmoho,Tlab))
            for i in range(0,len(lithos.geotherm.Z)):
                pressure = lithos.P_rho_DepthProfile.P[i]
                temp     = lithos.geotherm.T[i]
                if ( inside_lithosphere ):
                    if ( (temp>=Tmoho) & (temp<=Tlab) ):
                        sort_rho_ref.append(f([pressure/pmax,temp/tmax]))
                else:
                    if (temp>=Tlab):
                        sort_rho_ref.append(f([pressure/pmax,temp/tmax]))
        
    
    rho_ref_value = np.around(np.average(sort_rho_ref),1)
    rho_ref_std   = np.around(np.std(sort_rho_ref),1)
    
    return rho_ref_value, rho_ref_std, rho_ref

def compute_cp_D_ref_perplex(perplex,tlim=[573,1573], plim=[0.5e9,1e9],domainBased=True,lithos=None,printOut=True):
    """ Compute reference heat capacity and diffusivity at surface with its standard deviation.
        1- domain based estimation
        We assume that shallow rocks can be exhumed without any changes.
        So we use those default P,T conditions to compute 
        an average reference density at surface.
        2- T,P along a given geotherm based"""
        
    #Tref  = 273
    #Pref  = 0
    rho   = deepcopy(perplex.rho)
    T     = deepcopy(perplex.T)
    P     = deepcopy(perplex.P)
    #alpha = deepcopy(perplex.alpha)
    #beta  = deepcopy(perplex.beta)
    melt  = deepcopy(perplex.melt)
    
    #factor1 = 1 + np.multiply(beta,P-Pref)
    #factor2 = 1 - np.multiply(alpha,T-Tref)
    #factor  = np.multiply(factor1,factor2)
    #rho_ref = np.divide(rho,factor)

    cp_ref  = deepcopy(perplex.cp)
    D_ref   = np.multiply(rho,cp_ref)
    D_ref   = 2.25e6/D_ref

    if (domainBased):
        # domain selection
        sort_cp_ref  = cp_ref[ (melt==0) &
                                 (P<plim[1]) & (P>plim[0]) &
                                 (T<tlim[1]) & (T>tlim[0]) &
                                 (~np.isnan(cp_ref)) ]
        sort_D_ref  = D_ref[ (melt==0) &
                                 (P<plim[1]) & (P>plim[0]) &
                                 (T<tlim[1]) & (T>tlim[0]) &
                                 (~np.isnan(D_ref)) ]
    else:
        if ( not lithos ):
            print('Lithosphere lith125 used as default to compute rho_ref')
            lithos = geodyn1d.lithosphere('lith125')
        else:
            if ( isinstance(lithos, str) ):
                lithos = geodyn1d.lithosphere(lithos)
        # selection based on geotherm and lithospheric structure
        # if ( not lithos ):
        #     print('Lithosphere based rho_ref requires a lithopshere class')
        #     return 0
        # else:
            # we define the interpolation function
            pmax = np.max(perplex.P)
            tmax = np.max(perplex.T)
            T = np.linspace(np.min(perplex.T),
                            tmax,
                            perplex.nt)
            P = np.linspace(np.min(perplex.P),
                            pmax,
                            perplex.np)
            f   = RegularGridInterpolator((np.array(P)/pmax, np.array(T)/tmax), cp_ref.reshape(perplex.np,perplex.nt))
            f_D = RegularGridInterpolator((np.array(P)/pmax, np.array(T)/tmax), D_ref.reshape(perplex.np,perplex.nt))
            
            # we check for geotherm and pressure profile
            if ( not lithos.geotherm):
                lithos.get_steady_state_geotherm(printOut=printOut)
            if ( not lithos.P_rho_DepthProfile):
                lithos.get_pressure_density_profile(printOut=printOut)
            
            sort_cp_ref = []
            sort_D_ref  = []
            Tmoho = pos2value_1d(lithos.P_rho_DepthProfile.Z,lithos.geotherm.T,sum(lithos.thicknesses[0:3]))
            Tlab  = pos2value_1d(lithos.P_rho_DepthProfile.Z,lithos.geotherm.T,sum(lithos.thicknesses[0:5]))
            if (Tlab-Tmoho<=1.0e-4):
                Tlab  = 1569.15
            print('tmoho {} tlab {}'.format(Tmoho,Tlab))
            for i in range(0,len(lithos.geotherm.Z)):
                pressure = lithos.P_rho_DepthProfile.P[i]
                temp     = lithos.geotherm.T[i]
                if ( (temp>=Tmoho) & (temp<=Tlab) ):
                    sort_cp_ref.append(f([pressure/pmax,temp/tmax])[0])
                    sort_D_ref.append(f_D([pressure/pmax,temp/tmax])[0])
        
    
    cp_ref_value = np.around(np.average(sort_cp_ref),1)
    cp_ref_std   = np.around(np.std(sort_cp_ref),1)

    D_ref_value = np.around(np.average(sort_D_ref),1)
    D_ref_std   = np.around(np.std(sort_D_ref),2)
  
    return cp_ref_value, cp_ref_std, D_ref_value, D_ref_std


def compute_rho_ref_perplex(perplex,tlim=[573,1573], plim=[0.5e9,1e9],domainBased=True,lithos=None,inside_lithosphere=True,printOut=True):
    """ Compute reference density at surface with its standard deviation.
        1- domain based estimation
        We assume that shallow rocks can be exhumed without any changes.
        So we use those default P,T conditions to compute 
        an average reference density at surface.
        2- T,P along a given geotherm based"""
        
    Tref  = 273.15
    Pref  = 0
    rho   = deepcopy(perplex.rho)
    T     = deepcopy(perplex.T)
    P     = deepcopy(perplex.P)
    alpha = deepcopy(perplex.alpha)
    beta  = deepcopy(perplex.beta)
    melt  = deepcopy(perplex.melt)
    
    factor1 = 1 + np.multiply(beta,P-Pref)
    factor2 = 1 - np.multiply(alpha,T-Tref)
    factor  = np.multiply(factor1,factor2)
    rho_ref = np.divide(rho,factor)
    
    if (domainBased):
        # domain selection
        sort_rho_ref  = rho_ref[ (melt==0) &
                                 (P<plim[1]) & (P>plim[0]) &
                                 (T<tlim[1]) & (T>tlim[0]) &
                                 (~np.isnan(rho_ref)) ]
    else:
        if ( not lithos ):
            print('Lithosphere lith125 used as default to compute rho_ref')
            lithos = geodyn1d.lithosphere('lith125')
        else:
            if ( isinstance(lithos, str) ):
                lithos = geodyn1d.lithosphere(lithos)
        # selection based on geotherm and lithospheric structure
        # if ( not lithos ):
        #     print('Lithosphere based rho_ref requires a lithopshere class')
        #     return 0
        # else:
            # we define the interpolation function
            pmax = np.max(perplex.P)
            tmax = np.max(perplex.T)
            T = np.linspace(np.min(perplex.T),
                            tmax,
                            perplex.nt)
            P = np.linspace(np.min(perplex.P),
                            pmax,
                            perplex.np)
            f = RegularGridInterpolator((np.array(P)/pmax, np.array(T)/tmax), rho_ref.reshape(perplex.np,perplex.nt))
            
            # we check for geotherm and pressure profile
            if ( not lithos.geotherm):
                lithos.get_steady_state_geotherm(printOut=printOut)
            if ( not lithos.P_rho_DepthProfile):
                lithos.get_pressure_density_profile(printOut=printOut)
            
            sort_rho_ref = []
            sort_rho_ref_all = []
            Tmoho = pos2value_1d(lithos.P_rho_DepthProfile.Z,lithos.geotherm.T,sum(lithos.thicknesses[0:3]))
            Tlab  = pos2value_1d(lithos.P_rho_DepthProfile.Z,lithos.geotherm.T,sum(lithos.thicknesses[0:5]))
            if (Tlab-Tmoho<=1.0e-4):
                Tlab  = 1569.15
            print('Temperatures moho {} lab {} (degree Celsius)'.format(Tmoho-Tref,Tlab-Tref))
            for i in range(0,len(lithos.geotherm.Z)):
                pressure = lithos.P_rho_DepthProfile.P[i]
                temp     = lithos.geotherm.T[i]
                if ( inside_lithosphere ):
                    if ( (temp>=Tmoho) & (temp<=Tlab) ):
                        sort_rho_ref.append(f([pressure/pmax,temp/tmax]))
                else:
                    if (temp>=Tlab):
                        sort_rho_ref.append(f([pressure/pmax,temp/tmax]))
                if (temp>=Tmoho):
                    sort_rho_ref_all.append(f([pressure/pmax,temp/tmax]))
        
    
    rho_ref_value = np.around(np.average(sort_rho_ref),1)
    rho_ref_std   = np.around(np.std(sort_rho_ref),1)

    rho_ref_value_all = np.around(np.average(sort_rho_ref_all),1)
    rho_ref_std_all   = np.around(np.std(sort_rho_ref_all),1)
    print(' Average ref. density {} +/- {} kg/m3 (lithospheric + sub-lithospheric mantles together)'.format(rho_ref_value_all,rho_ref_std_all))

    return rho_ref_value, rho_ref_std, rho_ref

def compute_delta_rho_perplex(perplex,domainBased=False,lithos=None):
    
    if ( domainBased == False ):
        if ( not lithos ):
            lithos = geodyn1d.lithosphere('lith125')
        else:
            if ( isinstance(lithos, str) ):
                lithos = geodyn1d.lithosphere(lithos)
    rho_ref_value, rho_ref_std, rho_ref = compute_rho_ref_perplex(perplex,domainBased=domainBased,lithos=lithos)
    perplex.deltarho                    = rho_ref - rho_ref_value
    perplex.rho_ref                     = rho_ref_value
    perplex.rho_ref_std                 = rho_ref_std
    
    return perplex

def compute_delta_rho_vertical_gradient_perplex(perplex):
    
    deltarho_vert_gradient = np.zeros_like(perplex.rho).reshape(perplex.np,perplex.nt)
    deltarho               = deepcopy(perplex.deltarho).reshape(perplex.np,perplex.nt)
    P                      = deepcopy(perplex.P).reshape(perplex.np,perplex.nt)
    P                      = P / np.max(P)
    deltarho               = deltarho / np.max(deltarho)
    
    for i in range(0,perplex.nt):
        #dP               = np.gradient(P[:,i])
        deltarho_along_P = deltarho[:,i]
        deltarho_vert_gradient[:,i] = np.gradient(deltarho_along_P, P[:,i])
    
    return deltarho_vert_gradient

def extrapolate_lowP_perplex(perplex):
    """ Extrapolate low pressure values that have not been computed.
    The assumption is that alpha and beta are constant and pressure-independent
    at low pressure which is certainly wrong. """
    
    extrapolated_perplex = perplex_output()
    extrapolated_perplex = deepcopy(perplex)
    
    Tref = 273
    Pref = 0
    
    rho   = extrapolated_perplex.rho.reshape(  int(perplex.np), int(perplex.nt))
    T     = extrapolated_perplex.T.reshape(    int(perplex.np), int(perplex.nt))
    P     = extrapolated_perplex.P.reshape(    int(perplex.np), int(perplex.nt))
    alpha = extrapolated_perplex.alpha.reshape(int(perplex.np), int(perplex.nt))
    beta  = extrapolated_perplex.beta.reshape( int(perplex.np), int(perplex.nt))
    cp    = extrapolated_perplex.cp.reshape( int(perplex.np), int(perplex.nt))
    melt  = extrapolated_perplex.melt.reshape( int(perplex.np), int(perplex.nt))
    rhomelt   = extrapolated_perplex.rhomelt.reshape(  int(perplex.np), int(perplex.nt))
    rhoresidual   = extrapolated_perplex.rhoresidual.reshape(  int(perplex.np), int(perplex.nt))
    alphamelt = extrapolated_perplex.alphamelt.reshape(int(perplex.np), int(perplex.nt))
    alpharesidual = extrapolated_perplex.alpharesidual.reshape(int(perplex.np), int(perplex.nt))
    betamelt      = extrapolated_perplex.betamelt.reshape( int(perplex.np), int(perplex.nt))
    betaresidual  = extrapolated_perplex.betaresidual.reshape( int(perplex.np), int(perplex.nt))
    cpmelt      = extrapolated_perplex.cpmelt.reshape( int(perplex.np), int(perplex.nt))
    cpresidual  = extrapolated_perplex.cpresidual.reshape( int(perplex.np), int(perplex.nt))
    
    for i in range(0,int(perplex.nt)):
        are_values = False
        for j in range(0,int(perplex.np)):
             if ( np.isnan(rho[j,i]) ):
                 #print('None T {} P {}'.format(T[j,i],P[j,i]))
                 pass
             else:
                 alpha_used = alpha[j,i]
                 beta_used  = beta[j,i]
                 cp_used    = cp[j,i]
                 rho_ref    = rho[j,i] / ( (1+beta_used*(P[j,i]-Pref)) * (1-alpha_used*(T[j,i]-Tref)) )
                 melt_ref   = melt[j,i]
                 
                 alpha_used_melt = alphamelt[j,i]
                 beta_used_melt  = betamelt[j,i]
                 cp_used_melt    = cpmelt[j,i]
                 rho_ref_melt    = rhomelt[j,i] / ( (1+beta_used_melt*(P[j,i]-Pref)) * (1-alpha_used_melt*(T[j,i]-Tref)) )
                 
                 alpha_used_residual = alpharesidual[j,i]
                 beta_used_residual  = betaresidual[j,i]
                 cp_used_residual    = cpresidual[j,i]
                 rho_ref_residual    = rhoresidual[j,i] / ( (1+beta_used_residual*(P[j,i]-Pref)) * (1-alpha_used_residual*(T[j,i]-Tref)) )
                 
                 #print('{} {} {} {} {} {}'.format(T[j,i],P[j,i],alpha_used,beta_used,rho_ref,rho[j,i]))
                 are_values  = True
                 break
        if (are_values):
            for j in range(0,int(perplex.np)):
                 if ( np.isnan(rho[j,i]) ):
                     rho[j,i]   = rho_ref*(1+beta_used*(P[j,i]-Pref))*(1-alpha_used*(T[j,i]-Tref))
                     melt[j,i]  = melt_ref
                     alpha[j,i] = alpha_used
                     beta[j,i]  = beta_used
                     cp[j,i]  = cp_used
                     
                     rhomelt[j,i]   = rho_ref_melt*(1+beta_used_melt*(P[j,i]-Pref))*(1-alpha_used_melt*(T[j,i]-Tref))
                     alphamelt[j,i] = alpha_used_melt
                     betamelt[j,i]  = beta_used_melt
                     cpmelt[j,i]    = cp_used_melt
                     
                     rhoresidual[j,i]   = rho_ref_residual*(1+beta_used_residual*(P[j,i]-Pref))*(1-alpha_used_residual*(T[j,i]-Tref))
                     alpharesidual[j,i] = alpha_used_residual
                     betaresidual[j,i]  = beta_used_residual
                     cpresidual[j,i]    = cp_used_residual
                 else:
                     break

    return extrapolated_perplex

def extrapolate_highT_perplex(perplex):
    """ Extrapolate high temperature values that have not been computed.
    compute. The assumption is that alpha and beta are constant and temperature-independent. """
    
    extrapolated_perplex = perplex_output()
    extrapolated_perplex = deepcopy(perplex)
    
    Tref = 273
    Pref = 0
    
    rho   = extrapolated_perplex.rho.reshape(  int(perplex.np), int(perplex.nt))
    T     = extrapolated_perplex.T.reshape(    int(perplex.np), int(perplex.nt))
    P     = extrapolated_perplex.P.reshape(    int(perplex.np), int(perplex.nt))
    alpha = extrapolated_perplex.alpha.reshape(int(perplex.np), int(perplex.nt))
    beta  = extrapolated_perplex.beta.reshape( int(perplex.np), int(perplex.nt))
    cp    = extrapolated_perplex.cp.reshape( int(perplex.np), int(perplex.nt))
    melt  = extrapolated_perplex.melt.reshape( int(perplex.np), int(perplex.nt))
    rhomelt       = extrapolated_perplex.rhomelt.reshape(  int(perplex.np), int(perplex.nt))
    rhoresidual   = extrapolated_perplex.rhoresidual.reshape(  int(perplex.np), int(perplex.nt))
    alphamelt     = extrapolated_perplex.alphamelt.reshape(int(perplex.np), int(perplex.nt))
    alpharesidual = extrapolated_perplex.alpharesidual.reshape(int(perplex.np), int(perplex.nt))
    betamelt      = extrapolated_perplex.betamelt.reshape( int(perplex.np), int(perplex.nt))
    betaresidual  = extrapolated_perplex.betaresidual.reshape( int(perplex.np), int(perplex.nt))
    cpmelt        = extrapolated_perplex.cpmelt.reshape( int(perplex.np), int(perplex.nt))
    cpresidual    = extrapolated_perplex.cpresidual.reshape( int(perplex.np), int(perplex.nt))
    
    for j in range(0,int(perplex.np)):
        are_values       = False
        residual_values  = False
        betaresidual_values = False
        cpresidual_values   = False
        rhoresidual_values  = False
        alphamelt_values = False
        melt_values      = False
        rho_values       = False
        for i in range(int(perplex.nt)-1,-1,-1):
            
            if ( (~np.isnan(alphamelt[j,i]))&(alphamelt_values==False) ):
                alpha_used_melt = alphamelt[j,i]
                alphamelt_values = True
            if ( (~np.isnan(rhomelt[j,i]))&(melt_values==False) ):
                beta_used_melt  = betamelt[j,i]
                melt_values     = True
            if ( (alphamelt_values==True)&(melt_values==True)):
                rho_ref_melt    = rhomelt[j,i] / ( (1+beta_used_melt*(P[j,i]-Pref)) * (1-alpha_used_melt*(T[j,i]-Tref)) )
            
            if ( (~np.isnan(alpharesidual[j,i]))&(residual_values==False) ):                
                alpha_used_residual = alpharesidual[j,i]
                residual_values     = True
            if ( (~np.isnan(betaresidual[j,i]))&(betaresidual_values==False) ):                
                beta_used_residual  = betaresidual[j,i]
                betaresidual_values = True
            if ( (~np.isnan(cpresidual[j,i]))&(cpresidual_values==False) ):                
                cp_used_residual  = cpresidual[j,i]
                cpresidual_values = True
            if ( (~np.isnan(rhoresidual[j,i]))&(rhoresidual_values == False) ):
                rhoresidual_values = True
                
            if ( (betaresidual_values==True)&(residual_values==True)&(rhoresidual_values==True)):
                rho_ref_residual    = rhoresidual[j,i] / ( (1+beta_used_residual*(P[j,i]-Pref)) * (1-alpha_used_residual*(T[j,i]-Tref)) )
                 
            if ( (~np.isnan(rho[j,i]))&(rho_values == False) ):
                alpha_used = alpha[j,i]
                beta_used  = beta[j,i]
                cp_used    = cp[j,i]
                rho_ref    = rho[j,i] / ( (1+beta_used*(P[j,i]-Pref)) * (1-alpha_used*(T[j,i]-Tref)) )
                melt_ref   = melt[j,i]
                rho_values  = True
        
            if ( (residual_values == True)&(melt_values == True)&(alphamelt_values == True)&(rho_values==True)&(betaresidual_values == True)&
                 (rhoresidual_values==True)&(cpresidual_values==True)):
                are_values  = True
                break
        
        if (are_values):
            residual_values  = False
            betaresidual_values = False
            cpresidual_values   = False
            rhoresidual_values  = False
            alphamelt_values = False
            melt_values      = False
            rho_values       = False
            for i in range(int(perplex.nt)-1,-1,-1):
                if ( np.isnan(rho[j,i]) ):
                    rho[j,i]   = rho_ref*(1+beta_used*(P[j,i]-Pref))*(1-alpha_used*(T[j,i]-Tref))
                    melt[j,i]  = melt_ref
                    alpha[j,i] = alpha_used
                    beta[j,i]  = beta_used
                else:
                    rho_values = True
                
                if ( np.isnan(rhomelt[j,i]) ):
                    rhomelt[j,i]   = rho_ref_melt*(1+beta_used_melt*(P[j,i]-Pref))*(1-alpha_used_melt*(T[j,i]-Tref))
                    betamelt[j,i]  = beta_used_melt
                else:
                    melt_values = True
                    
                if ( np.isnan(alphamelt[j,i]) ):
                    alphamelt[j,i] = alpha_used_melt
                    rhomelt[j,i]   = rho_ref_melt*(1+beta_used_melt*(P[j,i]-Pref))*(1-alpha_used_melt*(T[j,i]-Tref))
                else:
                    alphamelt_values = True
                
                if ( np.isnan(betaresidual[j,i]) ):  
                    betaresidual[j,i]  = beta_used_residual
                    rhoresidual[j,i]   = rho_ref_residual*(1+beta_used_residual*(P[j,i]-Pref))*(1-alpha_used_residual*(T[j,i]-Tref))
                else:
                    betaresidual_values = True
                    
                if ( np.isnan(cpresidual[j,i]) ):  
                    cpresidual[j,i]   = cp_used_residual
                else:
                    cpresidual_values = True
                    
                if ( np.isnan(rhoresidual[j,i]) ):
                    rhoresidual[j,i]   = rho_ref_residual*(1+beta_used_residual*(P[j,i]-Pref))*(1-alpha_used_residual*(T[j,i]-Tref))
                else:
                    rhoresidual_values = True
                
                if ( np.isnan(alpharesidual[j,i]) ):    
                    rhoresidual[j,i]   = rho_ref_residual*(1+beta_used_residual*(P[j,i]-Pref))*(1-alpha_used_residual*(T[j,i]-Tref))
                    alpharesidual[j,i] = alpha_used_residual
                else:
                    residual_values = True
                    
                if ((residual_values == True)&(melt_values == True)&(alphamelt_values == True)&(rho_values==True)&(betaresidual_values == True)&(rhoresidual_values == True)):
                    break
                

    return extrapolated_perplex

def extrapolate_highP_perplex(perplex,n_smooth = 2):
    """ Extrapolate high pressure values that have not been computed.
    The assumption is that alpha and beta are constant and pressure-independent 
    at high pressure which is wrong"""
    
    extrapolated_perplex = perplex_output()
    extrapolated_perplex = deepcopy(perplex)
    
    Tref = 273
    Pref = 0
    
    rho   = extrapolated_perplex.rho.reshape(  int(perplex.np), int(perplex.nt))
    T     = extrapolated_perplex.T.reshape(    int(perplex.np), int(perplex.nt))
    P     = extrapolated_perplex.P.reshape(    int(perplex.np), int(perplex.nt))
    alpha = extrapolated_perplex.alpha.reshape(int(perplex.np), int(perplex.nt))
    beta  = extrapolated_perplex.beta.reshape( int(perplex.np), int(perplex.nt))
    cp    = extrapolated_perplex.cp.reshape( int(perplex.np), int(perplex.nt))
    melt  = extrapolated_perplex.melt.reshape( int(perplex.np), int(perplex.nt))
    rhomelt   = extrapolated_perplex.rhomelt.reshape(  int(perplex.np), int(perplex.nt))
    rhoresidual   = extrapolated_perplex.rhoresidual.reshape(  int(perplex.np), int(perplex.nt))
    alphamelt = extrapolated_perplex.alphamelt.reshape(int(perplex.np), int(perplex.nt))
    alpharesidual = extrapolated_perplex.alpharesidual.reshape(int(perplex.np), int(perplex.nt))
    betamelt  = extrapolated_perplex.betamelt.reshape( int(perplex.np), int(perplex.nt))
    betaresidual  = extrapolated_perplex.betaresidual.reshape( int(perplex.np), int(perplex.nt))
    cpmelt      = extrapolated_perplex.cpmelt.reshape( int(perplex.np), int(perplex.nt))
    cpresidual  = extrapolated_perplex.cpresidual.reshape( int(perplex.np), int(perplex.nt))
    
    # smoothing alpha and beta along the boundaries to avoid vertical discontinuities not suitable
    n_smooth = n_smooth
    rho_smooth           = []
    rhomelt_smooth       = []
    rhoresidual_smooth   = []
    alpha_smooth         = []
    beta_smooth          = []
    cp_smooth            = []
    alphamelt_smooth     = []
    betamelt_smooth      = []
    alpharesidual_smooth = []
    betaresidual_smooth  = []
    cpmelt_smooth        = []
    cpresidual_smooth    = []  
    
    i_smooth = 0
    i_int    = 0
    for i in range(0,int(perplex.nt)):
        are_values = False
        for j in range(int(perplex.np)-1,-1,-1):
             if ( np.isnan(rho[j,i]) ):
                 #print('None T {} P {}'.format(T[j,i],P[j,i]))
                 pass
             else:
                 if (i_smooth<n_smooth):
                     alpha_smooth.append(alpha[j,i])
                     beta_smooth.append(beta[j,i])
                     cp_smooth.append(cp[j,i])
                     alphamelt_smooth.append(alphamelt[j,i])
                     betamelt_smooth.append(betamelt[j,i])
                     alpharesidual_smooth.append(alpharesidual[j,i])
                     betaresidual_smooth.append(betaresidual[j,i])
                     cpmelt_smooth.append(cpmelt[j,i])
                     cpresidual_smooth.append(cpresidual[j,i])
                     rho_smooth.append(rho[j,i])
                     rhomelt_smooth.append(rhomelt[j,i])
                     rhoresidual_smooth.append(rhoresidual[j,i])
                     i_smooth = i_smooth + 1
                 else:
                     alpha_smooth[i_int]         = alpha[j,i]
                     beta_smooth[i_int]          = beta[j,i]
                     cp_smooth[i_int]            = cp[j,i]
                     alphamelt_smooth[i_int]     = alphamelt[j,i]
                     betamelt_smooth[i_int]      = betamelt[j,i]
                     cpmelt_smooth[i_int]        = cpmelt[j,i]
                     alpharesidual_smooth[i_int] = alpharesidual[j,i]
                     betaresidual_smooth[i_int]  = betaresidual[j,i]
                     cpresidual_smooth[i_int]    = cpresidual[j,i]
                     rho_smooth[i_int]           = rho[j,i]
                     rhomelt_smooth[i_int]       = rhomelt[j,i]
                     rhoresidual_smooth[i_int]   = rhoresidual[j,i]
                     i_int = i_int + 1
                     if (i_int>=n_smooth):
                         i_int = 0
                
                 alpha_used = sum(alpha_smooth)/len(alpha_smooth)
                 beta_used  = sum(beta_smooth)/len(beta_smooth)
                 cp_used    = sum(cp_smooth)/len(cp_smooth)
                 rho_ref    = sum(rho_smooth)/len(rho_smooth) / ( (1+beta_used*(P[j,i]-Pref)) * (1-alpha_used*(T[j,i]-Tref)) )
                 
                 alpha_used_melt = sum(alphamelt_smooth)/len(alphamelt_smooth)
                 beta_used_melt  = sum(betamelt_smooth)/len(betamelt_smooth)
                 cp_used_melt    = sum(cpmelt_smooth)/len(cpmelt_smooth)
                 rho_ref_melt    = sum(rhomelt_smooth)/len(rhomelt_smooth) / ( (1+beta_used_melt*(P[j,i]-Pref)) * (1-alpha_used_melt*(T[j,i]-Tref)) )
                 
                 alpha_used_residual = sum(alpharesidual_smooth)/len(alpharesidual_smooth)
                 beta_used_residual  = sum(betaresidual_smooth)/len(betaresidual_smooth)
                 cp_used_residual  = sum(cpresidual_smooth)/len(cpresidual_smooth)
                 rho_ref_residual    = sum(rhoresidual_smooth)/len(rhoresidual_smooth) / ( (1+beta_used_residual*(P[j,i]-Pref)) * (1-alpha_used_residual*(T[j,i]-Tref)) )
                 #print('{} {} {} {} {} {}'.format(T[j,i],P[j,i],alpha_used,beta_used,rho_ref,rho[j,i]))
                 are_values  = True
                 break
        if (are_values):
            for j in range(int(perplex.np)-1,-1,-1):
                 if ( np.isnan(rho[j,i]) ):
                     rho[j,i]   = rho_ref*(1+beta_used*(P[j,i]-Pref))*(1-alpha_used*(T[j,i]-Tref))
                     melt[j,i]  = 0.0e0
                     alpha[j,i] = alpha_used
                     beta[j,i]  = beta_used
                     cp[j,i]    = cp_used
                     
                     rhomelt[j,i]   = rho_ref_melt*(1+beta_used_melt*(P[j,i]-Pref))*(1-alpha_used_melt*(T[j,i]-Tref))
                     alphamelt[j,i] = alpha_used_melt
                     betamelt[j,i]  = beta_used_melt
                     cpmelt[j,i]    = cp_used_melt
                     
                     rhoresidual[j,i]   = rho_ref_residual*(1+beta_used_residual*(P[j,i]-Pref))*(1-alpha_used_residual*(T[j,i]-Tref))
                     alpharesidual[j,i] = alpha_used_residual
                     betaresidual[j,i]  = beta_used_residual
                     cpresidual[j,i]    = cp_used_residual
                 else:
                     break

    return extrapolated_perplex

def bin2csv_tectmod(perplex,tablename):
    # Units bar,ºC
    T   = np.linspace(min(perplex.T),max(perplex.T),perplex.nt) - 273
    P   = np.linspace(min(perplex.P),max(perplex.P),perplex.np)/1e5
    rho = perplex.rho.reshape(int(perplex.np), int(perplex.nt))
    print('Pmin {} dP {} Tmin {} dT {} '.format(np.min(perplex.P),perplex.dP,np.min(perplex.T),perplex.dT))
    print('Pmin {} Pmax {} Tmin {} Tmax {} '.format(min(perplex.P),max(perplex.P),min(perplex.T),max(perplex.T)))
    # Write ascii file
    path     = os.path.abspath(__file__)
    filename = path.partition('thermodyn.py')[0]+'data/'+tablename+'.csv'

    if not (os.path.isfile(filename)):
        sys.exit("ERROR (bin2csv_tectmod): Either this file does not exist or the path is uncorrect {} ".format(filename))

    f = open(filename, 'w')
    f.write('nan')
    for i in range(0,int(perplex.nt)-1):
        f.write(',')
        f.write(str(T[i]))
    f.write('\n')
    for j in range(0,int(perplex.np)-1):
        f.write(str(P[j]))
        for i in range(0,int(perplex.nt)-1):
             f.write(',')
             f.write(str(rho[j,i]))
        f.write('\n')
    f.close()
    return

def write_binary_fortran_file_perplex(perplex,tablename,path=None):
    if path==None:
        path     = os.path.abspath(__file__)
        filename = path.partition('thermodyn.py')[0]+'data/'+tablename+'_fortran.bin'
    else:
        filename = path+'/'+tablename+'_fortran.bin'

    f = FortranFile(filename, 'w')

    # binary file with 6 columns and one header
    # header = ndata nt np dT dP
    # and then
    # T P alpha beta rho melt
    ndata = perplex.nt * perplex.np
    array_header     = np.zeros((1,5))
    array_header[0,0]  = ndata
    array_header[0,1]  = deepcopy(perplex.nt)
    array_header[0,2]  = deepcopy(perplex.np)
    array_header[0,3]  = deepcopy(perplex.dT)
    array_header[0,4]  = deepcopy(perplex.dP)
    f.write_record(np.array(np.transpose(array_header),dtype=np.float64))
    array_perplex      = np.zeros((ndata,12))
    array_perplex[:,0] = deepcopy(perplex.T)
    array_perplex[:,1] = deepcopy(perplex.P)
    array_perplex[:,2] = deepcopy(perplex.alpha)
    array_perplex[:,3] = deepcopy(perplex.beta)
    array_perplex[:,4] = deepcopy(perplex.rho)
    array_perplex[:,5] = deepcopy(perplex.melt)
    array_perplex[:,6] = deepcopy(perplex.rhomelt)
    array_perplex[:,7] = deepcopy(perplex.rhoresidual)
    array_perplex[:,8] = deepcopy(perplex.alphamelt)
    array_perplex[:,9] = deepcopy(perplex.alpharesidual)
    array_perplex[:,10] = deepcopy(perplex.betamelt)
    array_perplex[:,11] = deepcopy(perplex.betaresidual)
    f.write_record(np.array(np.transpose(array_perplex),dtype=np.float64))
    f.close()
    return



def write_binary_file_perplex(perplex,tablename,path=None):
    if ( path == None ):
        path     = os.path.abspath(__file__)
        filename = path.partition('thermodyn.py')[0]+'data/'+tablename+'.bin' 
    else:
        filename = path+'/'+tablename+'.bin'
    binary_file=open(filename,mode="wb+")
    # binary file with 6 columns and one header
    # header = ndata nt np dT dP
    # and then
    # T P alpha beta rho melt
    ndata = perplex.nt * perplex.np
    array_header     = np.zeros((1,5))
    array_header[0,0]  = ndata
    array_header[0,1]  = deepcopy(perplex.nt)
    array_header[0,2]  = deepcopy(perplex.np)
    array_header[0,3]  = deepcopy(perplex.dT)
    array_header[0,4]  = deepcopy(perplex.dP)
    np.save(binary_file, array_header)
    array_perplex      = np.zeros((ndata,15))
    array_perplex[:,0] = deepcopy(perplex.T)
    array_perplex[:,1] = deepcopy(perplex.P)
    array_perplex[:,2] = deepcopy(perplex.alpha)
    array_perplex[:,3] = deepcopy(perplex.beta)
    array_perplex[:,4] = deepcopy(perplex.rho)
    array_perplex[:,5] = deepcopy(perplex.melt)
    array_perplex[:,6] = deepcopy(perplex.rhomelt)
    array_perplex[:,7] = deepcopy(perplex.rhoresidual)
    array_perplex[:,8] = deepcopy(perplex.alphamelt)
    array_perplex[:,9] = deepcopy(perplex.alpharesidual)
    array_perplex[:,10] = deepcopy(perplex.betamelt)
    array_perplex[:,11] = deepcopy(perplex.betaresidual)
    array_perplex[:,12] = deepcopy(perplex.cp)
    array_perplex[:,13] = deepcopy(perplex.cpmelt)
    array_perplex[:,14] = deepcopy(perplex.cpresidual)

    #print('write_binary_file_perplex ............')
    #print(perplex.cp)

    np.save(binary_file, array_perplex)
    binary_file.close()
    return

def read_binary_file_perplex(tablename,path=None):
    if path==None:
        path     = os.path.abspath(__file__)
        filename = path.partition('thermodyn.py')[0]+'data/'+tablename+'.bin'   
    else:
        filename = path+'/'+tablename+'.bin'

    if not (os.path.isfile(filename)):
        sys.exit("ERROR(read_binary_file_perplex): Either this file does not exist (maybe need to compute the grid using the notebook thermo.ipynb) or the path is uncorrect {} ".format(filename))

    perplex  = perplex_output()
    
    # we read the file
    binary_file=open(filename,mode="rb+")
    array_header  = np.load(binary_file)
    array_perplex = np.load(binary_file)
    binary_file.close()
    # we copy data to perplex
    perplex.nt = int(array_header[0,1])
    perplex.np = int(array_header[0,2])
    perplex.dT = array_header[0,3]
    perplex.dP = array_header[0,4]
    perplex.T  = array_perplex[:,0]
    perplex.P  = array_perplex[:,1]
    perplex.alpha = array_perplex[:,2]
    perplex.beta  = array_perplex[:,3]
    perplex.rho   = array_perplex[:,4]
    perplex.melt  = array_perplex[:,5]
    perplex.rhomelt       = array_perplex[:,6]
    perplex.rhoresidual   = array_perplex[:,7]
    perplex.alphamelt     = array_perplex[:,8]
    perplex.alpharesidual = array_perplex[:,9]
    perplex.betamelt      = array_perplex[:,10]
    perplex.betaresidual  = array_perplex[:,11]

    perplex.cp            = array_perplex[:,12]
    perplex.cpmelt        = array_perplex[:,13]
    perplex.cpresidual    = array_perplex[:,14]


    #print('read_binary_file_perplex ............')
    #print(perplex.cp)

    #compute_delta_rho_perplex(perplex)
    perplex.deltarho   = np.zeros_like(perplex.rho)
    perplex.rho_ref     = 0.0e0
    perplex.rho_ref_std = 0.0e0

    Tref  = 273
    Pref  = 0

    factor1        = 1 + np.multiply(perplex.beta,perplex.P-Pref)
    factor2        = 1 - np.multiply(perplex.alpha,perplex.T-Tref)
    factor         = np.multiply(factor1,factor2)
    perplex.rhoref = np.divide(perplex.rho,factor)

    return perplex

def read_matlab_file(tablename,nina_table=True):
    path     = os.path.abspath(__file__)
    filename = path.partition('thermodyn.py')[0]+'data/'+tablename+'.mat' 
    
    if (nina_table):
        perplex  = perplex_output()
        matfile_list = sio.loadmat(filename)
        P   = matfile_list['P_nina']
        T   = matfile_list['T_nina']
        rho = matfile_list['rho_nina']
        minT = np.min(T)
        maxT = np.max(T)
        minP = np.min(P)
        maxP = np.max(P)
        
        perplex.nt = int(T.shape[1])
        perplex.np = int(P.shape[1])
        ndata      = perplex.nt * perplex.np
        perplex.dT = (maxT-minT)/(perplex.nt-1)
        perplex.dP = (maxP-minP)/(perplex.np-1)
            
        #newT           = np.arange(minT,maxT+perplex.dT,perplex.dT)
        newT           = np.linspace(minT,maxT,perplex.nt)
        #newP           = np.arange(minP,maxP+perplex.dP,perplex.dP)
        newP           = np.linspace(minP,maxP,perplex.np)
        tt, pp         = np.meshgrid(newT, newP)
        perplex.T      = tt.ravel()
        perplex.P      = pp.ravel()
        perplex.rho    = np.transpose(rho).reshape(ndata)
        
        perplex.melt   = np.zeros_like(perplex.rho)
        perplex.alpha  = np.zeros_like(perplex.rho)
        perplex.beta   = np.zeros_like(perplex.rho)
        perplex.rhomelt       = np.zeros_like(perplex.rho)
        perplex.rhoresidual   = np.zeros_like(perplex.rho)
        perplex.alphamelt     = np.zeros_like(perplex.rho)
        perplex.alpharesidual = np.zeros_like(perplex.rho)
        perplex.betamelt      = np.zeros_like(perplex.rho)
        perplex.betaresidual  = np.zeros_like(perplex.rho)
        perplex.deltarho      = np.zeros_like(perplex.rho)
        perplex.rho_ref     = 0.0e0
        perplex.rho_ref_std = 0.0e0
        Tref  = 273
        Pref  = 0
    
        factor1        = 1 + np.multiply(perplex.beta,perplex.P-Pref)
        factor2        = 1 - np.multiply(perplex.alpha,perplex.T-Tref)
        factor         = np.multiply(factor1,factor2)
        perplex.rhoref = np.divide(perplex.rho,factor)

        # We define constant alpha and beta requested for extrapolation (no data from nina, check with Ritske?)
        # We use usual alpha and beta that are certainly a good approximation along continental geotherm
        # at thermal steady-state. Indeed, geotherms are almost parallel to alpha and beta changes.
        # However, these values are wrong for oceanic domains where temperature gradient are less important.
        # The pressure dependency of alpha,beta and then rho at high temperature is important.
        # 
        perplex.alpha.fill(3.1e-5)
        perplex.beta.fill(0.85e-11)
        
        # nina table------------------
        # __header__
        # __version__
        # __globals__
        # ph_1_T_nina
        # ph_1_P_nina
        # ph_1_rho_nina
        # P_nina
        # T_nina
        # rho_nina
        # ph_2_T_nina
        # ph_2_P_nina
        # ph_2_rho_nina
    else:
        print('NOT IMPLEMENTED. "read_matlab_file" is only implemented for Nina s files')
        return
    return perplex
