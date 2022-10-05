import numpy as np
import sys
import os
from ..tools import *

dz_init = 200 # 10 / 200  for fast low resolution perplex # vertical sampling (m)

class temperature_profile_object(object):
    
    def __init__(self,T=[],Z=[],k=[],dz=dz_init):
        self.T  = T   # temperature (ºK)
        self.Z  = Z   # depth (km)
        self.k  = k   # thermal conductivity profile
        self.dz = dz  # 10 / 200  for fast low resolution perplex # output sampling (m)


def read_imposed_geotherm(layers,filename,printOut=True,dz=dz_init):

    geotherm         = temperature_profile_object(dz=dz)
    nlayers          = len(layers)
    Column_thickness = sum([ layers[i].thickness for i in range(0,nlayers) ])

    #Read file -------------------------
    if (os.path.isfile(filename)):
        data = np.loadtxt(filename,dtype='float',delimiter=' ')
    else:
        print('Problem with unread file '+filename)
        quit()
    z_read = data[:,1]
    t_read = data[:,0]
    
    #Interpolate read values -----------
    #Initialization
    dz     = geotherm.dz # resolution for calculations (m)
    i      = 0
    T      = []
    Z      = []
    k_arr  = []
    Tbasal = layers[nlayers-1].thermalBc.temp_bottom
    k      = thermal_conductivity().get_value(depth=Column_thickness,layers=layers,temperature=Tbasal)
    T.append(Tbasal)
    k_arr.append(k)
    Z.append(Column_thickness)

    temp = Tbasal
    for zincr in np.arange(dz,Column_thickness+dz,dz):
        i = i + 1
        z    = Column_thickness - zincr
        k    = thermal_conductivity().get_value(depth=z,layers=layers,temperature=temp)
        #print('{0:8.2f} {1:8.2f}'.format(z,z_read[len(z_read)-1]))
        if ( z > z_read[len(z_read)-1] ):
            temp = t_read[len(z_read)-1]
        else:
            temp = pos2value_1d(z_read,t_read,z)
        #print('{0:8.2f} {1:8.2f}'.format(z,temp))
        k_arr.append(k)
        Z.append(z)
        T.append(temp)
    
    geotherm.T = np.asarray(T)+273.15
    geotherm.Z = np.asarray(Z)
    geotherm.k = np.asarray(k_arr)
    return geotherm

def apply_temperature_anomaly(Z,T,dT,mohodepth=35e3,labdepth=125e3):
#    i//dtAno = thermalparams->AnoDT * (1-fabs((thermalparams->AnoXcenter-x)/(thermalparams->AnoXwidth/2.0e0)))
#        //                             * (1-fabs((thermalparams->AnoYcenter-y)/(thermalparams->AnoYwidth/2.0e0)));
#            dtAno = thermalparams->AnoDT * (1-fabs((thermalparams->AnoXcenter-x)/(thermalparams->AnoXwidth/2.0e0)))
#                                             * (1-fabs((thermalparams->AnoYcenter-y)/thermalparams->AnoYwidth));
    zmax = Z.max()
    i    = 0
    newT = []
    for z in Z:
        dTAno = 0.0
        if (z>mohodepth and z <labdepth):
            dTAno = dT * (1-np.abs((labdepth-z)/(labdepth-mohodepth)))
        elif (z>=labdepth):
            dTAno = dT * (1-np.abs((labdepth-z)/(zmax-labdepth)))
        newT.append(T[i] + dTAno)
        i = i + 1
    return newT

def compute_steady_state_geotherm(layers,printOut=True,dz=dz_init):
    """ Compute the geotherm by iterative implicit finite difference solve
        looking for Qbasal that fits the temperatures at boundaries.
        The analytical solution is based on steady state equation type
        (constant thermal conductivity):
           T = -A/(2k) * dz^2 + Qtop/k * dz + Ttop
           Qbottom = Qtop - A * dz
           ( Chapman, 1986 )
        The solution is discretized in depth to deal with variable thermal
        conductivity.
        Layers represents a 1-D lithospheric structure.
        Each "layer" of "layers" is an object defined in core.py of geodyn1d."""

    geotherm = temperature_profile_object(dz=dz)
    nlayers  = len(layers)
   
    # Geotherm is computed from bottom
    # Qbasal (q_bottom) of the deepest layer is required
    if (layers[nlayers-1].thermalBc.q_bottom==None):
        layers[nlayers-1].thermalBc.q_bottom = 19.5e-3
        if ( printOut ):
            print('We start with a default value Qbasal = {0:7.5f}'.format(layers[nlayers-1].thermalBc.q_bottom))
    delta_q = 2.0e-3
    iloop   = 0
    
    # We look for the best Qbasal that fits
    # the Moho temperature.
    isnot_geotherm = True
    while (isnot_geotherm):
        iloop = iloop + 1
        isnot_geotherm = False
        dz     = geotherm.dz # resolution for calculations (m)
        T      = []
        Z      = []
        k_arr  = []
        i      = 0
        Tbasal = layers[nlayers-1].thermalBc.temp_bottom
        Tp     = layers[nlayers-1].thermalBc.temp_potential
        Qbasal = layers[nlayers-1].thermalBc.q_bottom
        Column_thickness = sum([ layers[i].thickness for i in range(0,nlayers) ])
        Grad_slm         = (Tbasal-Tp)/Column_thickness
        k                = thermal_conductivity().get_value(depth=Column_thickness,layers=layers,temperature=Tbasal)
        T.append(Tbasal)
        k_arr.append(k)
        Z.append(Column_thickness)
        
        Qdown = Qbasal;
        Tdown = Tbasal;
        
        for zincr in np.arange(dz,Column_thickness+dz,dz):
            i = i + 1
            z        = Column_thickness - zincr
            cumdepth = 0.0e0
            for ilayer in range(0,nlayers):
                #cumdepth = cumdepth + layers[ilayer].thickness
                if (z<cumdepth+ layers[ilayer].thickness):
                    break
                cumdepth = cumdepth + layers[ilayer].thickness
    
            A = layers[ilayer].material.thermalProperties.H
            k = thermal_conductivity().get_value(depth=z,layers=layers,temperature=T[i-1])
            
            k_arr.append(k)
            T.append(Tdown + (A/(2*k)) * dz**2 - (Qdown/k) * dz)
            Z.append(z)
            
            #print('z {} ilayer {} A {} k {} T {}'.format(z,ilayer,A,k,T[i]-273.15))
    
            Tdown = T[i]
            Qdown = Qdown + A*dz;
    
            #print('cumdepth {} z {} temp_bottom {} temp_top {}'.format(cumdepth,z,layers[ilayer].thermalBc.temp_bottom,layers[ilayer].thermalBc.temp_top))
    
            if (z<=cumdepth+dz/2):
                if ((layers[ilayer].thermalBc.temp_top!=None)&
                    (layers[ilayer].thermalBc.adiabat==False)&
                    (layers[ilayer].thickness!=0)):
                    thermalBc_Tdiff = Tdown-layers[ilayer].thermalBc.temp_top
                    if (abs(thermalBc_Tdiff)>0.1):
                        if (thermalBc_Tdiff>0):
                            if (delta_q>0):
                                delta_q = delta_q
                            else:
                                delta_q = -(1/2.0e0) * delta_q
                        else:
                            if (delta_q<0):
                                delta_q = delta_q
                            else:
                                delta_q = -(1/2.0e0) * delta_q
                        layers[nlayers-1].thermalBc.q_bottom = layers[nlayers-1].thermalBc.q_bottom + delta_q
                        isnot_geotherm = True
                        break

    if ( printOut ):
        print(' z bottom (km) = {0:6.1f}'.format(Z[0]))
        print(' k asth. = {0:5.2f}'.format(k_arr[0]))
        print(' Grad slm = {}'.format(Grad_slm*1e3))
        print(' Qbasal = {0:12.10f}'.format(Qbasal))
        print(' Tbasal = {0:6.1f}'.format(Tbasal-273.15))
        print(' Tp = {0:6.1f}'.format(Tp-273.15))
        i     = 0
        Qdown = Qbasal
        for zincr in np.arange(dz,Column_thickness+dz,dz):
            i = i + 1
            z = Column_thickness - zincr
            cumdepth = 0.0e0
            for ilayer in range(0,nlayers):
                #cumdepth = cumdepth + layers[ilayer].thickness
                if (z<cumdepth+ layers[ilayer].thickness):
                    break
                cumdepth = cumdepth + layers[ilayer].thickness

            A     = layers[ilayer].material.thermalProperties.H
            Qdown = Qdown + A*dz

            if (z<=cumdepth+dz/2):
                if ((layers[ilayer].thermalBc.temp_top!=None)&
                    (layers[ilayer].thermalBc.adiabat==False)&
                    (layers[ilayer].thickness!=0)):
                    #print('Tdown {} layers[ilayer].thermalBc.temp_top {}'.format(T[i],layers[ilayer].thermalBc.temp_top))
                    print('Layer {0:} T top = {1:6.1f} Q top = {2:5.2f} ({3:})'.format(ilayer,layers[ilayer].thermalBc.temp_top-273.15,Qdown*1e3,layers[ilayer].name_material))
                else:
                    print('Layer {0:} T top = {1:6.1f} Q top = {2:5.2f} ({3:})'.format(ilayer,T[i]-273.15,Qdown*1e3,layers[ilayer].name_material))
        print(' Qsurface = {0:5.2f} mW/m2'.format(Qdown*1e3))
    
    # output resolution
    # resampling at dz = geotherm.dz
                
    geotherm.T = np.asarray(T)
    geotherm.Z = np.asarray(Z)
    geotherm.k = np.asarray(k_arr) 
    return geotherm


class thermal_conductivity(object):
       
    def get_value(self,depth,layers,temperature):
        """Dispatch method"""
        nlayers  = len(layers)
        cumdepth = 0.0e0
        for ilayer in range(nlayers):
            if ((depth>=cumdepth) & (depth<cumdepth + layers[ilayer].thickness)):
                break
            cumdepth = cumdepth + layers[ilayer].thickness
        k_type = layers[ilayer].material.thermalProperties.k_type
        method_name = 'k_type' + str(k_type)
        # Get the method from 'self'. Default to a lambda.
        method = getattr(self, method_name, 0)
        return method(depth,layers,temperature)
    
    def k_type0(self,depth,layers,temperature):
        # Thermal conductivity constant
        # defined by layer properties
        nlayers  = len(layers)
        Column_thickness = sum([ layers[i].thickness for i in range(0,nlayers) ])
        cumdepth = 0.0e0
        for i in range(nlayers):
            if ((depth>=cumdepth) & (depth<cumdepth + layers[i].thickness)):
                break
            cumdepth = cumdepth + layers[i].thickness
        k = layers[i].material.thermalProperties.k
        if ( layers[i].thermalBc.adiabat ):
            Tbasal = layers[i].thermalBc.temp_bottom
            Tp     = layers[i].thermalBc.temp_potential
            Qbasal = layers[i].thermalBc.q_bottom
            if cumdepth == Column_thickness:
                Grad_slm = (Tbasal-Tp)/cumdepth
            else:
                Grad_slm = (Tbasal-Tp)/(cumdepth+layers[i].thickness)
            k        = Qbasal/Grad_slm
        return k
    def k_type1(self,depth,layers,temperature):
        # Thermal conductivity temperature dependent
        # defined by layer properties
        nlayers  = len(layers)
        Column_thickness = sum([ layers[i].thickness for i in range(0,nlayers) ])
        cumdepth = 0.0e0
        for i in range(nlayers):
            if ((depth>=cumdepth) & (depth<cumdepth + layers[i].thickness)):
                break
            cumdepth = cumdepth + layers[i].thickness
        k = layers[i].material.thermalProperties.k

        #if (depth>35e3):
        diffusivity_factor = 2.5 # mantle
        #else:
        #    diffusivity_factor = 1.5 # crust

        T1 = 273
        T2 = 876 + 273
        a  = (k-diffusivity_factor*k)/(T2-T1)
        b  = diffusivity_factor*k - a*T1
        if (temperature<T2):
            k = a*temperature + b

        if ( layers[i].thermalBc.adiabat ):
            Tbasal = layers[i].thermalBc.temp_bottom
            Tp     = layers[i].thermalBc.temp_potential
            Qbasal = layers[i].thermalBc.q_bottom
            if cumdepth == Column_thickness:
                Grad_slm = (Tbasal-Tp)/cumdepth
            else:
                Grad_slm = (Tbasal-Tp)/(cumdepth+layers[i].thickness)
            k        = Qbasal/Grad_slm
        return k
    
    def k_type2(self,depth,layers,temperature):
        # Thermal conductivity temperature dependent
        # defined by layer properties
        nlayers  = len(layers)
        Column_thickness = sum([ layers[i].thickness for i in range(0,nlayers) ])
        cumdepth = 0.0e0
        for i in range(nlayers):
            if ((depth>=cumdepth) & (depth<cumdepth + layers[i].thickness)):
                break
            cumdepth = cumdepth + layers[i].thickness
        k = layers[i].material.thermalProperties.k

        # Sekiguchi 1984
        km = 1.8418
        Tsurface = 5 # ºC
        Tm = 1473.15 + Tsurface
        T0 = 273.15 + Tsurface
        k  = km + ((Tm*T0)/(Tm-T0))*(k-km)*((1/temperature)-(1/Tm))

        if ( layers[i].thermalBc.adiabat ):
            Tbasal = layers[i].thermalBc.temp_bottom
            Tp     = layers[i].thermalBc.temp_potential
            Qbasal = layers[i].thermalBc.q_bottom
            if cumdepth == Column_thickness:
                Grad_slm = (Tbasal-Tp)/cumdepth
            else:
                Grad_slm = (Tbasal-Tp)/(cumdepth+layers[i].thickness)
            k        = Qbasal/Grad_slm
        return k


        
