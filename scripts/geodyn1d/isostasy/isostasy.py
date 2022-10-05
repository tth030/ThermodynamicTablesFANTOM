
import numpy as np
from ..tools import *
from ..thermodyn import *
from scipy.interpolate import RegularGridInterpolator

g = 9.80665 # (m/s2)

class pressure_density_profile(object):
    
    def __init__(self,P=[],Z=[],rho=[],melt_fraction=[],dz=200):
        self.P  = P     # pressure (MPa)
        self.Z  = Z     # depth (km)
        self.rho  = rho # density profile Kg/m3
        self.melt_fraction = melt_fraction # percentage
        self.dz = dz    # output sampling (m)

def read_imposed_density_profile(layers,filename,printOut=True):
    P_rho_DepthProfile = pressure_density_profile()
    nlayers            = len(layers)
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
    dz     = P_rho_DepthProfile.dz # resolution for calculations (m)
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
    
    geotherm.T = np.asarray(T)
    geotherm.Z = np.asarray(Z)
    geotherm.k = np.asarray(k_arr)
    return geotherm


def compute_pressure_density_profile(layers,geotherm,density_type=1,printOut=True,filename=None,drho=0,zdrho=[6500,125000],path=None):
    """ Compute the pressure profile based on properties
    given by the list of layers for a given temperature profile.
    Each "layer" is an object defined in core.py of geodyn1d """
    P_rho_DepthProfile = pressure_density_profile()
    nlayers  = len(layers)
    
    P             = np.zeros_like(geotherm.T).tolist()
    rho           = np.zeros_like(geotherm.T).tolist()
    melt_fraction = np.zeros_like(geotherm.T).tolist()
    T             = geotherm.T.tolist()
    Z             = geotherm.Z.tolist()
    
    nl  = len(Z)
    #print('min {} max {} len {}'.format(min(Z),max(Z),len(Z)))
    #print('min {} max {} len {}'.format(min(T),max(T),len(T)))

    if filename:
        if (os.path.isfile(filename)):
            data = np.loadtxt(filename,dtype='float',delimiter=' ')
        else:
            print('Problem with unread file '+filename)
            quit()
        z_read   = data[:,1]
        rho_read = data[:,0]
        #TODO melt_fraction_read = data[:,2]

    for ni in range(0,nl):
        z = Z[nl-ni-1]
        t = T[nl-ni-1]
        cumdepth = 0.0e0
        for ilayer in range(0,nlayers):
            #cumdepth = cumdepth + layers[ilayer].thickness
            if (z<cumdepth+ layers[ilayer].thickness):
                break
            cumdepth = cumdepth + layers[ilayer].thickness
        
        if (ni==0):
            dz       = 0.0e0
            p        = 0.0e0
            if filename:
                rho_used = rho_read[0]
                melt_fraction_used = 0 #TODO melt_fraction_read[0] need to be added inside the file
            else:
                rho_used = density().get_value(ilayer=ilayer,depth=z,pressure=p,temp=t,layers=layers,path=path)
                melt_fraction_used = melt_fraction_obj().get_value(ilayer=ilayer,depth=z,pressure=p,temp=t,layers=layers,path=path)
            if (z>zdrho[0] and z<zdrho[1]):
                rho_used = rho_used + drho
            rho[nl-ni-1] = rho_used
            melt_fraction[nl-ni-1] = melt_fraction_used
            P[nl-ni-1]   = p
        else:
            dz= (Z[nl-ni-1] - Z[nl-ni])
            p = P[nl-ni]
            if filename:
                if ( z > z_read[len(z_read)-1] ):
                    rho_used = rho_read[len(z_read)-1]
                else:
                    rho_used = pos2value_1d(z_read,rho_read,z)
                melt_fraction_used = 0 #TODO melt_fraction_read need to be added inside the file
            else:
                rho_used = density().get_value(ilayer=ilayer,depth=z,pressure=p,temp=t,layers=layers,path=path)
                melt_fraction_used = melt_fraction_obj().get_value(ilayer=ilayer,depth=z,pressure=p,temp=t,layers=layers,path=path)
            if (z>zdrho[0] and z<zdrho[1]):
                rho_used = rho_used + drho
            rho[nl-ni-1]       = rho_used
            melt_fraction[nl-ni-1] = melt_fraction_used
            p          = rho[nl-ni-1] * g * dz
            P[nl-ni-1] = P[nl-ni] + p
        #print('{} {} {}'.format(Z[nl-ni-1],z, rho[nl-ni-1]))
       
    P_rho_DepthProfile.P             = np.asarray(P)
    P_rho_DepthProfile.Z             = np.asarray(Z)
    P_rho_DepthProfile.rho           = np.asarray(rho)
    P_rho_DepthProfile.melt_fraction = np.asarray(melt_fraction)
    P_rho_DepthProfile.dz            = dz
    return P_rho_DepthProfile

#use_tidv triggers the use of thermally induced density variations
#use_tidv = 0 means no thermal expansion
#use_tidv = 1 means thermal expansion with only using coefficient thermal_expansion_i
#use_tidv = 2 means temperature dependent thermal expansion linearly scaled from thermal_expansion_1_i to thermal_expansion_2_i
#             within the temperature range thermal_expansion_t1_i to thermal_expansion_t2_i
#         = 3 temperature and pressure dependent thermal expansion
#         = 4 NOT IMPLEMENTED constant thermal expansion and compressibility
#         = 5 NOT IMPLEMENTED temp dep. thermal expansion and compressibility
#         = 6 NOT IMPLEMENTED temp. and P. dep. thermal expansion and compressibility
#         = 7 DENSITY READ FROM TABLE
#         = 8 DELTARHO that apply on ref. density READ FROM TABLE and constant thermal expansion similar to case "1"
#         = 9 DELTARHO that apply on ref. density READ FROM TABLE and thermal expansion(T) similar to case "2"
#         = 10 DELTARHO that apply on ref. density READ FROM TABLE and thermal expansion(T,P) similar to case "3"

class density(object):
       
    def get_value(self,ilayer,depth,pressure,temp,layers,path=None):
        """Dispatch method"""
        density_type=layers[ilayer].material.densityComposition.use_tidv
        method_name = 'density_type' + str(density_type)
        # Get the method from 'self'. Default to a lambda.
        method = getattr(self, method_name, 0)
        return method(ilayer,depth,pressure,temp,layers,path)
    
    def density_type0(self,ilayer,depth,pressure,temp,layers,path):
        # Density constant
        # defined by layer properties
        density = layers[ilayer].material.densityComposition.rho + layers[ilayer].material.densityComposition.deltarho
        return density

    def density_type1(self,ilayer,depth,pressure,temp,layers,path):
        # Density dependent of the thermal expansion only
        # defined by layer properties
        density0 = layers[ilayer].material.densityComposition.rho
        alpha    = thermal_expansion().get_value(ilayer=ilayer,depth=depth,pressure=pressure,temp=temp,layers=layers,alpha_type=1) 
        T0       = layers[ilayer].material.thermalProperties.T0 + 273
        density = density0 * ( 1 - alpha*(temp-T0)) + layers[ilayer].material.densityComposition.deltarho
        #print('density {} density0 {} alpha {} temp {} T0 {} depth {}'.format(density,density0,alpha,temp,T0,depth))
        return density
    
    def density_type2(self,ilayer,depth,pressure,temp,layers,path):
        # Density dependent of the thermal expansion and compressibility
        # defined by layer properties
        density0 = layers[ilayer].material.densityComposition.rho
        alpha1,alpha2,alphaT1,alphaT2    = thermal_expansion().get_value(ilayer=ilayer,depth=depth,pressure=pressure,temp=temp,layers=layers,alpha_type=2)
        T0       = layers[ilayer].material.thermalProperties.T0 + 273
        alphaT1  = alphaT1 + 273
        alphaT2  = alphaT2 + 273
        a        = (alpha2-alpha1)/(alphaT2-alphaT1)
        alpha    = a*temp + alpha1
        density = density0 * ( 1 - alpha*(temp-T0)) + layers[ilayer].material.densityComposition.deltarho
        return density

    def density_type7(self,ilayer,depth,pressure,temp,layers,path):
        # Density read in a table
        # defined by layer properties
        if layers[ilayer].material.densityComposition.perplex_name is None:
            print('You use use_tidv = 7 for this layer {} but perplex_name is None.'.format(ilayer)) 
            raise TypeError
        else:
            perplex_name = layers[ilayer].material.densityComposition.perplex_name
        
        if ( not layers[ilayer].material.densityComposition.perplex ):
            layers[ilayer].material.densityComposition.perplex = read_binary_file_perplex(perplex_name,path)
        
        interp_method = 0
        
        if (interp_method==0):
            pmax = np.max(layers[ilayer].material.densityComposition.perplex.P)
            tmax = np.max(layers[ilayer].material.densityComposition.perplex.T)
            if ( not layers[ilayer].f_interp):
                perplex_np = layers[ilayer].material.densityComposition.perplex.np
                perplex_nt = layers[ilayer].material.densityComposition.perplex.nt
                value = layers[ilayer].material.densityComposition.perplex.rho.reshape(perplex_np, perplex_nt)
                T = np.linspace(np.min(layers[ilayer].material.densityComposition.perplex.T),
                                tmax,
                                perplex_nt)
                P = np.linspace(np.min(layers[ilayer].material.densityComposition.perplex.P),
                                pmax,
                                perplex_np)            
                layers[ilayer].f_interp = RegularGridInterpolator((np.array(P)/pmax, np.array(T)/tmax), value)
            
            if ((pressure>pmax) | (temp>tmax)):
                print('WARNING: T {} ; Tmax {} P {} ; Pmax {}'.format(temp,tmax,pressure,pmax))
                if (pressure>pmax):
                    pressure = pmax
                if (temp>tmax):
                    temp = tmax
                density       = layers[ilayer].f_interp([pressure/pmax,temp/tmax])[0] + layers[ilayer].material.densityComposition.deltarho
                print('density {}'.format(density))
            else:
                density = layers[ilayer].f_interp([pressure/pmax,temp/tmax])[0] + layers[ilayer].material.densityComposition.deltarho

        elif (interp_method==1):
            # we use index to localize point inside the grid and we then interpolate between points
            pmin = np.min(layers[ilayer].material.densityComposition.perplex.P)
            tmin = np.min(layers[ilayer].material.densityComposition.perplex.T)
            # print('pmin {} tmin {} pressure {}'.format(pmin,tmin,pressure))
            ip = int(np.ceil((pressure-pmin)/layers[ilayer].material.densityComposition.perplex.dP))
            it = int(np.ceil((temp-tmin)/layers[ilayer].material.densityComposition.perplex.dT))
            #nt = layers[ilayer].material.densityComposition.perplex.nt
            perplex_np = layers[ilayer].material.densityComposition.perplex.np
            perplex_nt = layers[ilayer].material.densityComposition.perplex.nt
            # print('np {} ip {} it {} index {}/{}'.format(perplex_np,ip,it,perplex_nt*(ip-1)+it,perplex_np*perplex_nt))
            
            d1 = layers[ilayer].material.densityComposition.perplex.rho[perplex_nt*(ip-1)+it-1]
            p1 = layers[ilayer].material.densityComposition.perplex.P[perplex_nt*(ip-1)+it-1]
            t1 = layers[ilayer].material.densityComposition.perplex.T[perplex_nt*(ip-1)+it-1]
            
            d2 = layers[ilayer].material.densityComposition.perplex.rho[perplex_nt*(ip-1+1)+it-1]
            p2 = layers[ilayer].material.densityComposition.perplex.P[perplex_nt*(ip-1+1)+it-1]
            # t2 = layers[ilayer].material.densityComposition.perplex.T[perplex_nt*(ip-1+1)+it-1] # debugging
            
            d3 = layers[ilayer].material.densityComposition.perplex.rho[perplex_nt*(ip-1+1)+it-1+1]
            # p3 = layers[ilayer].material.densityComposition.perplex.P[perplex_nt*(ip-1+1)+it-1+1] # debugging
            # t3 = layers[ilayer].material.densityComposition.perplex.T[perplex_nt*(ip-1+1)+it-1+1] # debugging
            
            d4 = layers[ilayer].material.densityComposition.perplex.rho[perplex_nt*(ip-1)+it-1+1]
            # p4 = layers[ilayer].material.densityComposition.perplex.P[perplex_nt*(ip-1)+it-1+1] # debugging
            t4 = layers[ilayer].material.densityComposition.perplex.T[perplex_nt*(ip-1)+it-1+1]
            
            #               T
            #               ^
            #               |
            # 4 ------- 3   |
            # |         |   |
            # |         |   | 
            # |         |   |
            # 1---------2   |
            # -------------->P

            # print('DEBUG: T {} ; P {}'.format(temp,pressure))
            # print('1 {} {} {}'.format(t1,p1,d1))
            # print('2 {} {} {}'.format(t2,p2,d2))
            # print('3 {} {} {}'.format(t3,p3,d3))
            # print('4 {} {} {}'.format(t4,p4,d4))
            
            pp = (pressure-p1)/(p2-p1)
            tt = (temp-t1)/(t4-t1)
            density = d1*(1-tt)*(1-pp) + d4*tt*(1-pp) + d2*pp*(1-tt) + d3*pp*tt
            
            density = density + layers[ilayer].material.densityComposition.deltarho
            
            #print('density = {}'.format(density))

        return density

    def density_type8(self,ilayer,depth,pressure,temp,layers,path):
        # Density read in another file
        # Usually the file is from 2-D experiment
        if layers[ilayer].material.densityComposition.filename_densityprofile is None:
            print('You use use_tidv = 8 for this layer {} but filename_densityprofile is None.'.format(ilayer))
            raise TypeError
        else:
            filename_densityprofile = layers[ilayer].material.densityComposition.filename_densityprofile

        if ( layers[ilayer].material.densityComposition.densityprofile is None ):
            layers[ilayer].material.densityComposition.densityprofile = self.read_txt_file_densityprofile(filename_densityprofile)

        z_read   = layers[ilayer].material.densityComposition.densityprofile[:,1]
        rho_read = layers[ilayer].material.densityComposition.densityprofile[:,0]

        if ( depth > z_read[len(z_read)-1] ):
            density = rho_read[len(z_read)-1]
        else:
            density = pos2value_1d(z_read,rho_read,depth)

        return density

    def read_txt_file_densityprofile(self,filename_densityprofile):

        if (os.path.isfile(filename_densityprofile)):
            data = np.loadtxt(filename_densityprofile,dtype='float',delimiter=' ')
        else:
            print('Problem with unread file '+filename_densityprofile)
            quit()
            
        return data

class melt_fraction_obj(object):
       
    def get_value(self,ilayer,depth,pressure,temp,layers,path=None):
        """Dispatch method"""
        density_type=layers[ilayer].material.densityComposition.use_tidv
        method_name = 'density_type' + str(density_type)
        # Get the method from 'self'. Default to a lambda.
        method = getattr(self, method_name, 0)
        return method(ilayer,depth,pressure,temp,layers,path=path)
    
    def density_type0(self,ilayer,depth,pressure,temp,layers,path):
        return 0.0e0

    def density_type1(self,ilayer,depth,pressure,temp,layers,path):
        return 0.0e0
    
    def density_type2(self,ilayer,depth,pressure,temp,layers,path):
        return 0.0e0

    def density_type7(self,ilayer,depth,pressure,temp,layers,path):
        # Density read in a table
        # defined by layer properties
        if layers[ilayer].material.densityComposition.perplex_name is None:
            print('You use use_tidv = 7 for this layer {} but perplex_name is None.'.format(ilayer)) 
            raise TypeError
        else:
            perplex_name = layers[ilayer].material.densityComposition.perplex_name
        
        if ( not layers[ilayer].material.densityComposition.perplex ):
            layers[ilayer].material.densityComposition.perplex = read_binary_file_perplex(perplex_name,path)
        
        pmax = np.max(layers[ilayer].material.densityComposition.perplex.P)
        tmax = np.max(layers[ilayer].material.densityComposition.perplex.T)
                
        if ( not layers[ilayer].f_interp_melt_fraction):
            
            perplex_np = layers[ilayer].material.densityComposition.perplex.np
            perplex_nt = layers[ilayer].material.densityComposition.perplex.nt
            value = layers[ilayer].material.densityComposition.perplex.melt.reshape(perplex_np, perplex_nt)
            T = np.linspace(np.min(layers[ilayer].material.densityComposition.perplex.T),
                            tmax,
                            perplex_nt)
            P = np.linspace(np.min(layers[ilayer].material.densityComposition.perplex.P),
                            pmax,
                            perplex_np)            
            layers[ilayer].f_interp_melt_fraction = RegularGridInterpolator((np.array(P)/pmax, np.array(T)/tmax), value)
        
        if ((pressure>pmax) | (temp>tmax)):
            print('WARNING: T {} ; Tmax {} P {} ; Pmax {}'.format(temp,tmax,pressure,pmax))
            if (pressure>pmax):
                pressure = pmax
            if (temp>tmax):
                temp = tmax
                    
        melt_fraction = layers[ilayer].f_interp_melt_fraction([pressure/pmax,temp/tmax])[0]

    def density_type8(self,ilayer,depth,pressure,temp,layers):
        return 0.0e0

        # Melt fraction should be reduced to fit oceanic crust thickness - indeed no fractionation on residual composition
        
        return melt_fraction


 
class thermal_expansion(object):
    
    def get_value(self,ilayer,depth,pressure,temp,layers,alpha_type):
        """Dispatch method"""
        method_name = 'alpha_type' + str(alpha_type)
        # Get the method from 'self'. Default to a lambda.
        method = getattr(self, method_name, 0)
        return method(ilayer,depth,pressure,temp,layers)
    
    def alpha_type0(self,ilayer,depth,pressure,temp,layers):
        # Density constant
        # defined by layer properties
        alpha = 0.0e0
        print(alpha)
        return alpha
    
    def alpha_type1(self,ilayer,depth,pressure,temp,layers):
        # Density constant
        # defined by layer properties
        alpha = layers[ilayer].material.thermalProperties.alpha
        return alpha

    def alpha_type2(self,ilayer,depth,pressure,temp,layers):
        # Density constant
        # defined by layer properties
        alpha1 = layers[ilayer].material.thermalProperties.alpha1
        alpha2 = layers[ilayer].material.thermalProperties.alpha2
        alphaT1 = layers[ilayer].material.thermalProperties.alphaT1
        alphaT2 = layers[ilayer].material.thermalProperties.alphaT2
        return alpha1, alpha2, alphaT1, alphaT2

class compressibility(object):
    pass

class density_table(object):
    pass


