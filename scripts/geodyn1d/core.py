# ***********************************************************************************
# * Copyright © 2018-2100, <copyright holders: Thomas Theunissen>                   *
# *                                                                                 *
# * Permission is hereby granted, free of charge, to any person obtaining a copy of *
# * this software and associated documentation files (the “Software”), to deal in   *
# * the Software without restriction, including without limitation the rights to    *
# * use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies   *
# * of the Software, and to permit persons to whom the Software is furnished to do  *
# * so, subject to the following conditions:                                        *
# *                                                                                 *
# * The above copyright notice and this permission notice shall be included in all  *
# * copies or substantial portions of the Software.                                 *
# *                                                                                 *
# * The Software is provided “as is”, without warranty of any kind, express or      *
# * implied, including but not limited to the warranties of merchantability,        *
# * fitness for a particular purpose and noninfringement. In no event shall the     *
# * authors or copyright holders Thomas Theunissen be liable for any claim, damages *
# * or other liability, whether in an action of contract, tort or otherwise,        *
# * arising from, out of or in connection with the software or the use or other     *
# * dealings in the Software.                                                       *
# * Except as contained in this notice, the name of the <copyright holders: Thomas  *
# * Theunissen> shall not be used in advertising or otherwise to promote the sale,  *
# * use or other dealings in this Software without prior written authorization from *
# * the <copyright holders: Thomas Theunissen>.                                     *
# ***********************************************************************************

# ---------------------------------------------------------------
# -------Definition of global variables--------------------------
# ---------------------------------------------------------------

import sys
import copy
from scipy.interpolate import interp1d
from .layer import *
from .material import *
from .thermalBc import *

from .geotherm import *
from .isostasy import *
from .sedimento import *
from .mechanics import *
from .visualization import *

rho_water     = 1027.0e0 # ocean water density Kg/m3
rho_sediments = 2400.0e0 # average sediement density Kg/m3 
g             = 9.80665  # gravitational acceleration (m/s2)


# ---------------------------------------------------------------
# -------Load database-------------------------------------------
# ---------------------------------------------------------------

# specific properties for standard lithospheres
from .library import *

# ---------------------------------------------------------------
# -------Classes-------------------------------------------------
# ---------------------------------------------------------------

class PKIsFalseException(Exception):
    pass

def check_numlayers(n1,n2,n3):
        if ((n1 == n2) and (n2 == n3)):
            check_numlayers = True
        else:
            check_numlayers = False
        return check_numlayers

class lithosphere (object):
    
    # this is incompatible with BoundInnerClass
    #__slots__ = ['numlayers','nature_layers', 'thicknesses','thermalBc','layers','lith']

    def __init__(self,lith=None,
                      numlayers=7,
                      nature_layers=['matUC','matMC','matLC','matLM1','matLM2','matLM3','matSLM'],
                      thicknesses=[15e3,10e3,10e3,55e3,90e3,0.0e0,420e3],
                      thermalBc=['thermBcUC','thermBcMC','thermBcLC','thermBcLM1','thermBcLM2','thermBcLM3','thermBcSLM'],
                      printOut = False,
                      layers = None):
        try:
            self.lith          = lith
            if (lith):
                self.numlayers      = eval(lith)["numlayers"]
                self.nature_layers  = eval(lith)["nature_layers"]
                self.thicknesses    = eval(lith)["thicknesses"]
                self.thermalBc      = eval(lith)["thermalBc"]
            else:
                self.numlayers     = numlayers
                self.nature_layers = nature_layers
                self.thicknesses   = thicknesses
                self.thermalBc     = thermalBc
                self.lith          = 'lith_default'
                lith_default       = {  "numlayers":     self.numlayers,
                                        "nature_layers": self.nature_layers,
                                        "thicknesses":   self.thicknesses,
                                        "thermalBc":     self.thermalBc
                                     }
            correct_numlayers  = check_numlayers(self.numlayers,len(self.nature_layers),len(self.thicknesses))
            if not correct_numlayers: raise PKIsFalseException()            
        except (PKIsFalseException):
            print('The size of \'nature_layers\' and \'thicknesses\' should be identical to \'numlayers\'')
            sys.exit()

        self.layers             = []
        if (printOut):
            print('    num nature thickness')
        for ilayer in range(0,self.numlayers):
            newlayer = 'layer'+str(ilayer)
            if (printOut):
                print('{0:>7} {1:>6} {2:>9}'.format(newlayer,self.nature_layers[ilayer],self.thicknesses[ilayer]))
            # CERTAINLY USELESS !!! dynamically generate attribute for each layer
            # setattr(self, newlayer,layer(self.nature_layers[ilayer],self.thicknesses[ilayer],self.thermalBc[ilayer]))
            specific_book_properties = eval(self.lith)
            self.layers.append( layer( layerid       = ilayer,
                                       name_material = self.nature_layers[ilayer],
                                       thickness     = self.thicknesses[ilayer],
                                       thermalBc     = self.thermalBc[ilayer],
                                       specific_book_properties = specific_book_properties)
                              )
        if (printOut):
            print('\n')
        
        # we define those key variable at None
        self.geotherm           = None
        self.P_rho_DepthProfile = None
        self.StrengthProfile    = None
        self.sea_level          = None
        self.water_depth        = None
        self.sediment_thickness = None
    
    def update_materials(self,specific_book_properties={}):
        if (len(specific_book_properties)>0):
            #print(eval(self.lith))
            #new_book = merge_two_dicts(eval(self.lith),specific_book_properties)
            #new_book = eval(self.lith)

            # we initiatate book keys with lithosphere properties
            new_book = {}
            for item in self.nature_layers:
                new_book[item] = {}
            #TODO should be possible to update it
            new_book['thicknesses'] = self.thicknesses
            new_book['numlayers']   = self.numlayers
            #---------------------------------
            for item in self.thermalBc:
                new_book[item] = {}


            #for item in eval(self.lith).items():
            keynotupdated = []
            for item in new_book.items():
                key = item[0]
#                print("check key {}".format(key))
                # check if update necessary for this key
                key_is_found = False
                for item_new in specific_book_properties.items():
                    key_new = item_new[0]
                    if (key_new==key):
                        key_is_found = True
                        if (isinstance(item[1], dict) ):
                            new_subbook = merge_two_dicts(item[1],item_new[1])
                            new_book[key] = new_subbook
                            break
                        else:
                            new_book[key] = item_new[1]
                if not key_is_found:
                    for item_old in eval(self.lith).items():
                        key_old = item_old[0]
                        if (key_old==key):
                            key_is_found = True
 #                           print("we keep old content for this key {}".format(key))
                            if (isinstance(item[1], dict) ):
                                #print(item_old[1])
                                new_subbook = merge_two_dicts(item[1],item_old[1])
                                new_book[key] = new_subbook
                                break
                            else:
                                #print(item_old[1])
                                new_book[key] = item_old[1]
            
            #print(' update_materials --------------')
            #print(new_book)
            #print('--------------------------------')
            #TODO lithosphere dictionaries should be updated not only the layers
            # ===> update cannot be run two times consecutively !!

            self.layers             = []
            for ilayer in range(0,self.numlayers):
                self.layers.append(layer( layerid = ilayer,
                                    name_material = self.nature_layers[ilayer],
                                    thickness     = self.thicknesses[ilayer],
                                    thermalBc     = self.thermalBc[ilayer],
                                    specific_book_properties = new_book)
                                  )
        return
    
    def get_steady_state_geotherm(self,printOut=True,dz=200):
        self.geotherm = compute_steady_state_geotherm(layers=self.layers,printOut=printOut,dz=dz)
        return

    def get_imposed_geotherm(self,filename,printOut=False):
        self.geotherm = read_imposed_geotherm(layers=self.layers,filename=filename,printOut=printOut)
        return
    
    def get_pressure_density_profile(self,printOut=False, filename=None, drho=0, zdrho=[6500,125000],path=None):
        if (self.geotherm):
            pass
        else:
            self.get_steady_state_geotherm(printOut=printOut)
        self.P_rho_DepthProfile = compute_pressure_density_profile(layers=self.layers,geotherm=self.geotherm,printOut=printOut,filename=filename, drho=drho, zdrho=zdrho,path=path)
        return
   
    def get_strength_profile(self,strain_rate=1e-16,printOut=False):
        if (self.geotherm):
            pass
        else:
            self.get_steady_state_geotherm(printOut=printOut)
        if (self.P_rho_DepthProfile):
            pass
        else:
            self.P_rho_DepthProfile = compute_pressure_density_profile(layers=self.layers,geotherm=self.geotherm,printOut=printOut)

        self.StrengthProfile    = compute_strength_profile(layers=self.layers,geotherm=self.geotherm,P_rho_DepthProfile=self.P_rho_DepthProfile,strain_rate=strain_rate,printOut=printOut)

        return

    def get_reference_pressure(self,href,printOut=False,filename=None,path=None):
        Pref=0.0e0
        if (self.P_rho_DepthProfile):
            pass
        else:
            if (self.geotherm):
                pass
            else:
                self.get_steady_state_geotherm(printOut=printOut)
            self.get_pressure_density_profile(printOut=printOut,filename=filename,path=path)
                
        # a spline interpolation 2nd order
        f = interp1d(self.P_rho_DepthProfile.Z, self.P_rho_DepthProfile.P, kind='quadratic')

        Pref = f(href).tolist()
        
        return Pref

    def set_sea_level(self,sea_level,printOut = False):
        if (sea_level>0):
            if (printOut):
                print('Sea level is defined below top of reference lithosphere (the one used to compute "pref").\nIt must be a negative value in meters.\nSea level sets to zero.')
            sea_level = 0.0e0
        self.sea_level = sea_level
        
    def set_water_depth(self,water_depth,printOut = False):
        if (water_depth<=0):
            if (printOut):
                print('Water depth is defined on top of the current lithosphere.\nIt must be a positive value in meters.\nWater depth sets to zero.')
            water_depth = 1.0e-9
        self.water_depth = water_depth
        
    def set_sediment_thickness(self,sediment_thickness,printOut = False):
        if (sediment_thickness<=0):
            if (printOut):
                print('Sediment thickness is defined on top of the current lithosphere.\nIt must be a positive value in meters.\nSediment thickness sets to zero.')
            sediment_thickness = 1.0e-9
        self.sediment_thickness = sediment_thickness

    def get_subsidence(self,pref,href,SeaLevelFixed = False,printOut=False,filename=None,path=None):
        """ The module computes the subsidense for a fixed sea level or a fixed water depth.
        """
        h = 0.0e0
        if (self.P_rho_DepthProfile):
            pass
        else:
            if (self.geotherm):
                pass
            else:
                self.get_steady_state_geotherm(printOut=printOut)
            self.get_pressure_density_profile(printOut=printOut,filename=filename,path=path)
        
        if (SeaLevelFixed == True):
            if ( not self.sea_level):
                print('No results. You must define a sea_level before.\nRun "yourLithosphere.set_sea_level(-300)"\n')
                return 0,0,0
            else:
                if (printOut):
                    print('SEA LEVEL FIXED {0:7.1f} m below top of reference lithosphere (the one used to compute "pref")'.format(self.sea_level))
            
            # compute subsidence and water thickness based on sea_level
            sea_level         = self.sea_level
            Pwater            = 0.0e0
            delta_water_depth = 9.9e9
            water_depth_prev  = 0.0e0
            water_depth       = 0.0e0
            # a spline interpolation 2nd order on pressure profile
            f = interp1d(self.P_rho_DepthProfile.P, self.P_rho_DepthProfile.Z, kind='quadratic')
            # look for the effect of water depth
            while (delta_water_depth>5.0e-3):
                h  = f(pref-Pwater).tolist()
                dh = h - href
                if (dh<sea_level):
                    water_depth = sea_level - dh
                    Pwater      = (rho_water*g*water_depth)
                else:
                    water_depth = 0.0e0
                    Pwater      = 0.0e0
                delta_water_depth = abs(water_depth_prev - water_depth)
                water_depth_prev  = water_depth
            self.water_depth = water_depth
        else:
            if ( not self.water_depth):
                print('No results. You must define a water_depth before.\nRun "yourLithosphere.set_water_depth(3000)"\n')
                return 0,0,0
            else:
                if (printOut):
                    print('WATER DEPTH FIXED {0:7.1f} m on top of the current lithosphere'.format(self.water_depth))
            
            # compute subsidence and water thickness based on sea_level
            Pwater            = (rho_water*g*self.water_depth)
            delta_water_depth = 9.9e9
            water_depth_prev  = self.water_depth
            water_depth       = self.water_depth
            # a spline interpolation 2nd order on pressure profile
            f = interp1d(self.P_rho_DepthProfile.P, self.P_rho_DepthProfile.Z, kind='quadratic')
            # look for the effect of water depth
            while (delta_water_depth>5.0e-3):
                h  = f(pref-Pwater).tolist()
                dh = h - href
                if (dh>-1.0e0*water_depth):
                    if (dh>0):
                        water_depth = 0.0e0
                    else:
                        water_depth = -1.0e0*dh
                    Pwater      = (rho_water*g*water_depth)
                delta_water_depth = abs(water_depth_prev - water_depth)
                water_depth_prev  = water_depth
            if (water_depth != self.water_depth):
                if (printOut):
                    print('Water depth has been modified NEW WATER DEPTH = {0:7.1f} m'.format(water_depth))
                self.water_depth = water_depth
            if (dh>0):
                self.sea_level   = 0.0e0
            else:
                self.sea_level   = dh + self.water_depth
            
        return dh,self.water_depth,self.sea_level

    def get_subsidence_with_sediments(self,pref,href,density_sediments=rho_sediments,SeaLevelFixed = False,printOut=False,filename=None,path=None):
        """ The module computes the subsidense for a fixed sea level or for a fixed water depth.
            The sediment thickness is defined as follow:
                water_depth = sea_level - dh - sediment_thickness (with dh<0 and sea_level<0)
            TODO: consider using the module "sedimento" to consider sediment compaction instead of constant average density
        """
        h = 0.0e0
        if (self.P_rho_DepthProfile):
            pass
        else:
            if (self.geotherm):
                pass
            else:
                self.get_steady_state_geotherm(printOut=printOut)
            self.get_pressure_density_profile(printOut=printOut,filename=filename,path=path)
      
        if ( not self.sediment_thickness):
                print('No results. You must define a sediment_thickness before.\nRun "yourLithosphere.set_sediment_thickness(3000)"\n')
                return 0,0,0,0
        else:
            if (printOut):
                print('SEDIMENT THICKNESS FIXED {0:7.1f} m on top of the current lithosphere'.format(self.sediment_thickness))

        # sediment loading
        sediment_thickness = self.sediment_thickness
        psed = (density_sediments*g*sediment_thickness)
        if (printOut):
            print('Initial sed thickness {} psed {}'.format(sediment_thickness,psed))
        
        if (SeaLevelFixed == True):
            if ( not self.sea_level):
                print('No results. You must define a sea_level before.\nRun "yourLithosphere.set_sea_level(-300)"\n')
                return 0,0,0,0
            else:
                if (printOut):
                    print('SEA LEVEL FIXED {0:7.1f} m below top of reference lithosphere (the one used to compute "pref")'.format(self.sea_level))
            
            # compute subsidence and water thickness based on sea_level
            sea_level         = self.sea_level
            Pwater            = 0.0e0
            delta_water_depth = 9.9e9
            water_depth_prev  = 0.0e0
            water_depth       = 0.0e0
            delta_sediment_thickness = 9.9e9
            sediment_thickness_prev  = sediment_thickness
            # a spline interpolation 2nd order on pressure profile
            f = interp1d(self.P_rho_DepthProfile.P, self.P_rho_DepthProfile.Z, kind='quadratic')
            # look for the effect of water depth
            while ((delta_water_depth>5.0e-3)|(delta_sediment_thickness>5.0e-3)):
                if (printOut):
                    print('LOOP --------------------------------- water {} sed. {}'.format(water_depth,sediment_thickness))
                h  = f(pref-Pwater-psed).tolist()
                dh = h - href
                if (dh<sea_level):
                    water_depth = sea_level - dh - sediment_thickness
                    if (water_depth<0):
                        if (printOut):
                            print('too much sediment above sea level !!!')
                        # we allow sediment to be above sea_level
                        # sediment_thickness = sea_level - dh
                        water_depth        = 0.0e0
                        if (abs(dh)<sediment_thickness):
                            if (printOut):
                                print('too much sediment above href, we change sed. thick. !!!')
                            sediment_thickness = abs(dh)
                        else:
                            if (sediment_thickness<self.sediment_thickness):
                                if (printOut):
                                    print('check A !!!')
                                if (abs(dh)-sediment_thickness>self.sediment_thickness-sediment_thickness):
                                    sediment_thickness = self.sediment_thickness
                                    if (printOut):
                                        print('check B !!!')
                                else:
                                    sediment_thickness = abs(dh)
                                    if (printOut):
                                        print('check C !!!')

                else:
                    if (printOut):
                        print('not enough subsidence')
                    water_depth = 0.0e0
                    Pwater      = 0.0e0
                    if (abs(dh)<sediment_thickness):
                        if (printOut):
                            print('too much sediment above href, we change sed. thick. !!!')
                        sediment_thickness = abs(dh)
                    else:
                        if (sediment_thickness<self.sediment_thickness):
                            if (printOut):
                                print('check A !!!')
                            if (abs(dh)-sediment_thickness>self.sediment_thickness-sediment_thickness):
                                sediment_thickness = self.sediment_thickness
                                if (printOut):
                                    print('check B !!!')
                            else:
                                sediment_thickness = abs(dh)
                                if (printOut):
                                    print('check C !!!')
                                
                Pwater                   = (rho_water*g*water_depth)
                psed                     = (density_sediments*g*sediment_thickness)
                delta_water_depth        = abs(water_depth_prev - water_depth)
                water_depth_prev         = water_depth
                delta_sediment_thickness = abs(sediment_thickness_prev - sediment_thickness)
                sediment_thickness_prev  = sediment_thickness
                
            self.water_depth        = water_depth
            if (sediment_thickness != self.sediment_thickness):
                if (printOut):
                    print('Sediment thickeness has been modified NEW THICKNESS = {0:7.1f} m'.format(sediment_thickness))
                self.sediment_thickness = sediment_thickness
            
        else:
            if ( not self.water_depth):
                print('No results. You must define a water_depth before.\nRun "yourLithosphere.set_water_depth(3000)"\n')
                return 0,0,0,0
            else:
                if (printOut):
                    print('WATER DEPTH FIXED {0:7.1f} m on top of the current lithosphere'.format(self.water_depth))
                        
            # compute subsidence and water thickness based on sea_level
            Pwater            = (rho_water*g*self.water_depth)
            delta_water_depth = 9.9e9
            water_depth_prev  = self.water_depth
            water_depth       = self.water_depth
            delta_sediment_thickness = 9.9e9
            sediment_thickness_prev  = sediment_thickness
            # a spline interpolation 2nd order on pressure profile
            f = interp1d(self.P_rho_DepthProfile.P, self.P_rho_DepthProfile.Z, kind='quadratic')
            # look for the effect of water depth
            while ((delta_water_depth>5.0e-3)|(delta_sediment_thickness>5.0e-3)):
                if (pref-Pwater-psed>self.P_rho_DepthProfile.P.max()):
                    print("Change ref_depth to lower value (not bottom depth)")
                    quit()
                h  = f(pref-Pwater-psed).tolist()
                dh = h - href
                # water_depth = sea_level - dh - sediment_thickness
                # dh = -water_depth (sea = 0 sed = 0)
                if (dh>(-1.0e0*water_depth-sediment_thickness)):
                    if (dh>0):
                        water_depth = 0.0e0
                        if (abs(dh)<sediment_thickness):
                            if (printOut):
                                print('too much sediment above href, we change sed. thick. !!!')
                                sediment_thickness = abs(dh)
                            else:
                                if (sediment_thickness<self.sediment_thickness):
                                    if (printOut):
                                        print('check A !!!')
                                    if (abs(dh)-sediment_thickness>self.sediment_thickness-sediment_thickness):
                                        sediment_thickness = self.sediment_thickness
                                        if (printOut):
                                            print('check B !!!')
                                    else:
                                        sediment_thickness = abs(dh)
                                        if (printOut):
                                            print('check C !!!')
                    else:
                        water_depth = -1.0e0*dh-sediment_thickness
                        if (water_depth<0):
                            if (printOut):
                                print('too much sediment above sea level !!!')
                            # we allow sediment to be above sea_level
                            # sediment_thickness = sea_level - dh
                            water_depth        = 0.0e0
                            if (abs(dh)<sediment_thickness):
                                if (printOut):
                                    print('too much sediment above href, we change sed. thick. !!!')
                                sediment_thickness = abs(dh)
                            else:
                                if (sediment_thickness<self.sediment_thickness):
                                    if (printOut):
                                        print('check A !!!')
                                    if (abs(dh)-sediment_thickness>self.sediment_thickness-sediment_thickness):
                                        sediment_thickness = self.sediment_thickness
                                        if (printOut):
                                            print('check B !!!')
                                    else:
                                        sediment_thickness = abs(dh)
                                        if (printOut):
                                            print('check C !!!')
                        
                Pwater                   = (rho_water*g*water_depth)
                psed                     = (density_sediments*g*sediment_thickness)
                delta_water_depth = abs(water_depth_prev - water_depth)
                water_depth_prev  = water_depth
                delta_sediment_thickness = abs(sediment_thickness_prev - sediment_thickness)
                sediment_thickness_prev  = sediment_thickness
            
            if (sediment_thickness != self.sediment_thickness):
                if (printOut):
                    print('Sediment thickeness has been modified NEW THICKNESS = {0:7.1f} m'.format(sediment_thickness))
                self.sediment_thickness = sediment_thickness
            if (water_depth != self.water_depth):
                if (printOut):
                    print('Water depth has been modified NEW WATER DEPTH = {0:7.1f} m'.format(water_depth))
                self.water_depth = water_depth
            if (dh>0):
                self.sea_level   = 0.0e0
            else:
                # water_depth = sea_level - dh - sediment_thickness
                # dh = -water_depth (sea = 0 sed = 0)
                self.sea_level   = dh + self.water_depth + sediment_thickness
            
        return dh,self.water_depth,self.sea_level,self.sediment_thickness

