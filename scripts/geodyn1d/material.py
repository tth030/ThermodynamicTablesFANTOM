
from .flowlaws import *
from .plasticities import *
import copy

class material_object(object):
    
    def __init__(self,color = None,
                      name  = None,
                      mechanicalProperties = None,
                      thermalProperties    = None,
                      densityComposition   = None):
        self.color = color
        self.name  = name
        self.mechanicalProperties = mechanicalProperties
        self.thermalProperties    = thermalProperties
        self.densityComposition   = densityComposition

class mechanicalProperties(object):
    
    def __init__(self,
                      flowlaw    = None,
                      plasticity = None):
        self.flowlaw    = flowlaw
        self.plasticity = plasticity

class thermalProperties(object):
    
    def __init__(self,
                      alpha    = None, alpha1  = None, alpha2    = None, alphaT1   = None, alphaT2 = None,
                      alphaP1  = None, alphaP2 = None, alphaP1_f = None, alphaP2_f = None, # Thermal expansion K-1 =f(T,P)
                      beta  = None,  # Compressibility Pa-1  =f(T,P)
                      k     = None,  # Thermal conductivity W/m/K
                      k_type = None, # 0: constant thermal conductivities 1: temperature dependant conductivities
                      H     = None,  # Heat production W/m3
                      Cp    = None,  # Specific Heat Capacity J/K/Kg
                      D     = None,  # Thermal diffusivity m2/s
                      T0    = None): # temperature for lab measurement (=surface temp.))
        self.alpha    = alpha
        self.alpha1   = alpha1
        self.alpha2   = alpha2
        self.alphaT1   = alphaT1
        self.alphaT2   = alphaT2
        self.alphaP1   = alphaP1
        self.alphaP2   = alphaP2
        self.alphaP1_f = alphaP1_f
        self.alphaP2_f = alphaP2_f
        self.k        = k
        self.k_type   = k_type
        self.H        = H
        self.Cp       = Cp
        self.D        = D
        self.T0       = T0

#use_tidv triggers the use of thermally induced density variations
#use_tidv = 0 means no thermal expansion
#use_tidv = 1 means thermal expansion with only using coefficient thermal_expansion_i
#use_tidv = 2 means temperature dependent thermal expansion linearly scaled from thermal_expansion_1_i to thermal_expansion_2_i
#             within the temperature range thermal_expansion_t1_i to thermal_expansion_t2_i
#         = 3 temperature and pressure dependent thermal expansion
#         = 4 NOT IMPLEMENTED constant thermal expansion and compressibility
#         = 5 NOT IMPLEMENTED temp dep. thermal expansion and compressibility
#         = 6 NOT IMPLEMENTED temp. and P. dep. thermal expansion and compressibility
#         = 7 READ FROM TABLE
#         = 8 DELTARHO that apply on ref. density READ FROM TABLE and constant thermal expansion similar to case "1"
#         = 9 DELTARHO that apply on ref. density READ FROM TABLE and thermal expansion(T) similar to case "2"
#         = 10 DELTARHO that apply on ref. density READ FROM TABLE and thermal expansion(T,P) similar to case "3"

class densityComposition(object):
    
    def __init__(self,use_tidv      = None, # use_tidv triggers the use of thermally induced density variations
                      perplex_name  = None, # table name with density,alpha,beta,melt,deltarho
                      perplex       = None, # perplex table
                      T0    = None,         # temperature for lab measurement (=surface temp.))
                      deltarho = 0.0e0,     # density variation applied on material after P,T dependent density estimation (use_tidv)
                      rho   = None,         # density Kg/m3 =f(T,P,X)
                      beta  = None,         # Compressibility Pa-1  =f(T,P)
                      filename_densityprofile = None, # Density profile filename
                      densityprofile = None ):        # Density profile 
        self.rho  = rho
        self.beta = beta
        self.T0   = T0
        self.use_tidv = use_tidv
        self.perplex_name = perplex_name
        self.perplex  = perplex
        self.filename_densityprofile = filename_densityprofile
        self.densityprofile = densityprofile
        
def get_book_material(material,book=None):
    #list_attr = dir(material) ; list_attr = [s for s in list_attr if "__" not in s]   
    if ( not book):
        book = {}
    attrs = vars(material)
    for item in attrs.items():      
        if ( item[1].__class__.__name__ in globals()):
            if (isinstance(item[1],eval(item[1].__class__.__name__))):
                get_book_material(item[1],book)
            else:
                book[item[0]] = item[1]
        else:
            book[item[0]] = item[1]
    return book

def update_material(material,material_name,book_new_parameters,material_name_init=None,fullpath=None,printOut=False):
    #book_material = get_book_material(material)
    #material_name = 'matUC'
    if ( not fullpath ):
        fullpath=material_name
    if ( not material_name_init ):
        material_name_init=material_name
    attrs = vars(material)
    for item in attrs.items():
        if ( item[1].__class__.__name__ in globals()):
            if ( not material_name_init in fullpath):
                fullpath = material_name_init+'.'+fullpath
            if (isinstance(item[1],eval(item[1].__class__.__name__))):
                fullpath = fullpath+'.'+item[1].__class__.__name__
                update_material(item[1],item[0],book_new_parameters,fullpath=fullpath,material_name_init=material_name_init,printOut=printOut)
                fullpath=material_name
            else:
                for key,value in book_new_parameters.items():
                    if ( key == item[0]):
                        if ( value != item[1]):
                            if (printOut):
                                print('update {}: {} becomes {}'.format(key,item[1],value))
                            fullpathname = fullpath+'.'+key
                            setattr(eval(fullpath),key,value)
                        break
        else:
            for key,value in book_new_parameters.items():
                if ( key == item[0]):
                    if ( value != item[1]):
                        if (printOut):
                            print('update {}: {} becomes {}'.format(key,item[1],value))
                        fullpathname = fullpath+'.'+key
                        setattr(eval(fullpath),key,value)
                    break
        
    return
    
    

matUC         = material_object()
matUC.name    = "Material Upper-Crust"
matUC.color   = [1.0e0, 0.5490e0, 0.0]
matUC.mechanicalProperties = mechanicalProperties()
matUC.mechanicalProperties.flowlaw    = copy.deepcopy(WetQz)
matUC.mechanicalProperties.plasticity = copy.deepcopy(mohrCoulomb)
matUC.densityComposition = densityComposition()
matUC.densityComposition.beta = 8.0e-11
matUC.densityComposition.rho      = 2750
matUC.densityComposition.deltarho = 0.0e0
matUC.densityComposition.T0   = 0.0e0
matUC.thermalProperties = thermalProperties()
matUC.thermalProperties.alpha = 3.1e-5 
matUC.thermalProperties.k     = 2.25
matUC.thermalProperties.k_type = 0
matUC.thermalProperties.H     = 1.78e-6 #lith120 1.196e-6
matUC.thermalProperties.D     = 1.0e-6
matUC.thermalProperties.Cp    = matUC.thermalProperties.k / (matUC.densityComposition.rho*matUC.thermalProperties.D)
matUC.thermalProperties.T0    = 0.0e0
matUC.thermalProperties.alpha1    = 3.1e-5
matUC.thermalProperties.alpha2    = 3.1e-5 
matUC.thermalProperties.alphaT1   = 0.0e0
matUC.thermalProperties.alphaT2   = 0.0e0
matUC.thermalProperties.alphaP1   = 1e5
matUC.thermalProperties.alphaP2   = 1e5
matUC.thermalProperties.alphaP1_f = 1.0e0
matUC.thermalProperties.alphaP2_f = 1.0e0
matUC.densityComposition.use_tidv  = 1
matUC.densityComposition.perplex_name = None
matUC.densityComposition.perplex      = None


        
matMC         = material_object()
matMC.name    = "Material Middle-Crust"
matMC.color   = [0.9608e0, 0.9608e0, 0.9608e0]
matMC.mechanicalProperties = mechanicalProperties()
matMC.mechanicalProperties.flowlaw = copy.deepcopy(WetQz)
matMC.mechanicalProperties.plasticity = copy.deepcopy(mohrCoulomb)
matMC.densityComposition = densityComposition()
matMC.densityComposition.beta = 8.0e-11
matMC.densityComposition.rho      = 2750
matMC.densityComposition.deltarho = 0.0e0
matMC.densityComposition.T0   = 0.0e0
matMC.thermalProperties = thermalProperties()
matMC.thermalProperties.alpha = 3.1e-5 
matMC.thermalProperties.k     = 2.25
matMC.thermalProperties.k_type = 0
matMC.thermalProperties.H     = 1.78e-6 #lith120 1.196e-6
matMC.thermalProperties.D     = 1.0e-6
matMC.thermalProperties.Cp    = matMC.thermalProperties.k / (matMC.densityComposition.rho*matMC.thermalProperties.D)
matMC.thermalProperties.T0    = 0.0e0
matMC.thermalProperties.alpha1    = 3.1e-5 
matMC.thermalProperties.alpha2    = 3.1e-5 
matMC.thermalProperties.alphaT1   = 0.0e0
matMC.thermalProperties.alphaT2   = 0.0e0
matMC.thermalProperties.alphaP1   = 1e5
matMC.thermalProperties.alphaP2   = 1e5
matMC.thermalProperties.alphaP1_f = 1.0e0
matMC.thermalProperties.alphaP2_f = 1.0e0
matMC.densityComposition.use_tidv  = 1
matMC.densityComposition.perplex_name = None
matMC.densityComposition.perplex      = None


matLC         = material_object()
matLC.name    = "Lower-Crust"
matLC.color   = [0.8824e0, 0.8824e0, 0.8824e0]
matLC.mechanicalProperties = mechanicalProperties()
matLC.mechanicalProperties.flowlaw = copy.deepcopy(WetQz)
matLC.mechanicalProperties.plasticity = copy.deepcopy(mohrCoulomb)
matLC.densityComposition = densityComposition()
matLC.densityComposition.beta    = 8.0e-11
matLC.densityComposition.rho      = 2900
matLC.densityComposition.deltarho = 0.0e0
matLC.densityComposition.T0      = 0.0e0
matLC.thermalProperties = thermalProperties()
matLC.thermalProperties.alpha   = 3.1e-5 
matLC.thermalProperties.k       = 2.25
matLC.thermalProperties.k_type = 0
matLC.thermalProperties.H       = 0.79e-6 #lith120 0.483e-6
matLC.thermalProperties.D       = 1.0e-6
matLC.thermalProperties.Cp      = matLC.thermalProperties.k / (matLC.densityComposition.rho*matLC.thermalProperties.D)
matLC.thermalProperties.T0      = 0.0e0
matLC.thermalProperties.alpha1    = 3.1e-5 
matLC.thermalProperties.alpha2    = 3.1e-5 
matLC.thermalProperties.alphaT1   = 0.0e0
matLC.thermalProperties.alphaT2   = 0.0e0
matLC.thermalProperties.alphaP1   = 1e5
matLC.thermalProperties.alphaP2   = 1e5
matLC.thermalProperties.alphaP1_f = 1.0e0
matLC.thermalProperties.alphaP2_f = 1.0e0
matLC.densityComposition.use_tidv  = 1
matLC.densityComposition.perplex_name = None
matLC.densityComposition.perplex      = None


matLM1            = material_object()
matLM1.name       = "Lithospheric mantle 1"
matLM1.color      = [0.0, 0.7451e0, 0.5882e0]
matLM1.mechanicalProperties = mechanicalProperties()
matLM1.mechanicalProperties.flowlaw    = copy.deepcopy(WetOl)
matLM1.mechanicalProperties.flowlaw.Fc = 5.0
matLM1.mechanicalProperties.plasticity = copy.deepcopy(mohrCoulomb_mantle)
matLM1.densityComposition = densityComposition()
matLM1.densityComposition.beta       = 8.0e-11
matLM1.densityComposition.rho        = 3300
matLM1.densityComposition.deltarho   = 0.0e0
matLM1.densityComposition.T0         = 0.0e0
matLM1.thermalProperties = thermalProperties()
matLM1.thermalProperties.alpha      = 3.1e-5 
matLM1.thermalProperties.k          = 2.25
matLM1.thermalProperties.k_type = 0
matLM1.thermalProperties.H          = 0.0
matLM1.thermalProperties.D          = 1.0e-6
matLM1.thermalProperties.Cp         = matLM1.thermalProperties.k / (matLM1.densityComposition.rho*matLM1.thermalProperties.D)
matLM1.thermalProperties.T0         = 0.0e0
matLM1.thermalProperties.alpha1    = 3.0e-5 
matLM1.thermalProperties.alpha2    = 3.5e-5 
matLM1.thermalProperties.alphaT1   = 500.0e0
matLM1.thermalProperties.alphaT2   = 2000.0
matLM1.thermalProperties.alphaP1   = 1e5
matLM1.thermalProperties.alphaP2   = 19.8e9
matLM1.thermalProperties.alphaP1_f = 1.0e0
matLM1.thermalProperties.alphaP2_f = 0.5e0
matLM1.densityComposition.use_tidv  = 1
matLM1.densityComposition.perplex_name = None
matLM1.densityComposition.perplex      = None


matLM2            = material_object()
matLM2.name       = "Lithospheric mantle 2"
matLM2.color      = [0.0, 0.7451e0, 0.5882e0]
matLM2.mechanicalProperties = mechanicalProperties()
matLM2.mechanicalProperties.flowlaw    = copy.deepcopy(WetOl)
matLM2.mechanicalProperties.flowlaw.Fc = 5.0
matLM2.mechanicalProperties.plasticity = copy.deepcopy(mohrCoulomb_mantle)
matLM2.densityComposition = densityComposition()
matLM2.densityComposition.beta       = 8.0e-11
matLM2.densityComposition.rho        = 3300
matLM2.densityComposition.deltarho   = 0.0e0
matLM2.densityComposition.T0         = 0.0e0
matLM2.thermalProperties = thermalProperties()
matLM2.thermalProperties.alpha      = 3.1e-5 
matLM2.thermalProperties.k          = 2.25
matLM2.thermalProperties.k_type     = 0
matLM2.thermalProperties.H          = 0.0
matLM2.thermalProperties.D          = 1.0e-6
matLM2.thermalProperties.Cp         = matLM2.thermalProperties.k / (matLM2.densityComposition.rho*matLM2.thermalProperties.D)
matLM2.thermalProperties.T0         = 0.0e0
matLM2.thermalProperties.alpha1    = 3.0e-5 
matLM2.thermalProperties.alpha2    = 3.5e-5 
matLM2.thermalProperties.alphaT1   = 500.0e0
matLM2.thermalProperties.alphaT2   = 2000.0
matLM2.thermalProperties.alphaP1   = 1e5
matLM2.thermalProperties.alphaP2   = 19.8e9
matLM2.thermalProperties.alphaP1_f = 1.0e0
matLM2.thermalProperties.alphaP2_f = 0.5e0
matLM2.densityComposition.use_tidv  = 1
matLM2.densityComposition.perplex_name = None
matLM2.densityComposition.perplex      = None


matLM3            = material_object()
matLM3.name       = "Lithospheric mantle 3"
matLM3.color      = [0.0, 0.7451e0, 0.5882e0]
matLM3.mechanicalProperties = mechanicalProperties()
matLM3.mechanicalProperties.flowlaw    = copy.deepcopy(WetOl)
matLM3.mechanicalProperties.flowlaw.Fc = 5.0
matLM3.mechanicalProperties.plasticity = copy.deepcopy(mohrCoulomb_mantle)
matLM3.densityComposition = densityComposition()
matLM3.densityComposition.beta       = 8.0e-11
matLM3.densityComposition.rho        = 3300
matLM3.densityComposition.deltarho   = 0.0e0
matLM3.densityComposition.T0         = 0.0e0
matLM3.thermalProperties = thermalProperties()
matLM3.thermalProperties.alpha      = 3.1e-5 
matLM3.thermalProperties.k          = 2.25
matLM3.thermalProperties.k_type     = 0
matLM3.thermalProperties.H          = 0.0
matLM3.thermalProperties.D          = 1.0e-6
matLM3.thermalProperties.Cp         = matLM3.thermalProperties.k / (matLM3.densityComposition.rho*matLM3.thermalProperties.D)
matLM3.thermalProperties.T0         = 0.0e0
matLM3.thermalProperties.alpha1    = 3.0e-5 
matLM3.thermalProperties.alpha2    = 3.5e-5 
matLM3.thermalProperties.alphaT1   = 500.0e0
matLM3.thermalProperties.alphaT2   = 2000.0
matLM3.thermalProperties.alphaP1   = 1e5
matLM3.thermalProperties.alphaP2   = 19.8e9
matLM3.thermalProperties.alphaP1_f = 1.0e0
matLM3.thermalProperties.alphaP2_f = 0.5e0
matLM3.densityComposition.use_tidv  = 1
matLM3.densityComposition.perplex_name = None
matLM3.densityComposition.perplex      = None


matSLMd            = material_object()
matSLMd.name       = "Sub-Lithospheric mantle (depleted)"
matSLMd.color      = [1.0e0, 1.0e0, 0.0]
matSLMd.mechanicalProperties = mechanicalProperties()
matSLMd.mechanicalProperties.flowlaw    = copy.deepcopy(WetOl)
matSLMd.mechanicalProperties.flowlaw.Fc = 1.0
matSLMd.mechanicalProperties.plasticity = copy.deepcopy(mohrCoulomb_mantle)
matSLMd.densityComposition = densityComposition()
matSLMd.densityComposition.beta       = 8.0e-11
matSLMd.densityComposition.rho        = 3300
matSLMd.densityComposition.deltarho   = 0.0e0
matSLMd.densityComposition.T0         = 0.0e0
matSLMd.thermalProperties = thermalProperties()
matSLMd.thermalProperties.alpha      = 3.1e-5 
matSLMd.thermalProperties.k          = 2.25
matSLMd.thermalProperties.k_type     = 0
matSLMd.thermalProperties.H          = 0.0
matSLMd.thermalProperties.D          = 1.0e-6
matSLMd.thermalProperties.Cp         = matSLMd.thermalProperties.k / (matSLMd.densityComposition.rho*matSLMd.thermalProperties.D)
matSLMd.thermalProperties.T0         = 0.0e0
matSLMd.thermalProperties.alpha1    = 3.0e-5 
matSLMd.thermalProperties.alpha2    = 3.5e-5 
matSLMd.thermalProperties.alphaT1   = 500.0e0
matSLMd.thermalProperties.alphaT2   = 2000.0
matSLMd.thermalProperties.alphaP1   = 1e5
matSLMd.thermalProperties.alphaP2   = 19.8e9
matSLMd.thermalProperties.alphaP1_f = 1.0e0
matSLMd.thermalProperties.alphaP2_f = 0.5e0
matSLMd.densityComposition.use_tidv  = 1
matSLMd.densityComposition.perplex_name = None
matSLMd.densityComposition.perplex      = None


matSLM            = material_object()
matSLM.name       = "Sub-Lithospheric mantle"
matSLM.color      = [1.0e0, 1.0e0, 0.0]
matSLM.mechanicalProperties = mechanicalProperties()
matSLM.mechanicalProperties.flowlaw    = copy.deepcopy(WetOl)
matSLM.mechanicalProperties.flowlaw.Fc = 1.0
matSLM.mechanicalProperties.plasticity = copy.deepcopy(mohrCoulomb_mantle)
matSLM.densityComposition = densityComposition()
matSLM.densityComposition.beta       = 8.0e-11
matSLM.densityComposition.rho        = 3300
matSLM.densityComposition.deltarho   = 0.0e0
matSLM.densityComposition.T0         = 0.0e0
matSLM.thermalProperties = thermalProperties()
matSLM.thermalProperties.alpha      = 3.1e-5 
matSLM.thermalProperties.k          = 2.25
matSLM.thermalProperties.k_type     = 0
matSLM.thermalProperties.H          = 0.0
matSLM.thermalProperties.D          = 1.0e-6
matSLM.thermalProperties.Cp         = matSLM.thermalProperties.k / (matSLM.densityComposition.rho*matSLM.thermalProperties.D)
matSLM.thermalProperties.T0         = 0.0e0
matSLM.thermalProperties.alpha1    = 3.0e-5 
matSLM.thermalProperties.alpha2    = 3.5e-5 
matSLM.thermalProperties.alphaT1   = 500.0e0
matSLM.thermalProperties.alphaT2   = 2000.0
matSLM.thermalProperties.alphaP1   = 1e5
matSLM.thermalProperties.alphaP2   = 19.8e9
matSLM.thermalProperties.alphaP1_f = 1.0e0
matSLM.thermalProperties.alphaP2_f = 0.5e0
matSLM.densityComposition.use_tidv  = 1
matSLM.densityComposition.perplex_name = None
matSLM.densityComposition.perplex      = None
