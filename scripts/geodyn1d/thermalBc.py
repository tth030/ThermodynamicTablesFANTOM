
class thermal_boundary_conditions(object):
    # general description
    # name             = None
    # temp_bottom      = None
    # temp_top         = None
    # q_bottom         = None
    # q_top            = None
    # temp_potential   = None
    # adiabat          = False
    
    def __init__(self,name             = None,
                      temp_bottom      = None,
                      temp_top         = None,
                      q_bottom         = None,
                      q_top            = None,
                      temp_potential   = None,
                      adiabat          = False):
        self.name           = name
        self.temp_bottom    = temp_bottom
        self.temp_top       = temp_top
        self.q_bottom       = q_bottom
        self.q_top          = q_top
        self.temp_potential = temp_potential
        self.adiabat        = adiabat
        
def update_thermal_Bc(thermBc,specific_book_properties):
    book = thermBc.__dict__
    for specific_key in specific_book_properties.keys():
        for key in book.keys():
            if ( specific_key == key ):
                setattr(thermBc,key,specific_book_properties[key])
                break
    return

thermBcUC      = thermal_boundary_conditions()
thermBcUC.name = 'Thermal boundary conditions Upper-Crust' 
thermBcUC.temp_top = 273.15e0

thermBcMC      = thermal_boundary_conditions()
thermBcMC.name = 'Thermal boundary conditions Middle-Crust' 

thermBcLC      = thermal_boundary_conditions()
thermBcLC.name = 'Thermal boundary conditions Lower-Crust' 

thermBcLM1      = thermal_boundary_conditions()
thermBcLM1.name = 'Thermal boundary conditions Upper Lithospheric mantle' 

thermBcLM2      = thermal_boundary_conditions()
thermBcLM2.name = 'Thermal boundary conditions Medium Lithospheric mantle' 

thermBcLM3      = thermal_boundary_conditions()
thermBcLM3.name = 'Thermal boundary conditions Lower Lithospheric mantle'

thermBcSLMd         = thermal_boundary_conditions()
thermBcSLMd.name    = 'Thermal boundary conditions Sub-Lithospheric mantle (depleted)'
thermBcSLMd.adiabat = True

thermBcSLM      = thermal_boundary_conditions()
thermBcSLM.name = 'Thermal boundary conditions Sub-Lithospheric mantle' 
thermBcSLM.temp_bottom    = 1793.15e0
thermBcSLM.temp_potential = 1553.15e0
thermBcSLM.q_bottom       = 19.5e-3
thermBcSLM.adiabat        = True




