#!/usr/bin/env python
#
# Catalog of lithospheres with updated parameters
#

from .geodyn1d import flowlaws

######## lith125 ##########

book_lith125_GL = {
    'matUC': {
        'rho': 2860e0,
        'H': 0.9065e-6,
        'alpha': 3.0e-5,
        'Fc': 1
     },
    'matMC': {
        'rho': 2860e0,
        'H': 0.9065e-6,
        'alpha': 3.0e-5,
        'Fc': 1
     },
     'matLC': {
         'rho': 2860e0,
         'H': 0.9065e-6,
         'alpha': 3.0e-5,
         'flowlaw': flowlaws.Dia01,
         'Fc': 0.1
    },
    'matLM1': {'rho': 3311,
               'deltarho': -20.0,
               'use_tidv': 1,
               'k': 2.25,
               'H': 0.0,
               'k_type': 0,
               'alpha': 3.0e-5,
               'Fc': 5,
               'perplex_name': 'slm_sp2008'
              },
    'matLM2': {'rho': 3311,
               'deltarho': -20.0,
               'use_tidv': 1,
               'k': 2.25,
               'H': 0.0,
               'k_type': 0,
               'alpha': 3.0e-5,
               'Fc': 5,
               'perplex_name': 'slm_sp2008'
              },
    'matLM3': {'rho': 3311,
               'deltarho': -20.0,
               'use_tidv': 1,
               'k': 2.25,
               'H': 0.0,
               'k_type': 0,
               'alpha': 3.0e-5,
               'Fc': 5,
               'perplex_name': 'slm_sp2008'
              },
    'matSLM': {'rho': 3311,
               'deltarho': 0.0,
               'use_tidv': 1,
               'k': 2.25,
               'k_type': 0,
               'alpha': 3.0e-5,
               'Fc': 1,
               'perplex_name': 'slm_sp2008'
              },
    'thermBcSLM': {
        'temp_bottom': 1793.15e0,
        'temp_potential': 1553.15e0,
        'q_bottom': 19.5e-3
    }

}

######## lith125 ##########

book_lith125_GL_highTp = {
    'matUC': {
        'rho': 2860e0,
        'H': 0.9065e-6,
        'alpha': 3.0e-5,
        'Fc': 1
     },
    'matMC': {
        'rho': 2860e0,
        'H': 0.9065e-6,
        'alpha': 3.0e-5,
        'Fc': 1
     },
     'matLC': {
         'rho': 2860e0,
         'H': 0.9065e-6,
         'alpha': 3.0e-5,
         'flowlaw': flowlaws.Dia01,
         'Fc': 0.1
    },
    'matLM1': {'rho': 3311,
               'deltarho': -20.0,
               'use_tidv': 1,
               'k': 2.25,
               'H': 0.0,
               'k_type': 0,
               'alpha': 3.0e-5,
               'Fc': 5,
               'perplex_name': 'slm_sp2008'
              },
    'matLM2': {'rho': 3311,
               'deltarho': -20.0,
               'use_tidv': 1,
               'k': 2.25,
               'H': 0.0,
               'k_type': 0,
               'alpha': 3.0e-5,
               'Fc': 5,
               'perplex_name': 'slm_sp2008'
              },
    'matLM3': {'rho': 3311,
               'deltarho': -20.0,
               'use_tidv': 1,
               'k': 2.25,
               'H': 0.0,
               'k_type': 0,
               'alpha': 3.0e-5,
               'Fc': 5,
               'perplex_name': 'slm_sp2008'
              },
    'matSLM': {'rho': 3311,
               'deltarho': 0.0,
               'use_tidv': 1,
               'k': 2.25,
               'k_type': 0,
               'alpha': 3.0e-5,
               'Fc': 1,
               'perplex_name': 'slm_sp2008'
              },
    'thermBcSLM': {
        'temp_bottom': 1843.15e0,
        'temp_potential': 1603.15e0,
        'q_bottom': 19.5e-3
    }

}

######## lith125 ##########

book_lith125_GL_lowTp = {
    'matUC': {
        'rho': 2860e0,
        'H': 0.9065e-6,
        'alpha': 3.0e-5,
        'Fc': 1
     },
    'matMC': {
        'rho': 2860e0,
        'H': 0.9065e-6,
        'alpha': 3.0e-5,
        'Fc': 1
     },
     'matLC': {
         'rho': 2860e0,
         'H': 0.9065e-6,
         'alpha': 3.0e-5,
         'flowlaw': flowlaws.Dia01,
         'Fc': 0.1
    },
    'matLM1': {'rho': 3311,
               'deltarho': -20.0,
               'use_tidv': 1,
               'k': 2.25,
               'H': 0.0,
               'k_type': 0,
               'alpha': 3.0e-5,
               'Fc': 5,
               'perplex_name': 'slm_sp2008'
              },
    'matLM2': {'rho': 3311,
               'deltarho': -20.0,
               'use_tidv': 1,
               'k': 2.25,
               'H': 0.0,
               'k_type': 0,
               'alpha': 3.0e-5,
               'Fc': 5,
               'perplex_name': 'slm_sp2008'
              },
    'matLM3': {'rho': 3311,
               'deltarho': -20.0,
               'use_tidv': 1,
               'k': 2.25,
               'H': 0.0,
               'k_type': 0,
               'alpha': 3.0e-5,
               'Fc': 5,
               'perplex_name': 'slm_sp2008'
              },
    'matSLM': {'rho': 3311,
               'deltarho': 0.0,
               'use_tidv': 1,
               'k': 2.25,
               'k_type': 0,
               'alpha': 3.0e-5,
               'Fc': 1,
               'perplex_name': 'slm_sp2008'
              },
    'thermBcSLM': {
        'temp_bottom': 1743.15e0,
        'temp_potential': 1503.15e0,
        'q_bottom': 19.5e-3
    }

}


########## lith 160 #######

book_lith160_GL = {
    'thicknesses': [15e3,10e3,10e3,65e3,35e3,25e3,440e3],
    'matUC': {
        'H': 0.9065e-6,
        'rho':2860,
        'alpha': 3.0e-5,
        'k_type': 0,
        'Fc': 0.05 #0.5
     },
    'matMC': {
        'H': 0.9065e-6,
        'rho':2860,
        'alpha': 3.0e-5,
        'k_type': 0,
        'Fc': 0.05 #0.5
     },
     'matLC': {
         'H': 0.9065e-6,
         'rho':2860,
         'alpha': 3.0e-5,
         'k_type': 0,
         'flowlaw': flowlaws.Dia01,
         'Fc': 0.05
    },
    'matLM1': {'use_tidv': 1,
               'rho': 3311,
               'deltarho': -28.5,
               'alpha': 3.0e-5,
               'k': 2.25,
               'H': 0.0,
               'k_type': 0,
               'Fc': 0.5,
               'perplex_name': 'slm_sp2008'
              },
    'matLM2': {'use_tidv': 1,
               'rho': 3311,
               'deltarho': -28.5,
               'alpha': 3.0e-5,
               'k': 2.25,
               'H': 0.0,
               'k_type': 0,
               'Fc': 5,
               'perplex_name': 'slm_sp2008'
              },
    'matLM3': {'use_tidv': 1,
               'rho': 3311,
               'deltarho': -28.5,
               'alpha': 3.0e-5,
               'k': 2.25,
               'H': 0.0,
               'k_type': 0,
               'Fc': 5,
               'perplex_name': 'slm_sp2008'
              },
    'matSLM': {'use_tidv': 1,
               'rho': 3311,
               'deltarho': 0,
               'alpha': 3.0e-5,
               'k': 2.25,
               'k_type': 0,
               'Fc': 1,
               'perplex_name': 'slm_sp2008'
              },
    'thermBcSLM': {
        'temp_bottom': 1793.15e0,
        'temp_potential': 1553.15e0,
         'q_bottom': 15.4316406e-3
     }

}


######## lith180 #########

book_lith180_GL = {
    'thicknesses': [15e3,10e3,10e3,45e3,45e3,55e3,420e3],
    'matUC': {
        'H': 0.9065e-6,
        'rho':2860,
        'alpha':3.0e-5,
        'k_type': 0,
        'Fc': 0.05 #0.5
     },
    'matMC': {
        'H': 0.9065e-6,
        'rho':2860,
        'alpha':3.0e-5,
        'k_type': 0,
        'Fc': 0.05 #0.5
     },
     'matLC': {
         'H': 0.9065e-6,
         'rho':2860,
         'alpha':3.0e-5,
         'k_type': 0,
         'flowlaw': flowlaws.Dia01,
         'Fc': 0.05 #0.015
    },
    'matLM1': {'use_tidv': 1,
               'rho': 3311,
               'deltarho': -31.5,
               'k': 2.25,
               'H': 0.0,
               'alpha':3.0e-5,
               'k_type': 0,
               'Fc': 0.5,
               'perplex_name': 'slm_sp2008'
              },
    'matLM2': {'use_tidv': 1,
               'rho': 3311,
               'deltarho': -31.5,
               'k': 2.25,
               'H': 0.0,
               'alpha':3.0e-5,
               'k_type': 0,
               'Fc': 5,
               'perplex_name': 'slm_sp2008'
              },
    'matLM3': {
               'use_tidv': 1,
               'rho': 3311,
               'deltarho': -31.5,
               'k': 2.25,
               'H': 0.0,
               'alpha':3.0e-5,
               'k_type': 0,
               'Fc': 5,
               'perplex_name': 'slm_sp2008'
              },
    'matSLM': {
               'rho': 3311,
               'deltarho': 0,
               'use_tidv': 1,
               'k': 2.25,
               'alpha':3.0e-5,
               'k_type': 0,
               'Fc': 1,
               'perplex_name': 'slm_sp2008'
              },
    'thermBcSLM': {
        'temp_bottom': 1793.15e0,
        'temp_potential': 1553.15e0,
         'q_bottom': 13.8183593e-3
     }

}

########## lith200 ############

book_lith200_GL = {
    'thicknesses': [15e3,10e3,10e3,45e3,45e3,75e3,400e3],
    'matUC': {
        'rho': 2860,
        'H': 0.9065e-6,
        'alpha': 3.0e-5,
        'Fc': 0.02   #0.2
     },
    'matMC': {
        'rho': 2860,
        'H': 0.9065e-6,
        'alpha': 3.0e-5,
        'Fc': 0.02  #0.2
     },
     'matLC': {
         'rho': 2860,
         'H': 0.9065e-6,
         'alpha': 3.0e-5,
         'flowlaw': flowlaws.Dia01,
         'Fc': 0.02
    },
    'matLM1': {'rho': 3311,
               'deltarho': -34,
               'use_tidv': 1,
               'k': 2.25,
               'alpha': 3.0e-5,
               'H':0.0,
               'k_type': 0,
               'internalAngleOfFriction1': 15,
               'Fc': 0.2,
               'perplex_name': 'slm_sp2008'
              },
    'matLM2': {'rho': 3311,
               'deltarho': -34,
               'use_tidv': 1,
               'k': 2.25,
               'alpha': 3.0e-5,
               'H':0.0,
               'k_type': 0,
               'Fc': 1.5,
               'perplex_name': 'slm_sp2008'
              },
    'matLM3': {'rho': 3311,
               'deltarho': -34,
               'use_tidv': 1,
               'k': 2.25,
               'alpha': 3.0e-5,
               'H':0.0,
               'k_type': 0,
               'Fc': 1.5,
               'perplex_name': 'slm_sp2008'
              },
    'matSLM': {'rho': 3311,
               'deltarho': 0,
               'use_tidv': 1,
               'alpha': 3.0e-5,
               'H':0.0,
               'k_type': 0,
               'Fc': 1,
               'perplex_name': 'slm_sp2008'
              },
    'thermBcSLM': {
        'temp_bottom': 1793.15e0,
        'temp_potential': 1553.15e0,
        'q_bottom': 12.5253906e-3
    }
}

########## lith240 ############

book_lith240_GL = {
    'thicknesses': [15e3,10e3,10e3,80e3,10e3,115e3,360e3],
    'matUC': {'rho': 2860,
        'H': 0.9065e-6,
        'alpha': 3.0e-5,
        'Fc': 0.05   #0.2
     },
    'matMC': {
        'rho': 2860,
        'H': 0.9065e-6,
        'alpha': 3.0e-5,
        'Fc': 0.05  #0.2
     },
     'matLC': {
         'rho':2860,
         'H': 0.9065e-6,
         'alpha': 3.0e-5,
         'flowlaw': flowlaws.Dia01,
         'Fc': 0.05
    },
    'matLM1': {'deltarho':-37.5,
               'rho': 3311,
               'use_tidv': 1,
               'alpha': 3.0e-5,
               'k': 2.25,
               'H':0.0,
               'k_type': 0,
               'Fc': 0.5,
               'perplex_name': 'slm_sp2008'
              },
    'matLM2': {'deltarho':-37.5,
               'rho': 3311,
               'use_tidv': 1,
               'alpha': 3.0e-5,
               'k': 2.25,
               'H':0.0,
               'k_type': 0,
               'Fc': 5,
               'perplex_name': 'slm_sp2008'
              },
    'matLM3': {'deltarho':-37.5,
               'rho': 3311,
               'use_tidv': 1,
               'alpha': 3.0e-5,
               'k': 2.25,
               'H':0.0,
               'k_type': 0,
               'Fc': 5,
               'perplex_name': 'slm_sp2008'
              },
    'matSLM': {'deltarho':0,
               'rho': 3311,
               'use_tidv': 1,
               'alpha': 3.0e-5,
               'H':0.0,
               'k_type': 0,
               'Fc': 1,
               'perplex_name': 'slm_sp2008'
              },
    'thermBcSLM': {
        'temp_bottom': 1793.15e0,
        'temp_potential': 1553.15e0,
        'q_bottom': 10.5878906e-3,
    }
}

############ lith280 ##############

book_lith280_GL = {
    'thicknesses': [15e3,10e3,10e3,80e3,10e3,155e3,320e3],
    'matUC': {'rho': 2860,
        'H': 0.9065e-6,
        'alpha': 3.0e-5,
        'Fc': 0.05   #0.2
     },
    'matMC': {
        'rho': 2860,
        'H': 0.9065e-6,
        'alpha': 3.0e-5,
        'Fc': 0.05
     },
     'matLC': {
         'rho':2860,
         'H': 0.9065e-6,
         'alpha': 3.0e-5,
         'flowlaw': flowlaws.Dia01,
         'Fc': 0.05
    },
    'matLM1': {'deltarho':-40,
               'rho': 3311,
               'use_tidv': 1,
               'alpha': 3.0e-5,
               'k': 2.25,
               'H':0.0,
               'k_type': 0,
               'Fc': 0.5,
               'perplex_name': 'slm_sp2008'
              },
    'matLM2': {'deltarho':-40,
               'rho': 3311,
               'use_tidv': 1,
               'alpha': 3.0e-5,
               'k': 2.25,
               'H':0.0,
               'k_type': 0,
               'Fc': 5,
               'perplex_name': 'slm_sp2008'
              },
    'matLM3': {'deltarho':-40,
               'rho': 3311,
               'use_tidv': 1,
               'alpha': 3.0e-5,
               'k': 2.25,
               'H':0.0,
               'k_type': 0,
               'Fc': 5,
               'perplex_name': 'slm_sp2008'
              },
    'matSLM': {'deltarho':0,
               'rho': 3311,
               'use_tidv': 1,
               'alpha': 3.0e-5,
               'H':0.0,
               'k_type': 0,
               'Fc': 1,
               'perplex_name': 'slm_sp2008'
              },
    'thermBcSLM': {
        'temp_bottom': 1793.15e0,
        'temp_potential': 1553.15e0,
        'q_bottom': 9.2031250e-3,
    }
}


###### ridge #######


book_ridge   = {
               'matUC': {
                   'H': 0.0,
                   'alpha': 3.0e-5,
                   'rho': 2900
                },
               'matMC': {
                   'H': 0.0,
                   'alpha': 3.0e-5,
                   'rho': 2900
                },
                'matLC': {
                    'H': 0.0,
                    'alpha': 3.0e-5,
                    'rho': 2900
               },
               'matSLMd': {
                          'deltarho': -12,
                          'rho': 3311,
                          'use_tidv': 1,
                          'alpha': 3.0e-5,
                          'perplex_name': 'slm_sp2008'
                         },
              'matSLM': {
                          'deltarho': 0,
                          'rho': 3311,
                          'use_tidv': 1,
                          'alpha': 3.0e-5,
                          'perplex_name': 'slm_sp2008'
                         },
               'thermBcSLMd': {                  # applied LM1 is adiabatic
                   'temp_bottom': 1603.15e0,
                   'temp_potential': 1553.15e0,
                   'q_bottom': 437.2421875e-3
               },
               'thermBcSLM': {
                   'temp_bottom': 1793.15e0,
                   'temp_potential': 1553.15e0,
                   'q_bottom': 437.2421875e-3
               }
               }

