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

from .material import *
from .thermalBc import *

class layer (object):
    
    __slots__ = ['layerid','name_material','material','thickness', 'thermalBc','specific_book_properties','perplex','f_interp','T_top','T_down','f_interp_melt_fraction']

    def __init__(self,layerid,name_material,thickness,thermalBc,specific_book_properties=None):
        self.layerid = layerid
        self.name_material = name_material
        self.thickness = thickness
        if ( specific_book_properties ):
            if ( self.name_material ) in specific_book_properties.keys():
                specific_material_properties = specific_book_properties[self.name_material]
            else:
                specific_material_properties = None
            if ( specific_material_properties ):
                #TODO UGLY !!!! SHOULD BE UPDATED...
                book_init = get_book_material(eval(self.name_material))
                update_material(eval(self.name_material),self.name_material,specific_material_properties,printOut=False)
                
        self.material  = material_object()
        self.material.__dict__ = eval(self.name_material).__dict__.copy()
                
        self.material.mechanicalProperties = mechanicalProperties()
        self.material.densityComposition = densityComposition()
        self.material.thermalProperties = thermalProperties()
        
        self.material.mechanicalProperties.__dict__ = eval(self.name_material).mechanicalProperties.__dict__.copy()
        self.material.densityComposition.__dict__ = eval(self.name_material).densityComposition.__dict__.copy()
        self.material.thermalProperties.__dict__ = eval(self.name_material).thermalProperties.__dict__.copy()

        self.thermalBc = thermal_boundary_conditions()
        self.thermalBc.__dict__ = eval(thermalBc).__dict__.copy()
                
        if ( specific_book_properties ):
            #TODO UGLY !!!! SHOULD BE UPDATED... 
            if ( specific_material_properties ):
                update_material(eval(self.name_material),self.name_material,book_init,printOut=False)
                #material.mechanicalProperties
                if ( 'flowlaw' in specific_material_properties ):
                    self.material.mechanicalProperties.flowlaw = copy.deepcopy(specific_material_properties['flowlaw'])
                if ( 'Fc' in specific_material_properties ):
                    self.material.mechanicalProperties.flowlaw.Fc = specific_material_properties['Fc']
                if ( 'plasticity' in specific_material_properties ):
                    self.material.mechanicalProperties.flowlaw = copy.deepcopy(specific_material_properties['plasticity'])
                if ( 'internalAngleOfFriction1' in specific_material_properties ):
                    self.material.mechanicalProperties.plasticity.internalAngleOfFriction1 = specific_material_properties['internalAngleOfFriction1']
            if ( 'thicknesses' in specific_book_properties ):
                    self.thickness = specific_book_properties['thicknesses'][self.layerid]

            # !!!!
            if ( thermalBc ) in specific_book_properties.keys():
                specific_thermalBc_properties = specific_book_properties[thermalBc]
            else:
                specific_thermalBc_properties = None
            if ( specific_thermalBc_properties ):
                update_thermal_Bc(self.thermalBc,specific_thermalBc_properties)
        
        self.perplex  = None
        self.f_interp = None
        self.f_interp_melt_fraction = None
        self.T_top    = None
        self.T_down   = None
  

    def display_attributes(self,object):
        attrs = vars(object)
        print(', \n'.join("%s: %s" % item for item in attrs.items()))
        return

    def display(self):
        self.display_attributes(self.material)
        print('Mechanical Properties ---------------------')
        self.display_attributes(self.material.mechanicalProperties)
        self.display_attributes(self.material.mechanicalProperties.flowlaw)
        self.display_attributes(self.material.mechanicalProperties.plasticity)
        print('Density Composition -----------------------')
        self.display_attributes(self.material.densityComposition)
        print('Thermal Properties  -----------------------')
        self.display_attributes(self.material.thermalProperties)
        print('Thermal Boundary Conditions ---------------')
        self.display_attributes(self.thermalBc)
        return
        
    def get_book_material(self):
        book = get_book_material(self.material)
        return book


