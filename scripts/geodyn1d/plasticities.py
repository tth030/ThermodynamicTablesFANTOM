
# class plasticity(object):
#     kind = None

# class MohrCoulombCriterion(object):
#     name_plasticity = 'Mohr-Coulomb'
#     internalAngleOfFriction1 = None
#     cohesion1                = None
#     frictionalWeakening     = None
#     cohesionWeakening       = None
#     internalAngleOfFriction2 = None
#     cohesion2 = None
#     strain_phi_1 = None
#     strain_phi_2 = None
#     strain_C_1 = None
#     strain_C_2 = None

# mohrCoulomb           = plasticity()
# mohrCoulomb.kind      = MohrCoulombCriterion()
# mohrCoulomb.kind.name_plasticity          = MohrCoulombCriterion().name_plasticity
# mohrCoulomb.kind.internalAngleOfFriction1 = 15
# mohrCoulomb.kind.cohesion1                = 20e6
# mohrCoulomb.kind.frictionalWeakening      = True
# mohrCoulomb.kind.cohesionWeakening        = True
# mohrCoulomb.kind.internalAngleOfFriction2 = 2
# mohrCoulomb.kind.cohesion2                = 4e6
# mohrCoulomb.kind.strain_phi_1             = 0.05
# mohrCoulomb.kind.strain_phi_2             = 1.05
# mohrCoulomb.kind.strain_C_1               = 0.05
# mohrCoulomb.kind.strain_C_2               = 1.05

#class plasticity(object):
#    pass
#
#class MohrCoulombCriterion(plasticity):
class plasticity(object):
    name_plasticity = None # 'Mohr-Coulomb'
    internalAngleOfFriction1 = None
    cohesion1                = None
    frictionalWeakening     = None
    cohesionWeakening       = None
    internalAngleOfFriction2 = None
    PowerLawBreakdown = None
    cohesion2 = None
    strain_phi_1 = None
    strain_phi_2 = None
    strain_C_1 = None
    strain_C_2 = None

mohrCoulomb           = plasticity()
#mohrCoulomb      = MohrCoulombCriterion()
mohrCoulomb.name_plasticity          = 'Mohr-Coulomb' #MohrCoulombCriterion().name_plasticity
mohrCoulomb.internalAngleOfFriction1 = 15
mohrCoulomb.cohesion1                = 20e6
mohrCoulomb.frictionalWeakening      = True
mohrCoulomb.cohesionWeakening        = True
mohrCoulomb.internalAngleOfFriction2 = 2
mohrCoulomb.PowerLawBreakdown        = 300e6
mohrCoulomb.cohesion2                = 4e6
mohrCoulomb.strain_phi_1             = 0.05
mohrCoulomb.strain_phi_2             = 1.05
mohrCoulomb.strain_C_1               = 0.05
mohrCoulomb.strain_C_2               = 1.05

mohrCoulomb_mantle           = plasticity()
#mohrCoulomb_mantle      = MohrCoulombCriterion()
mohrCoulomb_mantle.name_plasticity          = 'Mohr-Coulomb' #MohrCoulombCriterion().name_plasticity
mohrCoulomb_mantle.internalAngleOfFriction1 = 15
mohrCoulomb_mantle.cohesion1                = 20e6
mohrCoulomb_mantle.frictionalWeakening      = True
mohrCoulomb_mantle.cohesionWeakening        = True
mohrCoulomb_mantle.internalAngleOfFriction2 = 4
mohrCoulomb_mantle.PowerLawBreakdown        = 300e6
mohrCoulomb_mantle.cohesion2                = 20e6
mohrCoulomb_mantle.strain_phi_1             = 0.05
mohrCoulomb_mantle.strain_phi_2             = 1.05
mohrCoulomb_mantle.strain_C_1               = 0.05
mohrCoulomb_mantle.strain_C_2               = 1.05
