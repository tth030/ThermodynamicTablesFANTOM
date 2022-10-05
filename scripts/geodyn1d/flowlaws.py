
class flowlaw(object):
    A = None
    n = None
    Ea = None
    Va = None
    Fc = None
    name_flowlaw = None

WetQz = flowlaw()
WetQz.name_flowlaw = "Wet Quartz"
WetQz.A = 8.574e-28
WetQz.n = 4.0
WetQz.Ea = 223e3
WetQz.Va = 0.
WetQz.Fc = 1.

WetOl = flowlaw()
WetOl.name_flowlaw = "Wet Olivine"
WetOl.A = 1.7578e-14
WetOl.n = 3.0
WetOl.Ea = 430e3
WetOl.Va = 15e-6
WetOl.Fc = 1.

Dia01 = flowlaw()
Dia01.name_flowlaw = "DRY DIABASE (Maryland) 0.1"
Dia01.A = 5.77904e-27
Dia01.n = 4.7
Dia01.Ea = 485.0e3
Dia01.Va = 0.
Dia01.Fc = 0.1

