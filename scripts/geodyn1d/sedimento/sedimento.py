#!/usr/bin/env python 
"""
  This module allows computing compaction for silico-clastic sediments
 
Releases:
    0.1: sedimento.py (included in pygeodyn1d package)

Req:           None
Version:       0.1
Date:          2020/05/26
Authors:       Thomas Theunissen, UiB 2020 version 0.0
Copyright:     Creative Commons CC0 - http://creativecommons.org/publicdomain/zero/1.0/legalcode
"""

import sys,os,time
import numpy as np

# Global variables

sedimento_surface_porosity_factor = 0.52e0   # no units :: sediment compaction following Athy, 1930 equation
sedimento_compaction_coefficient  = 4.7e-4   # m-1      :: sediment compaction following Athy, 1930 equation
rho_water                         = 1027.0e0 # ocean water density Kg/m3 | already defined in core.py
rho_salt                          = 2300     # mixed evapories 2150 = halite
rho_silicoclasticGrain            = 2700.0e0
dz                                = 10       # (m) vertical resolution

def main():
    """ Run outside from pygeodyn1d """

    # average values on distal margins
    # Gabon
    input_compactedthickness = 5600
    salt_thickness           = 600
    postrift_thickness       = 2500
    water_thickness          = 2500

    #South Congo
    input_compactedthickness = 600
    salt_thickness           = 700
    postrift_thickness       = 5500
    water_thickness          = 1500

    #Kwanza
    input_compactedthickness = 900
    salt_thickness           = 1500
    postrift_thickness       = 2000
    water_thickness          = 2500

    sed1                     = sedlayer(sedtype='silicoclastic',thickness=input_compactedthickness)
    salt                     = sedlayer(sedtype='salt',thickness=salt_thickness)
    postrift                 = sedlayer(sedtype='silicoclastic',thickness=postrift_thickness)
    water                    = sedlayer(sedtype='water',thickness=water_thickness)

    sed1.addLoadingLayers([water,postrift,salt]) # from top to bottom
    #sed1.removeLoadingLayers(-1)

    print("Compacted Thickness      = %s"%(input_compactedthickness))    
    print("UncompactedThickness     = %s"%(sed1.getUncompactedThickness()))
    print("SelfUncompactedThickness = %s"%(sed1.getSelfUnCompactedThickness(sed1.thickness,sed1.rho)))


#    #------------------------------------------
#    input_compactedthickness = 4580 # uncompacted thickness
#    sed1                     = sedlayer(sedtype='silicoclastic',thickness=input_compactedthickness,IsCompactedThickness=False,loadinglayers=[water,postrift,salt])
#    #salt                     = sedlayer(sedtype='salt',thickness=2000)
#    #postrift                 = sedlayer(sedtype='silicoclastic',thickness=4000)
#    #water                    = sedlayer(sedtype='water',thickness=1000)
#
#    #sed1.addLoadingLayers([water,postrift,salt]) # from top to bottom
#    #sed1.removeLoadingLayers(-1)
#
#    print("Compacted Thickness      = %s"%(sed1.thickness))
#    print("UncompactedThickness     = %s"%(sed1.getUncompactedThickness()))
#    print("SelfUncompactedThickness = %s"%(sed1.getSelfUnCompactedThickness(sed1.thickness,sed1.rho)))
#
#    #------------------------------------------
#    input_compactedthickness = 3000 # uncompacted thickness
#    salt                     = sedlayer(sedtype='salt',thickness=2500)
#    sed1                     = sedlayer(sedtype='silicoclastic',thickness=input_compactedthickness,IsCompactedThickness=False,loadinglayers=[salt])
#    print("Compacted Thickness      = %s"%(sed1.thickness))
#    print("UncompactedThickness     = %s"%(sed1.getUncompactedThickness()))
#    print("SelfUncompactedThickness = %s"%(sed1.getSelfUnCompactedThickness(sed1.thickness,sed1.rho)))
#

class sedlayer(object):

    __slots__ = ['sedtype','thickness','IsCompactedThickness','uncompactedthickness','selfUncompactedthickness','loadinglayers',
                 'averageDensity','z','rho']

    def __init__(self,sedtype,thickness,IsCompactedThickness=True,loadinglayers=None):
        self.sedtype                    = sedtype
        self.thickness                  = thickness
        self.IsCompactedThickness       = IsCompactedThickness
        self.loadinglayers              = loadinglayers   # list of layers loading the sed. layer
        self.uncompactedthickness       = None   #
        self.selfUncompactedthickness   = None   #
        self.z                          = np.arange(0,self.thickness+dz,dz)
        self.rho                        = self.getRho()
        self.averageDensity             = np.sum(self.rho)/len(self.z)
        #print('%s %f'%(self.sedtype,self.averageDensity))
        if (not self.IsCompactedThickness and self.loadinglayers):
            self.uncompactedthickness   = self.thickness
            self.thickness              = self.getCompactedThickness(self.uncompactedthickness)
            self.z                      = np.arange(0,self.thickness+dz,dz)
            self.rho                    = self.getRho()
            self.averageDensity         = np.sum(self.rho)/len(self.z)
            self.IsCompactedThickness   = True

    def getRho(self):
        if ( self.sedtype == 'silicoclastic' ):
            rho = rho_silicoclasticGrain - ( rho_silicoclasticGrain - rho_water) * sedimento_surface_porosity_factor * np.exp(-1.e0*sedimento_compaction_coefficient*self.z)
        elif ( self.sedtype == 'salt' ):
            rho = np.zeros_like(self.z) ; rho[:] = rho_salt
        elif ( self.sedtype == 'water' ):
            rho = np.zeros_like(self.z) ; rho[:] = rho_water
        return rho

    def getUncompactedThickness(self):
        if (not self.IsCompactedThickness):
            self.uncompactedthickness   = self.thickness
            self.thickness              = self.getCompactedThickness(self.uncompactedthickness)
            self.z                      = np.arange(0,self.thickness+dz,dz)
            self.rho                    = self.getRho()
            self.averageDensity         = np.sum(self.rho)/len(self.z)
            self.IsCompactedThickness   = True
        if ( self.sedtype != 'silicoclastic' ):
            self.uncompactedthickness     = self.thickness
            self.selfUncompactedthickness = self.uncompactedthickness
        else:
            ''' Decompaction layers '''
            self.uncompactedthickness = self.thickness
            if (self.IsCompactedThickness and not self.loadinglayers ):
                """ No loading layer """
                self.uncompactedthickness = self.thickness
            if (self.IsCompactedThickness and self.loadinglayers ):
                ''' Substract Weight from Loading layers '''
                self.uncompactedthickness = self.CalculateUncompactedThickness()
            self.selfUncompactedthickness = self.getSelfUnCompactedThickness(self.uncompactedthickness,self.rho)

        return self.uncompactedthickness

    def getCompactedThickness(self,uncompactedthickness):
        if (self.loadinglayers==None):
            thickness = uncompactedthickness
        else:
            SelfUncompactedLoadingLayerThickness = 0
            for sedimlayer in self.loadinglayers:
                pseudolayer                          = self.getEquivalentSilicoClasticLayer(sedimlayer.averageDensity*sedimlayer.thickness)
                tmp                                  = pseudolayer.getUncompactedThickness()
                SelfUncompactedLoadingLayerThickness = SelfUncompactedLoadingLayerThickness + pseudolayer.selfUncompactedthickness
            #print("Total SelfUncompactedLoadingLayerThickness = %f (getCompactedThickness)"%(SelfUncompactedLoadingLayerThickness))
            rho0        = rho_silicoclasticGrain - ( rho_silicoclasticGrain - rho_water) * sedimento_surface_porosity_factor
            pseudolayer = self.getEquivalentSilicoClasticLayer(rho0*SelfUncompactedLoadingLayerThickness)
            tmp         = pseudolayer.getUncompactedThickness()

            tmplayer    = sedlayer(sedtype='silicoclastic',thickness=uncompactedthickness)
            mass        = tmplayer.averageDensity*uncompactedthickness

            thicknesses = np.arange(0+dz,30000+dz,dz)
            for thickness in thicknesses:
                z   = np.arange(pseudolayer.thickness,pseudolayer.thickness+thickness+dz,dz)
                rho = rho_silicoclasticGrain - ( rho_silicoclasticGrain - rho_water) * sedimento_surface_porosity_factor * np.exp(-1.e0*sedimento_compaction_coefficient*z)
                averageDensity = np.sum(rho)/len(z)
                if (mass<averageDensity*thickness):
                    break
        return thickness-dz

    def addLoadingLayers(self,loadinglayers):
        self.loadinglayers = loadinglayers
        if (not self.IsCompactedThickness):
            self.uncompactedthickness   = self.thickness
            self.thickness              = self.getCompactedThickness(self.uncompactedthickness)
            self.z                      = np.arange(0,self.thickness+dz,dz)
            self.rho                        = self.getRho()
            self.averageDensity             = np.sum(self.rho)/len(self.z)
            self.IsCompactedThickness   = True
        return 1
        
    def removeLoadingLayers(self,index):
        if (index<0):
            self.loadinglayers = None
        else:
            try:
                self.loadinglayers = self.loadinglayers.pop(index)
            except IndexError:
                print("Index Out of Range")

    def getSelfUnCompactedThickness(self,thickness,rho):
        uncompactedthickness = 0
        rho0                 = rho_silicoclasticGrain - ( rho_silicoclasticGrain - rho_water) * sedimento_surface_porosity_factor
        for i in np.arange(0,len(rho)):
            mass                 = rho[i] * dz
            uncompactedthickness = uncompactedthickness + mass/rho0
        return uncompactedthickness

    def CalculateUncompactedThickness(self):
        SelfUncompactedLoadingLayerThickness = 0
        for sedimlayer in self.loadinglayers:
            pseudolayer                          = self.getEquivalentSilicoClasticLayer(sedimlayer.averageDensity*sedimlayer.thickness)
            tmp                                  = pseudolayer.getUncompactedThickness()
            #print("%f %s %s"%(pseudolayer.thickness,pseudolayer.sedtype,pseudolayer.uncompactedthickness))
            SelfUncompactedLoadingLayerThickness = SelfUncompactedLoadingLayerThickness + pseudolayer.selfUncompactedthickness
        rho0        = rho_silicoclasticGrain - ( rho_silicoclasticGrain - rho_water) * sedimento_surface_porosity_factor
        pseudolayer = self.getEquivalentSilicoClasticLayer(rho0*SelfUncompactedLoadingLayerThickness)
        #print("Total SelfUncompactedLoadingLayerThickness = %f (CalculateUncompactedThickness)"%(SelfUncompactedLoadingLayerThickness))
        uncompactedthickness = pseudolayer.getUncompactedThickness()
        #print("%f %s %f %f"%(pseudolayer.thickness,pseudolayer.sedtype,pseudolayer.selfUncompactedthickness,pseudolayer.rho[len(pseudolayer.rho)-1]))

        z   = np.arange(pseudolayer.thickness,pseudolayer.thickness+self.thickness+dz,dz)
        rho = rho_silicoclasticGrain - ( rho_silicoclasticGrain - rho_water) * sedimento_surface_porosity_factor * np.exp(-1.e0*sedimento_compaction_coefficient*z)
        self.rho = rho
        #print(rho[0])
        averageDensity       = np.sum(rho)/len(self.z)
        pseudolayer          = self.getEquivalentSilicoClasticLayer(averageDensity*self.thickness)
        uncompactedthickness = pseudolayer.getUncompactedThickness()
        return uncompactedthickness

    def getEquivalentSilicoClasticLayer(self,mass):
        thicknesses = np.arange(0+dz,20000+dz,dz)
        pseudolayer = sedlayer(sedtype='silicoclastic',thickness=0)
        for thickness in thicknesses:
            pseudolayer = sedlayer(sedtype='silicoclastic',thickness=thickness)
            if (mass<=pseudolayer.averageDensity*thickness):
                break
        return pseudolayer

if __name__ == '__main__':
    print(time.strftime("%a, %d %b %Y %H:%M:%S", time.localtime()))
    t1 = time.time()
    main()
    dt = time.time()-t1
    m, s = divmod(dt, 60)
    h, m = divmod(m, 60)
    print("Process completed in %d hours %02d minutes and %02d seconds" % (h, m, s))

