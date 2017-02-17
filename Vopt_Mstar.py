#!/usr/bin/env python

# This code is meant to reproduce at least the late type relations presented in
# Dutton et al. 2011

#from scipy.special import iv,kv;
import numpy as np;
from math import *;
from pylab import *;
import csv
from matplotlib.ticker import MultipleLocator

def interpolateSMHM(stelmass,SM,HM):
    size = len(SM)
    for i in range (0,size):
        if((stelmass)>=SM[i] and (stelmass)<SM[i+1]):
            slope = (HM[i+1]-HM[i])/(SM[i+1]-SM[i])             
            halomass = HM[i] + slope * (stelmass-SM[i])
        elif((stelmass<SM[0])or(stelmass>SM[size-1])):
            return -1
    return halomass

class Galaxy(object):
    def __init__(self,galtype,xmstel):
        self.galtype=galtype;
        self.xmstel=xmstel;
        #self.HaloMandelb=HaloMandelb
        
        # Equation 20, D11
        self.xmstel_bell=10.130+0.922*(xmstel-10.);

        if(self.galtype=="ET"):
            # Constants for <M200>/Mstar (Mstar) relation D10
            self.Ms_M_alpha=-0.15;
            self.Ms_M_beta=0.85;
            self.Ms_M_gamma=2.0;
            self.Ms_M_x0=10.8;
            self.Ms_M_y0=1.97;

            # Constants for <R50>/kpc (Mstar) relation D10
            self.R50_Ms_alpha=0.27;
            self.R50_Ms_beta=0.59;
            self.R50_Ms_gamma=1.5;
            self.R50_Ms_x0=9.97;
            self.R50_Ms_y0=0.32;

            # Constants for fdev 
            self.fdev0=0.67
            self.fdev_alpha=0.19
            self.fdev_beta=-0.22
            self.fdev_gamma=1.75
            self.fdev_M0=10.77

            # Constants for Faber Jackson relation
            self.fj_y0=2.23;
            self.fj_x0=10.9;
            self.fj_alpha=0.37;
            self.fj_beta=-0.19;

        elif (self.galtype=="LT"):
            # Constants for <M200>/Mstar (Mstar) relation D10
            self.Ms_M_alpha=-0.50;
            self.Ms_M_beta=0.0;
            self.Ms_M_gamma=1.0;
            self.Ms_M_x0=10.4;
            self.Ms_M_y0=1.61;

            # Constants for <R50>/kpc (Mstar) relation D10
            self.R50_Ms_alpha=0.20;
            self.R50_Ms_beta=0.46;
            self.R50_Ms_gamma=1.95;
            self.R50_Ms_x0=10.39;
            self.R50_Ms_y0=0.75;

            # Constants for fdev 
            self.fdev2=0.34
            self.fdev1=0.08
            self.fdev_gamma=2.5
            self.fdev_M0=10.67

    def getMhalo_chalo(self):
        # Eq 3 D10
        xbyx0=10.0**(self.xmstel_bell-self.Ms_M_x0);
        avM200=self.xmstel_bell+self.Ms_M_y0+self.Ms_M_alpha*log10(xbyx0)+(self.Ms_M_beta-self.Ms_M_alpha)/self.Ms_M_gamma*log10(0.5+0.5*xbyx0**self.Ms_M_gamma);
        # Eq 40 D11
        c200=10.0**(0.830-0.098*(avM200-12.0));
        return avM200,c200;

    def getR50(self):
        # Eq 21 D11
        xbyx0=10.0**(self.xmstel-self.R50_Ms_x0);
        R50=self.xmstel+self.R50_Ms_y0+self.R50_Ms_alpha*log10(xbyx0)+(self.R50_Ms_beta-self.R50_Ms_alpha)/self.R50_Ms_gamma*log10(0.5+0.5*xbyx0**self.R50_Ms_gamma);
        return R50;

    # This requires us to know R_d,R 
    def getRdstar_corr(self):
        if(self.galtype=="LT"):
        # Eq 22 D11
            return 1./(1.+2.5*0.683*(1-0.87));
        else:
            return 1.;

    def getfracdev(self):
        xbyx0=10.0**(self.xmstel-self.fdev_M0);
        if(self.galtype=="LT"):
            return self.fdev2+(self.fdev1-self.fdev2)/(1.+xbyx0**self.fdev_gamma);
        elif(self.galtype=="ET"):
            return self.fdev0*10.0**(self.fdev_alpha*log10(xbyx0)+(self.fdev_beta-self.fdev_alpha)/self.fdev_gamma*log10(0.5+0.5*(xbyx0**self.fdev_gamma)));

    def getTFvelocity(self):
        if(self.galtype=="LT"):
            return 10.0**(2.179+0.259*(self.xmstel-10.3));
        else:
            print "Tully fisher relation should not be used for early types"
            return 0;

    def getFJsigma(self):
        if(self.galtype=="ET"):
            return 10.0**(self.fj_y0+self.fj_alpha*(self.xmstel-self.fj_x0)+self.fj_beta*log10(0.5+0.5*(10**self.xmstel/10**self.fj_x0)))
        else:
            print "Faber Jackson relation should not be used for late types"

    def getVc_re(self):
        if(self.galtype=="LT"):
            return getTFvelocity(self);
        else:
            return 1.54*getFJsigma(self);

    def getMgas(self):
       if(self.galtype=="LT"):
           return self.xmstel-0.27-0.47*(self.xmstel-10.)
       else:
           return 0.0;

    def getGallazi(self):
        if(self.galtype=="ET"):
            return 2.051+0.286*(self.xmstel-10.3)
        else:
            "Nope"
   # This requires you to know the Rd,R
    def getRdHI(self):
       if(self.galtype=="LT"):
           return 0.19;
       else:
           return 0;
           
    def getSMHM(self,xmstel):

       SMred=[]
       HMred=[]
       SMblue=[]
       HMblue=[]
       with open("SMHM_red.txt",'r') as redfile:
           next(redfile) #to skip headings
           reader = csv.reader(redfile,delimiter='\t')
           for SM,HM in reader:
               SMred.append(float(SM))
               HMred.append(float(HM))
               #print SMred 
               #print HMred
       
       with open("SMHM_blue.txt",'r') as bluefile:
		 next(bluefile) #to skip headings
		 reader = csv.reader(bluefile,delimiter='\t')
		 for SM,HM in reader:
		     SMblue.append(float(SM))
		     HMblue.append(float(HM))
		#print SMblue
		#print HMblue	
       if(self.galtype=="LT"):
           return interpolateSMHM(xmstel,SMblue,HMblue)
       elif(self.galtype=="ET"):
		return interpolateSMHM(xmstel,SMred,HMred) 
#----------------------------------------------------------------------------------
           
#----------------------------------------------------------------------------------    

    
    
           
def set_locators(ax,xmaj,xmin,ymaj=0,ymin=0):
    ax.xaxis.set_major_locator(MultipleLocator(xmaj));
    ax.xaxis.set_minor_locator(MultipleLocator(xmin));
    if(ymaj!=0):
        ax.yaxis.set_major_locator(MultipleLocator(ymaj));
    if(ymin!=0):
        ax.yaxis.set_minor_locator(MultipleLocator(ymin));
           
if __name__ == "__main__":
    xstelmass=np.arange(10.28,11.68,0.01);
    #xstelmass=[10.99,11.1];
    #galtype=["ET"];
    xstelmass2=np.arange(10.28,11.68,0.01)
    output=""
    f=open("SMHM.txt",'w')
    for i in range (0,141):
        galaxy_list=Galaxy("ET",xstelmass[i]);
        galaxy_list2=Galaxy("LT",xstelmass[i]);
    #sigma_list=galaxy_list.getFJsigma();
        print xstelmass2[i],"\t",galaxy_list2.getSMHM(xstelmass[i])
        output= str(xstelmass[i])+"\t"+str(galaxy_list.getSMHM(xstelmass[i]))+"\t"+str(xstelmass2[i])+"\t"+str(galaxy_list2.getSMHM(xstelmass2[i]))+"\n"
        f.write(output)
    f.close()
