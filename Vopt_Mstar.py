#!/usr/bin/env python

# This code is meant to reproduce at least the late type relations presented in
# Dutton et al. 2011
#-------------------------------------------------------------------------------
"""
Make sure these files are there in the directory before compiling:
SMHM_red.txt        SMHM relation for ET galaxies
SMHM_blue.txt       SMHM relation for LT galaxies
Vsq.txt             Velocity components for LT galaxy

These files are created by this code:
SMHM.txt            SMHM relation data
VsqUpdate.txt       Square of Vhalo obtained from mandelbaum data

Plots 
FJ.pdf
TF.pdf
"""
#-------------------------------------------------------------------------------
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

h=0.673

class Galaxy(object):
    def __init__(self,galtype,xmstel):
        self.galtype=galtype;
        self.xmstel=xmstel;
        
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
        #gives SMHM relation from Mandelbaum 2015
        #SM in Msolar, HM in h^-1 Msolar
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
       
        with open("SMHM_blue.txt",'r') as bluefile:
            next(bluefile) #to skip headings
            reader = csv.reader(bluefile,delimiter='\t')
            for SM,HM in reader:
                SMblue.append(float(SM))
                HMblue.append(float(HM))

        if(self.galtype=="LT"):
            return interpolateSMHM(xmstel,SMblue,HMblue)
        elif(self.galtype=="ET"):
            return interpolateSMHM(xmstel,SMred,HMred) 

    def getHM_critical_overdensity(self,xmstel):
    #converts halo mass from mean overdensity to critical overdensity scale
    #get mean overdensity halo mass from a given stellar mass
        SMred=[]
        HMred=[]
        SMblue=[]
        HMblue=[]
        SM_hinv2Msolarred=[]
        SM_hinv2Msolarblue=[]
        with open("SMHM_red.txt",'r') as redfile:
            next(redfile) #to skip headings
            reader = csv.reader(redfile,delimiter='\t')
            for SM,HM in reader:
               SMred.append(float(SM))
               HMred.append(float(HM))
               SM_hinv2Msolarred.append(float(SM) - 2 * np.log10(h))

        with open("SMHM_blue.txt",'r') as bluefile:
            next(bluefile) #to skip headings
            reader = csv.reader(bluefile,delimiter='\t')
            for SM,HM in reader:
                SMblue.append(float(SM))
                HMblue.append(float(HM))
                SM_hinv2Msolarblue.append(float(SM) - 2 * np.log10(h))

        if(self.galtype=="LT"):
            halom = interpolateSMHM(xmstel,SMblue,HMblue)
        elif(self.galtype=="ET"):
            halom = interpolateSMHM(xmstel,SMred,HMred) 

        Halo_m = []
        Halo_c = []
        with open("Conversion_between_M200m_M200c.txt",'r') as file:
            next(file)
            reader = csv.reader(file,delimiter=' ')
            for h_m,h_c in reader:
                Halo_m.append(float(h_m))
                Halo_c.append(float(h_c))

        return interpolateSMHM(halom,Halo_m,Halo_c)
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
    xstelmass2=np.arange(10.28,11.68,0.01)
    output=""
    output+= "SMred" + "\t" + "HMred" + "\t" + "SMblue" + "\t" + "HMblue" + "\n"
    f=open("SMHM.txt",'w')
    f.write(output)
    for i in range (0,141):
        galaxy_list1=Galaxy("ET",xstelmass[i]);
        galaxy_list2=Galaxy("LT",xstelmass[i]);
        #sigma_list=galaxy_list.getFJsigma();
        print xstelmass2[i],"\t",galaxy_list2.getSMHM(xstelmass[i])
        output= str(xstelmass[i])+"\t"+str(galaxy_list1.getSMHM(xstelmass[i]))+"\t"+str(xstelmass2[i])+"\t"+str(galaxy_list2.getSMHM(xstelmass2[i]))+"\n"
        f.write(output)
    f.close()
    #writes SM and HM values for red and blue galaxies. SM in Msolar, HM in h^-1 Msolar
    #-----------------------------------------------------------------------------------
    
    xms=np.arange(9,12,0.1)                         #an array of Stellar masses
    galaxy_list=Galaxy("ET",xms)                    #a list of objects
    ax=subplot(111)                                 #subplot(abc) creates axb grid and c is the index of the plot
    #ax.set_xlim([9,12])
    #ax.set_ylim([1.5,2.8]);
    set_locators(ax,0.5,0.1,0.2,0.02)
    ax.autoscale(enable=True, axis='both', tight=None)
    ax.set_xlabel("log_10 $M_* (M_\odot)$");
    ax.set_ylabel("$log_{10}$ ($\sigma_{e}$ /[km/s])");
    ax.plot(galaxy_list.xmstel,np.log10(galaxy_list.getFJsigma()),label="FJ");
    #ax.plot(galaxy_list.xmstel,galaxy_list.getGallazi(),label="Gallazi");
    legend(fontsize=6,ncol=2);   
    tight_layout();
    savefig("FJ.pdf")
    close()
    #----------------------------------------------------------------------------------------
    
    #testobj = Galaxy("ET",11)
    #print testobj.getHM_critical_overdensity(11)

    SM_TF=[]
    HM_TF=[]
    V200sq=[]
    
    with open("SMHM_blue.txt",'r') as file:
        next(file)
        rdr = csv.reader(file,delimiter='\t')
        for SM,HM in rdr:
            SM_TF.append(float(SM))         #in Msolar
            HM_TF.append(float(HM))         #in h^-1 Msolar

    G = 4.301 * 10**(-6)  # in km^2 s^-2 kpc Msolar^-1
    
    output2=""
    
    for i in range(0,141):
        sm=10.24 + (i * 0.01)
        #Eq 2, D10
        log10_V200 = (log10(G) + interpolateSMHM(sm,SM_TF,HM_TF)) / 3       #requires some parameter addition to account for scale change of HM
        output2+= str(sm) + "\t" + str(log10_V200) + "\n"
        if((sm*100)%10 == 0 ):
            V200sq.append ( 10 ** ((2 * log10_V200)))
        print sm, 10**log10_V200, interpolateSMHM(sm,SM_TF,HM_TF)
    ##print Vhalosq   #, Vgassq, Vdisksq, Vdisksq

    #slicing for plotting purposes
    Vhalosq = Vhalosq[13:]
    Vbulgesq = Vbulgesq[13:]
    Vdisksq = Vdisksq[13:]
    Vgassq = Vgassq[13:]
    V200sq= V200sq[:8]
    Vsumsq = []
    
    plotmass = np.arange(10.3,11.1,0.1)
    
    for i in range (0,8):
        Vsumsq.append(np.log10(Vbulgesq[i]+Vdisksq[i]+Vgassq[i]+V200sq[i]))#+Vhalosq[i]))
    #print   (Vsumsq)
     
    
    for i in range (0,8):
        Vsumsq[i] = Vsumsq[i] - 2.0                #transposing just to compare slopes
        

    galaxy_list=Galaxy("LT",xms)
    ax2=subplot(111)     #subplot(abc) creates axb grid and c is the index of the plot
    #ax.set_xlim([9,12])
    #ax.set_ylim([1.5,2.8]);
    set_locators(ax,0.5,0.1,0.2,0.02)
    ax2.autoscale(enable=True, axis='both', tight=None)
    ax2.set_xlabel("log_10 $M_* (M_\odot)$");
    ax2.set_ylabel("$log_{10}$ ($V_{2.2}$ /[km/s])");
    ax2.plot(galaxy_list.xmstel,np.log10(galaxy_list.getTFvelocity()),label="TF");
    ax2.plot(plotmass,Vsumsq,label="model");
    ax2.plot(sm,log10_V200,label="V200")
    #ax.plot(galaxy_list.xmstel,galaxy_list.getGallazi(),label="Gallazi");
    legend(fontsize=6,ncol=2);
    tight_layout();
    savefig("TF.pdf")
    
