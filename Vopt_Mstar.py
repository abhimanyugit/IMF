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

def interpolate(stelmass,SM,HM):
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
           
    def getHM_critical_overdensity(self,xmstel):
    #takes xmstel in h^-2 Msolar units and outputs in h^-1 Msolar units
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
               SMred.append(float(SM))      #in Msolar
               HMred.append(float(HM))      #in h^-1 Msolar
               SM_hinv2Msolarred.append(float(SM) - 2 * log10(h))       #in h^-2 Msolar

        with open("SMHM_blue.txt",'r') as bluefile:
            next(bluefile) #to skip headings
            reader = csv.reader(bluefile,delimiter='\t')
            for SM,HM in reader:
                SMblue.append(float(SM))    #in Msolar
                HMblue.append(float(HM))    #in h^-1 Msolar
                SM_hinv2Msolarblue.append(float(SM) - 2 * log10(h))     #in h^-2 Msolar

        if(self.galtype=="LT"):
            SM=list(SM_hinv2Msolarblue)
            halom = interpolate(xmstel,SM_hinv2Msolarblue,HMblue)   #gives HM(h^-1 Msol) for xmstel(h^-2 Msol)
        elif(self.galtype=="ET"):
            SM=list(SM_hinv2Msolarred)
            halom = interpolate(xmstel,SM_hinv2Msolarred,HMred) 

        Halo_m = []
        Halo_c = []
        with open("Conversion_between_M200m_M200c.txt",'r') as file:
            next(file)
            reader = csv.reader(file,delimiter=' ')
            for h_m,h_c in reader:
                Halo_m.append(float(h_m))
                Halo_c.append(float(h_c))

        return interpolate(halom,Halo_m,Halo_c)                     #returns HM as critical overdensity
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
    
    xstelmass=np.arange(10.63,12.02,0.01);
    xstelmass2=np.arange(10.63,12.02,0.01)
    output=""
    output+= "SMred(h^-2Ms)" + "\t" + "HMred(h^-1Ms)" + "\t" + "SMblue(h^-2Ms)" + "\t" + "HMblue(h^-1Ms)" + "Halo masses are w.r.t critical density of universe" +"\n"
    f=open("SMHM.txt",'w')
    f.write(output)
    for i in range (0,139):
        galaxy_list1=Galaxy("ET",xstelmass[i]);
        galaxy_list2=Galaxy("LT",xstelmass2[i]);
        #sigma_list=galaxy_list.getFJsigma();
        #print xstelmass2[i],"\t",galaxy_list2.getHM_critical_overdensity(xstelmass2[i])
        output= str(xstelmass[i])+"\t"+str(galaxy_list1.getHM_critical_overdensity(xstelmass[i]))+"\t"+str(xstelmass2[i])+"\t"+str(galaxy_list2.getHM_critical_overdensity(xstelmass2[i]))+"\n"
        f.write(output)
    f.close()
    #writes SM and HM values for red and blue galaxies. SM in h^-2 Msolar, HM in h^-1 Msolar
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
            SM_TF.append(float(SM) - 2 *log10(h))         #in h^-2 Msolar
            Halo_m=[]
            Halo_c=[]
            with open("Conversion_between_M200m_M200c.txt",'r') as file2:
                next(file2)
                reader = csv.reader(file2,delimiter=' ')
                for h_m,h_c in reader:
                    Halo_m.append(float(h_m))
                    Halo_c.append(float(h_c))
            HM_TF.append(interpolate(float(HM),Halo_m,Halo_c))         #in h^-1 Msolar critical overdensity

    #print SM_TF
    #print HM_TF

    G = 4.301 * 10**(-6)  # in km^2 s^-2 kpc Msolar^-1
    
    output2="SM(hinv2 Msol)"+"\t"+"V200^2"+"\t"+"V200"+"\n"
    
    for i in range(0,139):
        sm=10.59 + (i * 0.01)           #sm in Hinv2 Msolar
        #Eq 2, D10
        log10_V200 = (log10(G) + interpolate(sm,SM_TF,HM_TF)) / 3       #HM returned is w.r.t critical overdensity 
        V200sq.append ( 10 ** ((2 * log10_V200)))
        print sm, V200sq[i], sqrt(V200sq[i]),log10_V200
        output2+= str(sm) + "\t" + str(V200sq[i]) + "\t" + str(sqrt(V200sq[i])) + "\t" + str(log10_V200) + "\n"
        #print sm, 10**(2*log10_V200), interpolate(sm,SM_TF,HM_TF)
    #print output2

    print len(V200sq)    
    f=open("V200sqUpdate.txt",'w')
    f.write(output2)
    f.close()

    Vhalosq=[] 
    Vbulgesq=[]
    Vdisksq=[]
    Vgassq=[]

    with open("Vsq.txt",'r') as file:
        next(file)
        reader = csv.reader(file,delimiter='\t')
        for xm,vbsq,vdsq,vgsq,vhsq in reader:
            Vhalosq.append(float(vhsq))
            Vdisksq.append(float(vdsq))
            Vbulgesq.append(float(vbsq))
            Vgassq.append(float(vgsq))

    #print Vgassq


    #slicing for plotting purposes
    Vhalosq = Vhalosq[16:]
    Vbulgesq = Vbulgesq[16:]
    Vdisksq = Vdisksq[16:]
    Vgassq = Vgassq[16:]
    V200sq= V200sq[1:43:10]
    Vsumsq = []
    Vsumsq2 =[]
    plotmass = np.arange(10.6,11.0,0.1)
    print V200sq
    print Vgassq
    print Vbulgesq
    print Vdisksq
    print Vhalosq

    for i in range (0,5):
        Vsumsq.append(0.5*np.log10(Vbulgesq[i]+Vdisksq[i]+Vgassq[i]+Vhalosq[i]))
        Vsumsq2.append(0.5*np.log10(Vbulgesq[i]+Vdisksq[i]+Vgassq[i]+V200sq[i]))
    print   (Vsumsq)
           

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
    ax2.plot(plotmass,Vsumsq2,label="model with V200")
    #ax2.plot(sm,log10_V200,label="V200")
    #ax.plot(galaxy_list.xmstel,galaxy_list.getGallazi(),label="Gallazi");
    legend(fontsize=6,ncol=2);
    tight_layout();
    savefig("TF.pdf")
    
