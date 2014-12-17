from pylab import *
from pymiecoated import Mie
import numpy as np
from scipy.special import erf
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt
from scipy.stats import lognorm


def ddlogn(D,Dm,sigma,N ):
#e.g normalised log normal mode (uncomment or add mode) 
 dND = np.zeros(len(D)) 
 for n in range(len(Dm)):
   print range(len(Dm))
   A =  - (np.log(D) -  np.log(Dm[n]) )**2 / (2 * (np.log([sigma[n]]))**2)
   B = N[n] / ( (2* math.pi)**0.5 * D * np.log([sigma[n]]))
   dND = dND  +   B * np.exp(A)


 return dND

# user defined parameter
#define the number (or mass) distribution (e.g from Kok et al). 
#diameter micron
D = np.arange(0.01,65, 0.05)
cN = 0.9539
Ds=3.4
sigs=3.0
lamb=12
A= np.log(D/Ds)/ (np.sqrt(2) * np.log(sigs))
#care dND/dD or dND/d(lnD) !!
dND =(1/D)* (1/(cN*D*D)) * (1 + erf(A)) * np.exp(-1.*(D/lamb)**3)
#dND = (D/12.62) * (1 + erf(A)) * np.exp(-1.*(D/lamb)**3)
 

# here use a 3 log distribution 
# e.g. alfaro Gomez
if(1==2):
 sig = [1.75, 1.75, 1.75]
 Dm = [0.64, 3.45, 8.67]
 Nf = [0.74, 0.20, 0.06]
 dND = ddlogn(D,Dm,sig,Nf )





#define the 12 new size bin cut of diam from LISA optimised distribution
#you can take also one big bin encompassing the whole distrib (that 's what we do for BC/OC, in this case it is equivalent to integrate over the full distrib mode). 

#here is LISA optimal distrib for dust
sbin =  np.array([0.09,0.18, 0.6, 1.55, 2.50, 3.75, 4.70, 5.70, 7.50, 14.50, 26.0, 41.0, 63.0])
#sbin = [0.01 ,1.,2.5,5.,20. ]
#sbin =[0.98,1.2,2.5,5.,20.]
# set aerosol density in g/m3
rhop = 2650E3 

# define the refractive index, vary spectrally acoording to Wagner et al., ACP 2012 
# I added  an extrapolation point for wavlenght belo 0.305 
#real part doe snot change
# imaginary part 
# define the refractive index from biblio:
# wbib are the waveleght in micrometer
# kbib second the value of k taken form bib 
# e.gfor dust  vary spectrally acoording to Wagner et al., ACP 2012 
# added  an extrapolation point for wavlenght belo 0.305 and above 0.9 
# should complete this data with IR indices 

ndbib= [1.53] * 20 

wbib= [0.2, 0.305, 0.355, 0.405, 0.455, 0.505, 0.555, 0.605, 0.655, 0.705, 0.755, 0.805, 0.855, 0.905, 0.955,\
       2,3, 5,8,13]
kbib = [0.04, 0.027866667,	0.020111111,	0.014933333,	0.010222222,	0.007955556,	0.005133333,	0.003788889,	0.003177778,	0.003033333,	0.003011111,	0.002911111,	0.003111111,	0.003111111,	0.003066667 , 0.0028,0.0025,0.002, 0.00195, 0.0019]



#kbib = [0.04, 0.027866667,     0.020111111,    0.014933333,    0.001,    0.001,    0.001,    0.001,    0.003177778,    0.003033333,    0.003011111,    0.002911111,    0.003111111,    0.003111111,    0.003066667 , 0.0028,0.0025,0.002, 0.00195, 0.0019]

#with this function we can simply interpolate kbib for every value of wavelenght !
kint = interp1d(wbib,kbib,kind='linear')
ndint = interp1d(wbib,ndbib,kind='linear')

#########
# define the spectral band ( wavmin and wavemax) of the radiation model
# here from radiation RegCM standard regcm code  

wavmin = [0.2000 , 0.2450 , 0.2650 , 0.2750 , 0.2850 , \
              0.2950 , 0.3050 , 0.3500 , 0.6400 , 0.7000 ,\
              0.7010 , 0.7010 , 0.7010 , 0.7010 , 0.7020 ,\
              0.7020 , 2.6300 , 4.1600 , 4.1600 ]

wavmax = [0.2450 , 0.2650 , 0.2750 , 0.2850 , 0.2950 ,\
              0.3050 , 0.3500 , 0.6400 , 0.7000 , 5.0000 ,\
              5.0000 , 5.0000 , 5.0000 , 5.0000 , 5.0000 ,\
              5.0000 , 2.8600 , 4.5500 , 4.5500]
#num = 7 for visible band standard scheme
num = 7
# si RRTM set 1 == 1 to overwrite
if(1==1):
    # RRTM Shortwave spectral band limits (wavenumbers)
    #wavenum are in cm-1 / take the inverse for wavelenght 
    wavenum1 =  np.array([2600., 3250., 4000., 4650., 5150., 6150., 7700., \
                          8050.,12850.,16000.,22650.,29000.,38000.,  820.])
    wavenum2 =  np.array([3250., 4000., 4650., 5150., 6150., 7700., 8050., \
                         12850.,16000.,22650.,29000.,38000.,50000., 2600.])

# Longwave spectral band limits (wavenumbers)
    if(1==1): 
       wavenum1 =  np.array([ 10., 350., 500., 630., 700. , 820. , \
                      980. ,1080. ,1180. ,1390. ,1480. ,1800. , \
                     2080. ,2250. ,2380. ,2600. ])
       wavenum2 =  np.array([350. , 500. , 630. , 700. , 820. , 980. , \
                     1080. ,1180. ,1390. ,1480. ,1800. ,2080. , \
                     2250. ,2380. ,2600. ,3250. ])
# determined from indica marc for the longwave           
       wbib=  [3, 4, 5, 6, 7, 8, 9, 10, 12, 14, 16, 20, 30, 1100]
       kbib = [7.54E-2, 5.14E-2, 5.94E-2, 1.41E-1, 1.46E-1, 2.18E-1, 5.35E-1, 3.70E-1, 1.18E-1,  2.28E-1, 2.79E-1, 6.12E-1, 4.95E-1, 6.50E-1 ]
       nbib =   [1.49,  1.51,  1.58, 1.45,  1.46, 1.19,  1.86, 1.84, 1.78,  1.65,  1.55,  2.25, 2.40, 2. ]
       kint = interp1d(wbib,kbib,kind='linear')
       ndint = interp1d(wbib,nbib,kind='linear')
    #convert to wavelenght in micron
    # visible is index 
    wavmin  = np.power(wavenum2*1.E-4,-1.)
    wavmax  = np.power(wavenum1*1.E-4,-1.)

    #num = 7 for visible band standard scheme
    num = 9 # for vis RRTM



print wavmin
print wavmax

# end of user deined parameter
##################################################################
####################################################################3



# 1) calculate the effective radius of each bin (used in regcm )
Sv= np.zeros(len(sbin)-1)
Ss= np.zeros(len(sbin)-1)
normb= np.zeros(len(sbin)-1)

for b in range(len(sbin)-1):
 for dd in range(len(D)):
  if (D[dd] >= sbin[b] and D[dd] <sbin[b+1]):
   Sv[b] = Sv[b] + D[dd]**3 * dND[dd]
   Ss[b] = Ss[b] + D[dd]**2 * dND[dd]
   normb[b] = normb[b] +   dND[dd]
Deff = Sv/Ss
print Deff

# visu

fig = plt.figure()
ax = fig.add_subplot(1,3,1)
#here plot dN/d(logD) !!
ax.plot(D,dND *D , color='blue')
ax.set_yscale('log')
ax.set_xscale('log')
ax.set_xlim([0.08, 65])
ax.set_ylim([1.E-4, 1])
ax.set_title("distribution and bins, dot= Deff")
ax.set_xlabel('D in micrometer')
#plot also bin end effrad 
for x in sbin:
 ax.plot([x, x],[1.E-4, 1], color='black')
for x in Deff:
 ax.scatter(x,2E-4, color='green')


#visu also the ref index spectral variation interp + values 
ax2 = fig.add_subplot(1,3,2)
#xx = np.linspace(0.2,13,500)
xx = np.linspace(3,40,500)
ax2.plot(xx, ndint(xx), color='green')
#ax2.scatter(wbib,kbib)
ax2.set_title(" im ref index, dot = data, line = interp used for oppt calc. ")
ax2.set_xlabel('wavelenght in micrometer')



#2  mie calculation 
######################################
#calculate optical properties per bin
#######################################3
# methode 1 considering only bin effective diam calculated before
if(1==1):
    extb =  np.zeros((len(sbin)-1,len(wavmin)))
    ssab = np.zeros((len(sbin)-1,len(wavmin)))
    asyb =  np.zeros((len(sbin)-1,len(wavmin)))
    for nband in range(len(wavmin)):
        specbin = np.linspace(wavmin[nband],wavmax[nband],50)    
        for db in range(len(Deff)):
          qext = 0.
          ssa = 0.
          asy =0.
          
          for wl in range(len(specbin)):  
            sp = np.pi * Deff[db] / specbin[wl]
            k = kint(specbin[wl])
            nd = ndint(specbin[wl])
            mie = Mie(x=sp,m=complex(nd,k))
            qext= qext + mie.qsca() +  mie.qabs()
            qabs= mie.qabs()
            ssa=  ssa + mie.qsca() / (mie.qsca() +  mie.qabs())
            asy =  asy +  mie.asy()
          extb[db,nband] = extb[db,nband] +  (qext / len(specbin)) / (2./3. * rhop * Deff[db]*1E-6)
    #cross section in m2/g
          ssab [db,nband] = ssab[db,nband] +  (ssa / len(specbin))
          asyb [db,nband] = asyb[db,nband] +  (asy / len(specbin))



if (1==2):
#method 2(slower) double integration  /more precise prblem mie return nan for extrem coarse particles
# calculate the mie parameter for every diameter of the range, averaged on spectral band

  extb =  np.zeros((len(sbin)-1,len(wavmin)))
  ssab = np.zeros((len(sbin)-1,len(wavmin)))
  asyb =  np.zeros((len(sbin)-1,len(wavmin)))


  for nband in range(len(wavmin)):
         specbin = np.linspace(wavmin[nband],wavmax[nband],10)   
         extd = np.zeros(len(D))
         ssad  = np.zeros(len(D))
         asyd  = np.zeros(len(D))
         for db in range(len(D)):
           qext = 0.
           ssa = 0.
           asy =0.
           for wl in range(len(specbin)):  
             sp = np.pi * D[db] / specbin[wl]
             k = kint(specbin[wl])
             nd = ndint(specbin[wl])
             mie = Mie(x=sp,m=complex(nd,k))
             qext= qext + mie.qsca() +  mie.qabs()
             ssa=  ssa + mie.qsca() / (mie.qsca() +  mie.qabs())
             asy=  asy + mie.asy()
           extd[db] = qext / len(specbin) / (2./3. * rhop * D[db]*1E-6)  
           ssad [db] = ssa / len(specbin)
           asyd [db] = asy / len(specbin)
        # perform the bin wighting av according to distibution
           for b in range(len(sbin)-1):  
             if (D[db] >= sbin[b] and D[db] <sbin[b+1]): 
               extb[b,nband] = extb[b,nband] +  extd[db] * dND[db] / normb[b]
               ssab [b,nband] = ssab[b,nband]  + ssad [db]  * dND[db] / normb[b]  
               asyb [b,nband] = asyb[b,nband] + asyd [db]  * dND[db] / normb[b]
 
# 3  nan out for the big bin which can have nan values .... 
# the kext is set  very low which means that the bin won't matter much in the rad 
extb[np.isnan(extb)] = 1.E-20
asyb[np.isnan(asyb)] = 0.99
ssab[np.isnan(ssab)] = 0.5


print "extb vis", extb[:,num]
print "ssab vis", ssab[:,num]
print "asyb vis", asyb[:,num]
print "absb vis", extb[:,num] * ssab[:,num]


#4 visu  oppt / bin

ax3 = fig.add_subplot(1,3,3)
xx = np.linspace(0.2,4,500)
ax3.set_xlim([0.08,65 ])
ax3.set_ylim([0, 5])
for n in range(len(sbin)-1):
 ax3.plot([sbin[n], sbin[n+1]],[extb[n,num],extb[n,num]], color='black')
 ax3.plot([sbin[n], sbin[n+1]],[ssab[n,num],ssab[n,num]], color='blue')
 ax3.plot([sbin[n], sbin[n+1]],[asyb[n,num],asyb[n,num]], color='red')

ax3.set_title(" bin oppt vis: blk:kext, blue:ssa, red:asy ")
ax3.set_xscale('log')
ax3.set_xlabel('D in micrometer')



# finally  save in afile in a format close to mod_rad_aerosol block data
# will just need copy paste + a few edit for fortan 
# file = open("aeroppt_blocl.txt", "w")
np.savetxt('Deff.txt', Deff, fmt = '%7.5E',delimiter = ', ')
np.savetxt('extb.txt', extb.T, fmt = '%7.5E',delimiter = ', ')
np.savetxt('ssab.txt', ssab.T, fmt = '%7.5E',delimiter = ', ')
np.savetxt('asyb.txt', asyb.T, fmt = '%7.5E',delimiter = ', ')
np.savetxt('absb.txt', (extb*ssab).T, fmt = '%7.5E',delimiter = ', ')

plt.show()
exit()
