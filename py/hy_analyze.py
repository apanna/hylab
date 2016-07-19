# -*- coding: utf-8 -*-
"""
"""
__author__          =   "Stephan Schlamminger and Alireza Panna"
__email__           =   "apanna1236@gmail.com"
__status__          =   "Stable"
__date__            =   "06/26/13"
__version__         =   "0.1"

import pylab as py
import numpy as np
import hyreader 
import os

os.chdir('..')
bp = os.getcwd()
if not os.path.exists('Figures'):
    os.makedirs('Figures')
    
rp = "1021 Steel Annealed/400mHz_10mVpp-1Vpp"
rpn = "1021 Steel Not Annealed/400mHz_10mVpp-1Vpp"
fin = "Hysteresis_Integral_10mVpp-1Vpp (0.4Hz)"

suf = '.xlsx'
fd = bp + '/Figures'
HyAn = []
HyNan = []
ax = []
fig = []

myxlabels = ['H (A/m)',   \
             'Hmax (A/m)',\
             'Hmax (A/m)',\
             'Hmax (A/m)',\
             'H (A/m)',   \
             'H (A/m)',   \
             'Bmax (T)',  \
             'H (A/m)',   \
             'H (A/m)',   \
             'H (A/m)',   \
             'Hmax (A/m)',\
             'Hmax (A/m)',\
             'H (A/m)',   \
             'Havg(A/m)',
             ]
           
myylabels = ['B (T)',       \
             'Temp (C)',    \
             'Hcoerc (A/m)',\
             'Breten (T)',  \
             'mu_r',        \
             'mu_diff',     \
             'Br/Bm',       \
             'B (T)',       \
             'M (A/m)',     \
             'B (T)',       \
             'Hcoerc (A/m)',\
             'Breten (T)',  \
             'mu',          \
             'mu_diff',
             ]
           
mytitles = ['MAIN CURVE (Bmax vs. Hmax)',               \
            'Temperature Vs. Hmax',                     \
            'Coercivity Vs. Hmax',                      \
            'Retentivity Vs. Hmax',                     \
            'mu_r vs. H (Main Curve)',                  \
            'mu_diff Vs. H (Main Curve)',               \
            'Retentivity ratio vs. Bmax',               \
            'B Vs. H (Sample with main curve)',         \
            'M Vs. H (Sample with main curve)',         \
            'Mean of the main curve',                   \
            'Mean of the coercivity curve',             \
            'Mean of the retentivity curve',            \
            'Mean of the permeability curve',           \
            'Mean of the differential permeability curve',
            ]
          
tit = fin.split('_')
ann = tit[-1]

for i in range(5):
    fn = fin + '_{0}'.format(i)+suf    
    print fn
    HyAn.append(hyreader.hyreader(bp, rp, fn))
    HyNan.append(hyreader.hyreader(bp, rpn, fn))
    
for i in range(15):
    fig.append(py.figure(i))
    ax.append(fig[i].add_subplot(111))
    py.grid()
    py.suptitle(ann)

for b,c in zip(HyNan,HyAn):
    # get all mean curves
    Hmax, Bmax = c.getallMax()           
    Hmaxn, Bmaxn = b.getallMax()
    # get all temperature data
    temp = c.getallTemp()                 
    tempn = b.getallTemp()
    # get coercivity curves
    Hc = c.getallCoercivity()             
    Hcn = b.getallCoercivity()
    # get retentivity curves
    Br = c.getallRetentivity()            
    Brn = b.getallRetentivity()
    # get mean permeability of mean curve
    muM, muHM = c.getMumean()              
    muMn, muHMn = b.getMumean()
    # get differential permeability of mean curve
    hdiff, mudiff = c.getallMudiffMean()   
    hdiffn, mudiffn = b.getallMudiffMean()
    # get retentivity ratio for all curves
    brratio, bm = c.getBrRatio()           
    brration, bmn = b.getBrRatio()
    
    ax[0].plot(Hmax, Bmax, 'ro-')
    ax[0].plot(Hmaxn, Bmaxn, 'bo-')
    
    ax[1].plot(Hmax, temp, 'ro-')
    ax[1].plot(Hmaxn, tempn, 'bo-')
    
    ax[2].plot(Hmax, Hc, 'ro-')
    ax[2].plot(Hmaxn, Hcn, 'bo-')
    
    ax[3].plot(Hmax, Br, 'ro-') 
    ax[3].plot(Hmaxn, Brn, 'bo-')
    
    ax[4].plot(muHM, muM, 'ro-')
    ax[4].plot(muHMn, muMn, 'bo-')
    
    ax[5].plot(hdiff, mudiff, 'ro-')
    ax[5].plot(hdiffn, mudiffn, 'bo-')
    
    ax[6].plot(bm, brratio, 'ro-')
    ax[6].plot(bmn, brration, 'bo-')

    muMinvn = []
    for i in (muMn):
        muMinvn.append(1/i)
        
    ax[14].plot(tempn, muMinvn, 'bo-')

for i, c in enumerate(HyAn):
    
    if i==0:
        Hmm = np.array(Hmax)
        Bmm = np.array(Bmax)
        Hcm = np.array(Hc)
        Brm = np.array(Br)
        muMm = np.array(muM)
        muHMm = np.array(muHM)
        hdiffm = np.array(hdiff)
        mudiffm = np.array(mudiff)
        
        Hmmn = np.array(Hmaxn)
        Bmmn = np.array(Bmaxn)
        Hcmn = np.array(Hcn)
        Brmn = np.array(Brn)
        muMmn = np.array(muMn)
        muHMmn = np.array(muHMn)
        hdiffmn = np.array(hdiffn)
        mudiffmn = np.array(mudiffn)
        
    else:
        Hmm = np.vstack((Hmm, np.array(Hmax)))
        Bmm = np.vstack((Bmm, np.array(Bmax)))
        Hcm = np.vstack((Hcm, np.array(Hc)))
        Brm = np.vstack((Brm, np.array(Br)))
        muMm = np.vstack((muMm, np.array(muM)))
        muHMm = np.vstack((muHMm, np.array(muHM)))
        hdiffm = np.vstack((hdiffm, np.array(hdiff)))
        mudiffm = np.vstack((mudiffm, np.array(mudiff)))
        
        Hmmn = np.vstack((Hmmn, np.array(Hmaxn)))
        Bmmn = np.vstack((Bmmn, np.array(Bmaxn)))
        Hcmn = np.vstack((Hcmn, np.array(Hcn)))
        Brmn = np.vstack((Brmn, np.array(Brn)))
        muMmn = np.vstack((muMmn, np.array(muMn)))
        muHMmn = np.vstack((muHMmn, np.array(muHMn)))
        hdiffmn = np.vstack((hdiffmn, np.array(hdiffn)))
        mudiffmn = np.vstack((mudiffmn, np.array(mudiffn)))
    
LL = np.shape(Hmm)[0]
norm = np.sqrt(LL)
LLmu = np.shape(muHMm)[0]
normmu = np.sqrt(LLmu)
LLmudiff = np.shape(hdiffm)[0]
normmudiff = np.sqrt(LLmudiff)

ax[9].errorbar(np.mean(Hmm, axis=0), np.mean(Bmm, axis=0),\
            xerr=np.std(Hmm, axis=0)/norm,\
            yerr=np.std(Bmm, axis=0)/norm, fmt='r-')    
ax[10].errorbar(np.mean(Hmm, axis=0), np.mean(Hcm, axis=0),\
            xerr=np.std(Hmm, axis=0)/norm,\
            yerr=np.std(Hcm, axis=0)/norm, fmt='r-')    
ax[11].errorbar(np.mean(Hmm, axis=0),np.mean(Brm, axis=0),\
            xerr=np.std(Hmm, axis=0)/norm,\
            yerr=np.std(Brm, axis=0)/norm, fmt='r-')    
ax[12].errorbar(np.mean(muHMm, axis=0), np.mean(muMm, axis=0),\
            xerr=np.std(muHMm, axis=0)/normmu,\
            yerr=np.std(muMm, axis=0)/normmu, fmt='r-')
ax[13].errorbar(np.mean(hdiffm, axis=0), np.mean(mudiffm, axis=0),\
            xerr=np.std(hdiffm, axis=0)/normmudiff,\
            yerr=np.std(mudiffm, axis=0)/normmudiff, fmt='r-')
                       
LL = np.shape(Hmmn)[0]
norm = np.sqrt(LL)
LLmu = np.shape(muHMmn)[0]
normmu = np.sqrt(LLmu)
LLmudiff = np.shape(hdiffmn)[0]
normmudiff = np.sqrt(LLmudiff)

ax[9].errorbar(np.mean(Hmmn, axis=0), np.mean(Bmmn, axis=0),\
            xerr=np.std(Hmmn, axis=0)/norm,\
            yerr=np.std(Bmmn, axis=0)/norm,fmt='b-')    
ax[10].errorbar(np.mean(Hmmn, axis=0), np.mean(Hcmn, axis=0),\
            xerr=np.std(Hmmn, axis=0)/norm,\
            yerr=np.std(Hcmn, axis=0)/norm,fmt='b-')    
ax[11].errorbar(np.mean(Hmmn, axis=0), np.mean(Brmn, axis=0),\
            xerr=np.std(Hmmn, axis=0)/norm,\
            yerr=np.std(Brmn, axis=0)/norm, fmt='b-')    
ax[12].errorbar(np.mean(muHMmn, axis=0), np.mean(muMmn, axis=0),\
            xerr=np.std(muHMmn, axis=0)/normmu,\
            yerr=np.std(muMmn, axis=0)/normmu, fmt='b-')
ax[13].errorbar(np.mean(hdiffmn, axis=0), np.mean(mudiffmn, axis=0),\
            xerr=np.std(hdiffmn, axis=0)/normmudiff,\
            yerr=np.std(mudiffmn, axis=0)/normmudiff, fmt='b-')
            
# Sample Curves 
Htot = []
Btot = []
Ht = []
Mt = []

Htot,Btot = HyAn[-1].getallHB()
Ht,Mt = HyAn[-1].getallHM()

ax[7].plot(Hmax,Bmax, 'k--', linewidth = 3)
for i, j in zip(Htot, Btot):
    ax[7].plot(i, j)
for i, j in zip(Ht, Mt):
    ax[8].plot(i, j)
    
    
for i in range(14):
    ax[i].set_xlabel(myxlabels[i])
    ax[i].set_ylabel(myylabels[i])
    ax[i].set_title(mytitles[i])

fig1 = py.figure(16)
ax1 = fig1.add_subplot(121)
ax2 = fig1.add_subplot(122)

# Energy Product Curve:
for c in HyAn:
     bh2, b2 = c.getEproduct()
     hd, bd = c.getDemagCurve()
     for i, j in zip(hd, bd):
         ax1.plot(i, j)    
     for i, j in zip(bh2, b2):
         ax2.plot(i, j)
         
#fig1.show()
#for i in range(14):
#    fig[i].show()

r = rp.split('_')
for i in range(14):
  sname = mytitles[i] + '_' + r[-1]  
  fig[i].savefig(fd+'/' + sname + '.png')
  
  
  

