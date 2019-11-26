import numpy as np
import os
import scipy
import matplotlib.pylab as plt
from astropy.table import Table
from astropy.coordinates import SkyCoord
import astropy.units as u
from astropy.table import Column
import matplotlib
import scipy.stats as scipystats
from astropy.cosmology import WMAP9 as cosmo
import astropy.io.ascii as ascii
from importlib import reload   # reload if I accidentally assigned plot keywords to something
import matplotlib.gridspec as gridspec 
from matplotlib.colorbar import Colorbar # For dealing with Colorbars the proper way 
from scipy.optimize import curve_fit
from matplotlib import lines
from scipy.stats import norm
import matplotlib.mlab as mlab
from scipy.optimize import leastsq
from matplotlib.ticker import ScalarFormatter

import matplotlib.ticker as ticker

Radians=206264.8065         #number of arcsecs in a radian


def separate(umv, vmj):
    qui = []; sf =[]; dsf = []
    
    for j in range(0,len(umv)):
        x = vmj[j]
        y = umv[j]

        if y > ya and y > s1 * x + yint:
            qui.append(j)
        else:
            sf.append(j)

    return qui, sf



goodsn = Table.read('/Users/kalina/Work/Mass_Size/GOODSN_catalog')
goodss = Table.read('/Users/kalina/Work/Mass_Size/GOODSS_catalog')
cosmos = Table.read('/Users/kalina/Work/Mass_Size/COSMOS_catalog')
uds = Table.read('/Users/kalina/Work/Mass_Size/UDS_catalog')


# USE: Non-overlap completeness instead
def PrelimCuts_CANDELS(cat, completeness90): #, completeness90_reg2):
    # GALAPAGOS CUTS
    ind = range(len(cat['N_GALFIT_BAND']))
    ind_ = [i for i in ind if (cat['FLAG_GALFIT'][i] == 2)]
    ind_ = [i for i in ind_ if (cat['use_phot'][i] == 1)]
    ind_ = [i for i in ind_ if (cat['z_2'][i] > 0.2)]
    ind_ = [i for i in ind_ if (cat['star_flag'][i] != 1)]

    ind_ = [i for i in ind_ if all(cat['N_GALFIT_BAND'][i] > 0.205) and all(cat['N_GALFIT_BAND'][i] < 7.95)]
    ind_ = [i for i in ind_ if all(cat['RE_GALFIT_BAND'][i] > 0.305) and all(cat['RE_GALFIT_BAND'][i] < 395)]
    ind_ = [i for i in ind_ if all(cat['Q_GALFIT_BAND'][i] > 0.001) and all(cat['Q_GALFIT_BAND'][i] <= 1)]
    ind_ = [i for i in ind_ if all(cat['MAG_GALFIT_BAND'][i] > 0) and all(cat['MAG_GALFIT_BAND'][i] < 40)]
    ind_ = [i for i in ind_ if all(cat['MAG_GALFIT_BAND'][i] > cat['MAG_BEST'][i] -5)] 
    ind_ = [i for i in ind_ if all(cat['MAG_GALFIT_BAND'][i] < cat['MAG_BEST'][i] +5)]


    ind_ = [i for i in ind_ if ( -2.5*np.log10(cat['f_f160w'][i])+25 <=completeness90)]


    print('number of sources left:', len(ind_))
    return ind_


goodsn = goodsn[PrelimCuts_CANDELS(goodsn, 24)]
goodss = goodss[PrelimCuts_CANDELS(goodss, 24)]
cosmos = cosmos[PrelimCuts_CANDELS(cosmos, 24)]
uds = uds[PrelimCuts_CANDELS(uds, 24)]



#define UVJ lines
s1=.99
s2=1.43 #slope dusty line
xa=0.854
xb=2.5
yint=0.45
xc=1.55
ya=s1*xa+yint
yb=s1*xb+yint
yc=xc*s1+yint
barx = [0.25, xc]
bary = [0, yc]
bar2x = [-2,xa]
bar2y = [ya,ya]
bar3x = [xa,xb]
bar3y = [ya,yb]

# CANDELS UVJ 
s1=.99
s2=1.43 #slope dusty line
xa=0.854
xb=2.5
yint=0.45
xc=1.55
ya=s1*xa+yint
yb=s1*xb+yint
yc=xc*s1+yint
barx = [0.25, xc]
bary = [0, yc]
bar2x = [-2,xa]
bar2y = [ya,ya]
bar3x = [xa,xb]
bar3y = [ya,yb]




def plot_UVJ(ax, minimum, maximum, colormap, red_upper, red_lower):
    redshift_limit_lower=red_lower; redshift_limit_upper=red_upper
    U_band_HFF=[]; V_band_HFF=[]; J_band_HFF=[]
    z_indexes =[]; sfr_indexes=[]
    
    catalog_mags = [goodsn, goodss, cosmos, uds]
    
    for c in range (len(catalog_mags)): 
    
        sfr_indexes.append([catalog_mags[c]['lssfr'][i] for i in range(len(catalog_mags[c]['lssfr'])) if catalog_mags[c]['lssfr'][i] < -1])
        z_indexes.append([catalog_mags[c]['z_2'][i] for i in range(len(catalog_mags[c]['lssfr'])) if catalog_mags[c]['lssfr'][i] < -1])


        #print(sfr_indexes)
        U_band_HFF.append(np.array([-2.5*np.log10(catalog_mags[c]['l153'][i]) + 25 - catalog_mags[c]['dm'][i] for i in range(len(catalog_mags[c]['lssfr'])) if catalog_mags[c]['lssfr'][i] < -1]))
        #print(U_band_HFF)
        V_band_HFF.append(np.array([-2.5*np.log10(catalog_mags[c]['l155'][i]) + 25 - catalog_mags[c]['dm'][i] for i in range(len(catalog_mags[c]['lssfr'])) if catalog_mags[c]['lssfr'][i] < -1]))
        J_band_HFF.append(np.array([-2.5*np.log10(catalog_mags[c]['l161'][i]) + 25 - catalog_mags[c]['dm'][i] for i in range(len(catalog_mags[c]['lssfr'])) if catalog_mags[c]['lssfr'][i] < -1]))

    U_HFF = [item for sublist in U_band_HFF for item in sublist]
    V_HFF = [item for sublist in V_band_HFF for item in sublist]
    J_HFF = [item for sublist in J_band_HFF for item in sublist]
    # print(np.shape(J_HFF))
    sfr_indexes = [item for sublist in sfr_indexes for item in sublist]
    z_indexes = [item for sublist in z_indexes for item in sublist]

    inds_ = range(len(sfr_indexes))
    indexes= [k for k in inds_ if (z_indexes[k] <= redshift_limit_upper and z_indexes[k] >= redshift_limit_lower)]

    im_uvj = ax.scatter(np.array(V_HFF)[indexes]- np.array(J_HFF)[indexes], np.array(U_HFF)[indexes]-np.array(V_HFF)[indexes], c=np.array(sfr_indexes)[indexes], \
                    vmin=minimum, vmax=maximum,  cmap=colormap, s = 5)

    circ1 = ax.scatter(-10,-10, color ='white',  s = 20, label='%i objects' %(len(indexes)))
    ax.legend(handles=[circ1], markerscale=0.01, markerfirst=False, edgecolor = 'white', fontsize=12,loc=2, borderpad=0.5)

    
    # --------------------------------------------------------
    # THE COLORBAR
    cbax = plt.subplot(gs[:,2]) # Place it where it should be.
    # --------------------------------------------------------
    cb = Colorbar(ax = cbax, mappable=im_uvj, ticklocation = 'right')
    cb.set_label(r'$\log_{10}$ (sSFR $\ [$yr$^{-1}]) $', labelpad=10)
    #plt.savefig('/Users/kalina/Work/Mass_Size/UVJ_CANDELS.pdf')







plt.close(1)
fig = plt.figure(1,figsize=(14.16,8))
gs = gridspec.GridSpec(2,3, height_ratios=[1,1], width_ratios=[1,1,0.05],left=0.15, right=0.85, bottom=0.02, top=0.93, wspace=0.05, hspace=0.05)
ax1 = plt.subplot(gs[0,0]); ax2 = plt.subplot(gs[0,1])
ax3 = plt.subplot(gs[1,0]); ax4 = plt.subplot(gs[1,1])

circ1 = ax1.scatter(-1,-1, color ='white',  s = 20, alpha=0.0)
leg1 = ax1.legend([circ1],['0.2 $\leq$ z $\leq$ 0.5'], loc=4, handletextpad=0.3, borderpad=0.2, ncol=1, frameon=False)
ax1.legend(handles=[circ1], loc=4, handletextpad=0, borderpad=0.4, ncol=1, framealpha=0.0)
ax1.add_artist(leg1)
leg1 = ax2.legend([circ1],['0.5 $\leq$ z $\leq$ 1.0'], loc=4, handletextpad=0.3, borderpad=0.2, ncol=1, frameon=False)
ax2.legend(handles=[circ1], loc=4, handletextpad=0, borderpad=0.4, ncol=1, framealpha=0.0)
ax2.add_artist(leg1)
leg1 = ax3.legend([circ1],['1.0 $\leq$ z $\leq$ 1.5'], loc=4, handletextpad=0.3, borderpad=0.2, ncol=1, frameon=False)
ax3.legend(handles=[circ1], loc=4, handletextpad=0, borderpad=0.4, ncol=1, framealpha=0.0)
ax3.add_artist(leg1)
leg1 = ax4.legend([circ1],['1.5 $\leq$ z $\leq$ 2.0'], loc=4, handletextpad=0.3, borderpad=0.2, ncol=1, frameon=False)
ax4.legend(handles=[circ1], loc=4, handletextpad=0, borderpad=0.4, ncol=1, framealpha=0.0)
ax4.add_artist(leg1)



plot_UVJ(ax=ax1, maximum = -8, minimum=-13, colormap='jet_r', red_lower=0.2, red_upper=0.5)

# # Create a legend for the first line.
# circ1 = ax1.scatter(-1,-1, color ='white',  s = 20, alpha=0.0 , label='0.2 $\leq$ z $\leq$ 0.5')
# first_legend = ax1.legend(handles=[circ1], loc=4,handletextpad=0, borderpad=0.4, ncol=1, framealpha=0.0)

ax1.set_xticklabels([])   
ax1.set_ylabel('U $-$ V'); ax3.set_ylabel('U $-$ V')
ax3.set_xlabel('V $-$ J'); ax4.set_xlabel('V $-$ J')


circ1 = ax1.scatter(-1,-1, color ='white',  s = 20, alpha=0.0)
plot_UVJ(ax=ax2, maximum = -8, minimum=-13, colormap='jet_r', red_lower=0.5, red_upper=1)



ax2.set_yticklabels([])
ax2.set_xticklabels([])


plot_UVJ(ax=ax3, maximum = -8, minimum=-13, colormap='jet_r', red_lower=1, red_upper=1.5)


plot_UVJ(ax=ax4, maximum = -8, minimum=-13, colormap='jet_r', red_lower=1.5, red_upper=2)
ax4.set_yticklabels([])

for ax in [ax1, ax2, ax3, ax4]:
    #ax.plot(barx, bary, color='black', linewidth=2, linestyle='--')
    ax.plot(bar2x, bar2y, color='black', linewidth=2)
    ax.plot(bar3x, bar3y, color='black', linewidth=2)
    ax.set_xlim(left=-0.6,right=2.5)
    ax.set_ylim(bottom=0,top=2.5)
    ax.grid(linestyle='--',linewidth=0.2)
    



plt.show()


