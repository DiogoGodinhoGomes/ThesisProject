#!/usr/bin/env python

import copy as cp
import glob as gb
import func as fc
import numpy as np
import argparse as ap
import matplotlib.pyplot as pl
from astropy.io import fits
from scipy import interpolate
from scipy import stats


if __name__ == '__main__':
    # --------------------
    # -CORRELATION PLOTS -
    # --------------------
    flux = np.copy(palig[0, :, 10])

    print('\nAIRMASS CORRELATIONS (AM VS. FL):')
    linr_lin = stats.linregress(am, flux)
    print('LIN-LIN: '+str(linr_lin[2]**2.0))
    linr_log = stats.linregress(np.log10(am), np.log10(flux))
    print('LOG-LOG: '+str(linr_log[2]**2.0))
    linr_linlog = stats.linregress(am, np.log10(flux))
    print('LIN-LOG: '+str(linr_linlog[2]**2.0))
    linr_loglin = stats.linregress(np.log10(am), flux)
    print('LOG-LIN: '+str(linr_loglin[2]**2.0))

    print('\nTEMPERATURE CORRELATIONS (TP VS. FL):')
    linr_lin = stats.linregress(tp, flux)
    print('LIN-LIN: '+str(linr_lin[2]**2.0))
    linr_log = stats.linregress(np.log10(tp), np.log10(flux))
    print('LOG-LOG: '+str(linr_log[2]**2.0))
    linr_linlog = stats.linregress(tp, np.log10(flux))
    print('LIN-LOG: '+str(linr_linlog[2]**2.0))
    linr_loglin = stats.linregress(np.log10(tp), flux)
    print('LOG-LIN: '+str(linr_loglin[2]**2.0))

    print('\nSEEING CORRELATIONS (SI VS. FL):')
    linr_lin = stats.linregress(si, flux)
    print('LIN-LIN: '+str(linr_lin[2]**2.0))
    linr_log = stats.linregress(np.log10(si), np.log10(flux))
    print('LOG-LOG: '+str(linr_log[2]**2.0))
    linr_linlog = stats.linregress(si, np.log10(flux))
    print('LIN-LOG: '+str(linr_linlog[2]**2.0))
    linr_loglin = stats.linregress(np.log10(si), flux)
    print('LOG-LIN: '+str(linr_loglin[2]**2.0))

    iname = 'Images/p0_am_fl.png'
    pl.plot(am, flux, 'rx')
    pl.xlabel('Airmass')
    pl.ylabel('Flux')
    pl.grid(True)
    pl.savefig(iname, dpi=500)
    pl.close()

    iname = 'Images/p1_sn_am.png'
    pl.plot(range(len(am)), flux, 'rx')
    pl.xlabel('Spectra number')
    pl.ylabel('Airmass')
    pl.grid(True)
    pl.savefig(iname, dpi=500)
    pl.close()

    iname = 'Images/p2_tp_fl.png'
    pl.plot(tp, flux, 'rx')
    pl.xlabel('Temperature')
    pl.ylabel('Flux')
    pl.grid(True)
    pl.savefig(iname, dpi=500)
    pl.close()

    iname = 'Images/p3_sn_tp.png'
    pl.plot(range(len(tp)), tp, 'rx')
    pl.xlabel('Spectra number')
    pl.ylabel('Temperature')
    pl.grid(True)
    pl.savefig(iname, dpi=500)
    pl.close()

    iname = 'Images/p4_si_fl.png'
    pl.plot(si, flux, 'rx')
    pl.xlabel('Seeing')
    pl.ylabel('Flux')
    pl.grid(True)
    pl.savefig(iname, dpi=500)
    pl.close()

    iname = 'Images/p5_sn_si.png'
    pl.plot(range(len(si)), si, 'rx')
    pl.xlabel('Spectra number')
    pl.ylabel('Seeing')
    pl.grid(True)
    pl.savefig(iname, dpi=500)
    pl.close()

# ----------------------
# - SPECTRA CORRECTION -
# ----------------------
    '''
            if i == 3 and j == 6:
                pl.plot(xold[450:550], yold[450:550], 'kx',
                        xnew, ynew, 'rx')
                pl.title('CHIP %d - Spectrum %d' % (i+1, j))
                pl.xlabel('Pixel number')
                pl.ylabel('Flux')
    pl.axis([500, 525, 0, 35])
    pl.grid(True)

    par = '03'
    iname = 'Data/Temp/ipol_' + par + '.png'
    pl.savefig(iname, dpi=500)
    pl.close()
    '''

# -----------------
# - COMPLEX PLOTS -
# -----------------
chip_num, spec_num = 1, 12
rows, cols = 33, 1024
ratio = float(rows)/float(cols)
bot, top, lef, rig = 0.08, 0.16, 0.06, 0.06
quant = (1.0-(bot+top))/float(spec_num)
fig_w = 18
fig_h = ratio*float(fig_w*spec_num)*((1.0-(lef+rig))/(1.0-(bot+top)))

fig = pl.figure(figsize=(fig_w, fig_h))
fig.suptitle(r'$Chip\,\,%d-Waterfall\,\,plot$' % chip_num, fontsize=24)

gb = fig.add_axes([lef, bot, 1.0-(lef+rig), 1.0-(bot+top)])
gb.axis([0, cols-1, 0, rows-1])
gb.xaxis.set_ticks(np.arange(15, cols, 100))
gb.set_xticklabels(map(str, gb.get_xticks()*0.00022+3.6367), fontsize=8)
gb.set_xlabel(r'$Wavelength-\lambda\,\,(\mu m)$', fontsize=12)

gt = gb.twiny()
gt.axis([0, cols-1, 0, rows-1])
gt.xaxis.set_ticks(np.arange(0, cols, 200))
gt.set_xticklabels(map(str, gt.get_xticks()), fontsize=8)
gt.set_xlabel(r'$Pixel\,\,number$', fontsize=12)

gb.yaxis.set_ticks([])
gb.set_ylabel(r'$Spectrum\,\,number$', fontsize=12, labelpad=20)

sub = []
for i in range(spec_num):
    a = np.random.random((rows, cols))
    sub.append(fig.add_axes([lef, bot+i*quant, 1.0-(lef+rig), quant]))
    sub[i].imshow(a**(i+1), cmap='gray', interpolation='nearest')
    sub[i].xaxis.set_ticks([])
    sub[i].yaxis.set_ticks(np.arange(0, 33, 11, dtype=int))
    sub[i].set_yticklabels(map(str, sub[i].get_yticks()), fontsize=8)
    sub[i].invert_yaxis()

pl.show()
pl.close()

# -------------------
# - WATERFALL PLOTS -
# -------------------
f, axarr = pl.subplots(4, 4, gridspec_kw={'wspace': 0, 'hspace': 0})

for i, ax in enumerate(f.axes):
    ax.set_xticklabels([])
    ax.set_yticklabels([])

pl.show()
pl.close()

fig = pl.figure(figsize=(18, 2))

ax1 = fig.add_subplot(3, 1, 1)
que = pextr[1]
dta = np.std(que)
ax1.imshow(que, cmap='gray').set_clim(1.0 - dta, 1.0 + dta)
ax1.xaxis.set_visible(False)
ax1.yaxis.set_visible(False)

ax2 = fig.add_subplot(3, 1, 2)
que = paraa[1]
dta = np.std(que)
ax2.imshow(que, cmap='gray').set_clim(1.0 - dta, 1.0 + dta)
ax2.xaxis.set_visible(False)
ax2.yaxis.set_visible(False)

ax3 = fig.add_subplot(3, 1, 3)
que = psraa[1]
dta = np.std(que)
ax3.imshow(que, cmap='gray').set_clim(1.0 - dta, 1.0 + dta)
ax3.xaxis.set_visible(False)
ax3.yaxis.set_visible(False)

pl.tight_layout(w_pad=0.0, h_pad=0.0)
pl.show()

# ----------------------
# - STRONG LINES PLOTS -
# ----------------------
color = ['k', 'r', 'g', 'b']
for i in range(int(len(palig)/2)):
    iname = 'Data/Temp/c0' + str(i+1) + 'lines.png'
    a = np.mean(palig[2*i+1], axis=0)
    pl.figure(figsize=(18, 4))
    pl.plot(range(len(a)), a, color[i] + '-', markersize=1)
    pl.axis([0, 1023, 0.8, 1.1])
    pl.subplots_adjust(bottom=0.075, top=0.925, left=0.075, right=0.925)
    pl.grid(True)
    pl.savefig(iname, dpi=500)
    pl.close()

    a = (np.diff(a)**2)
    a = a/np.std(a)
    a[a < 0.05] = 0.045
    a[a >= 0.05] = 0.005

    iname = 'Data/Temp/c0' + str(i+1) + 'diffs.png'
    pl.figure(figsize=(18, 4))
    pl.plot(range(len(a)), a, color[i] + '.', markersize=2)
    pl.axis([0, 1023, 0, 0.05])
    pl.subplots_adjust(bottom=0.075, top=0.925, left=0.075, right=0.925)
    pl.grid(True)
    pl.savefig(iname, dpi=500)
    pl.close()

# ---------------
# - MODEL PLOTS -
# ---------------
    nam = 'Data/Models/spec_1-1-1-1-1-1.dat'
    iname = 'Data/Temp/model.png'
    mod = np.swapaxes(np.array(np.genfromtxt(nam), dtype=None), 0, 1)
    xps = mod[0]
    yps = mod[1]/np.median(mod[1])
    pl.figure(figsize=(18, 4))
    pl.plot(xps, yps, 'k-.', linewidth=0.5, markersize=1)
    yps = psraa[1][0]
    xps = pwave[0][0][yps != 0]*10.0
    yps = yps[yps != 0]
    pl.plot(xps, yps, 'r.', linewidth=0.5, markersize=2)
    pl.axis([pwave[0][0][0]*10.0, pwave[0][0][-1]*10.0, 0.9, 1.1])
    pl.grid(True)
    pl.savefig(iname, dpi=500)
    pl.close()

# ---------------------------
# - FLUX VS WAVELENGTH PLOT -
# ---------------------------
    spc, wla, fla = 6, [], []
    for i in range(len(pwave)):
        for j in range(len(pwave[i])):
            if i == 0:
                wla.append([])
                fla.append([])
            wla[j] = np.insert(wla[j], len(wla[j]), np.copy(pwave[i][j]))
            fla[j] = np.insert(fla[j], len(fla[j]), np.copy(palig[2*i+1][j]))
    xfw = wla[spc][fla[spc] != 0]
    yfw = fla[spc][fla[spc] != 0]
    pl.figure(figsize=(18, 5))
    pl.plot(xfw, yfw, 'rx', markersize=3)
    pl.grid(True)
    pl.savefig('Data/Temp/fvsw_s' + str(spc) + '_n.png', dpi=500)
    pl.close()

# -------------------------
# - WAVELENGTH DIFFS PLOT -
# -------------------------
    w = np.diff(model[0])
    z = np.diff(pcomp[0][0])
    pl.figure(figsize=(18, 5))
    pl.plot(model[0][1:], w, 'rx', markersize=3)
    pl.plot(pcomp[0][0][1:], z, 'kx', markersize=2)
    pl.axis([3.435, 3.56, 1.1e-5, 1.9e-5])
    pl.grid(True)
    pl.savefig('Data/Temp/wldiffs.png')
    pl.close()

# -----------------------------
# - INJECTION BEFORE TELLURIC -
# -----------------------------
    pl.figure(figsize=(18, 5))
    pl.plot(ximp, yimp, 'rx', ximp, yimp**rdt, 'kx', markersize=1)
    pl.show()
    pl.close()
    pl.figure(figsize=(18, 5))
    pl.plot(ximp, (1.0-yimp**rdt)/yimp, 'rx', markersize=1)
    pl.show()
    pl.close()

# ----------------------------
# - INJECTION AFTER TELLURIC -
# ----------------------------
    pcomp = [[], []]
    rdt = 0.0025
    for i in range(len(wlen[0])):
        pcomp[0].append([])
        pcomp[1].append([])
        for j in range(len(wlen)):
            ins = np.copy(10.0*wlen[j][i])
            pcomp[0][i] = np.insert(pcomp[0][i], len(pcomp[0][i]), list(ins))
            mflux = interpolate.splev(ins + slist[i], rep, der=0)
            ins = np.copy((mflux**rdt)*data[2*j+1][i])
            pcomp[1][i] = np.insert(pcomp[1][i], len(pcomp[1][i]), list(ins))
    pcomp = np.array(pcomp)

# ------------------------
# - FIRST MODEL RECOVERY -
# ------------------------
    pdiv = 1.816e-5
    pamp = 100.0*pdiv
    slist = np.arange(-pamp+pdiv, pamp, pdiv)
    vlist = slist*(sol/wl)

    other = []
    for i in range(len(pcomp[0])):
        other.append([])
        for s in slist:
            mflux = interpolate.splev(pcomp[0][i] + wldel[i] + s, rep, der=0)
            other[i].append(sum((mflux-pcomp[1][i])**2.0))
    for i in range(len(other)):
        other[i] /= np.median(other[i])
    other = np.array(other)

    imptnt = np.sum(other, axis=0)
    pl.plot(range(len(imptnt)), imptnt, 'kx')
    pl.grid(True)
    pl.savefig('Data/Temp/alione.png')
    pl.close()

    temp = fits.HDUList(fits.PrimaryHDU())
    temp.append(cp.deepcopy(fits.ImageHDU(data=other, name='REC1')))
    temp.writeto('Data/Temp/t07reco.fits', overwrite=True)

    fig, ax = pl.subplots(nrows=1, ncols=1, figsize=(18, 3))
    ax.imshow(other, interpolation='nearest', cmap='gray', aspect='auto')
    ax.xaxis.set_ticks(np.arange(9, len(other[0]), 15))
    gent = np.round(vlist[ax.get_xticks()], 0).astype(int)
    ax.set_xticklabels(map(str, gent), fontsize=12)
    ax.set_xlabel(r'$Radial\,\,velocity-RV_p$', fontsize=16)
    ax.yaxis.set_ticks(np.arange(0, 33, 5, dtype=int))
    ax.set_yticklabels(map(str, ax.get_yticks()), fontsize=12)
    ax.set_ylabel(r'$Spectrum\,\,number$', fontsize=16)
    ax.invert_yaxis()
    pl.tight_layout()
    pl.savefig('Data/Temp/recone.png')
    pl.close()

# -----------------------
# - OLD WATERFALL PLOTS -
# -----------------------
    pwfal = []
    for i in range(int(len(pextr)/2)):
        pwfal.append([])
        pwfal[i] = np.copy(pextr[2*i])
        pwfal[i] = np.insert(pwfal[i], len(pwfal[i]), np.nan, axis=0)
        pwfal[i] = np.insert(pwfal[i], len(pwfal[i]),
                             np.copy(pextr[2*i+1]), axis=0)
        pwfal[i] = np.copy(pextr[2*i+1])
        pwfal[i] = np.insert(pwfal[i], len(pwfal[i]), np.nan, axis=0)
        pwfal[i] = np.insert(pwfal[i], len(pwfal[i]),
                             np.copy(pcorr[2*i+1]), axis=0)
        pwfal[i] = np.insert(pwfal[i], len(pwfal[i]), np.nan, axis=0)
        pwfal[i] = np.insert(pwfal[i], len(pwfal[i]),
                             np.copy(palig[2*i+1]), axis=0)
        pwfal[i] = np.insert(pwfal[i], len(pwfal[i]), np.nan, axis=0)
        pwfal[i] = np.insert(pwfal[i], len(pwfal[i]),
                             np.copy(paraa[2*i+1]), axis=0)
        pwfal[i] = np.insert(pwfal[i], len(pwfal[i]), np.nan, axis=0)
        pwfal[i] = np.insert(pwfal[i], len(pwfal[i]),
                             np.copy(ptraa[2*i+1]), axis=0)
        pwfal[i] = np.insert(pwfal[i], len(pwfal[i]), np.nan, axis=0)
        pwfal[i] = np.insert(pwfal[i], len(pwfal[i]),
                             np.copy(psraa[2*i+1]), axis=0)
    fname = 'Data/Temp/Old/P' + tag + '/oldwfal' + tag + '.fits'
    fc.flis(pwfal, 'NORM').writeto('Data/Temp/.fits', overwrite=True)
