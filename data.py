#!/usr/bin/env python

import glob as gb
import func as fc
import numpy as np
import argparse as ap
import matplotlib.pyplot as pl
import matplotlib.ticker as tk
from astropy.io import fits
from scipy import interpolate
from scipy import ndimage
from scipy import stats
from astropy.analytic_functions.blackbody import blackbody_lambda


if __name__ == '__main__':
    # Opens the software and gives simple instructions to the user.
    print('--- SPECTRA ANALYSIS ---')
    print('Type 1 to accept a step or anything else to reject it.')

    # ---------------------------
    # - FILE'S NAMES EXTRACTION -
    # ---------------------------
    # #
    # Selects the source files and creates a list with all their names {1} or
    # simply uses the previously created list of file's names {other strings}.
    # #
    v = input(' - File\'s names extraction... ')
    fname = 'Data/Temp/flist.txt'
    if v is '1':
        flist = np.array(sorted(gb.glob('Data/*.fits')))
        with open(fname, 'w') as doc:
            [doc.write(name + '\n') for name in flist]
        print('     Process complete (%d files).' % len(flist))
    else:
        flist = np.array(np.genfromtxt(fname, dtype=str))
        print('     Read-in complete (%d files).' % len(flist))

    # ---------------------------
    # - WAVELENGTH EXTRACTION -
    # ---------------------------
    # #
    # Extracts the wavelength solution from the provided .fits file {1} or
    # simply uses the previously extracted values {other strings}.
    # #
    v = input(' - Wavelength extraction... ')
    fname = 'Data/Notes/data3.fits'
    tname = 'Data/Temp/t01wlen.fits'
    if v is '1':
        with fits.open(fname) as doc:
            pwave = 10*doc[1].data['WLEN'][0]
        fc.flis(pwave, 'NORM').writeto(tname, overwrite=True)
        print('     Process complete.')
    else:
        pwave = fc.frea(tname)
        print('     Read-in complete.')

    # ----------------------
    # - SPECTRA EXTRACTION -
    # ----------------------
    # #
    # Extracts the optimized flux values from the chips and generates the .fits
    # file with two stacked spectra (original and normalized versions) per chip
    # {1} or simply uses the previously extracted values {other strings}.
    # #
    v = input(' - Spectra extraction... ')
    fname = 'Images/s01extr.fits'
    if v is '1':
        pextr = fc.fima(flist, 'Extracted_OPT')
        fc.flis(pextr, 'NORM').writeto(fname, overwrite=True)
        print('     Process complete.')
    else:
        pextr = fc.frea(fname)
        print('     Read-in complete.')

    # -------------------------
    # - FLUX ERROR EXTRACTION -
    # -------------------------
    # #
    # Extracts the optimized flux errors from the chips and generates the .fits
    # file with all the stacked errors {1} or simply uses the previously
    # extracted values {other strings}.
    # #
    v = input(' - Flux errors extraction... ')
    fname = 'Data/Temp/t02erro.fits'
    if v is '1':
        perro = fc.ferr(flist, 'Error_OPT')
        fc.flis(perro, 'NORM').writeto(fname, overwrite=True)
        print('     Process complete.')
    else:
        perro = fc.frea(fname)
        print('     Read-in complete.')

    # ------------------------
    # - HISTOGRAM GENERATION -
    # ------------------------
    # #
    # Generates a histogram with all the flux values per chip {1} or does
    # nothing {other strings}.
    # #
    v = input(' - Histograms generation... ')
    if v is '1':
        fc.fhis(pextr, [0, 1.6e5, 0, 2500], 'c0')
        print('     Process complete.')
    else:
        print('     Process skipped.')

    # ------------------
    # - SIGMA CLIPPING -
    # ------------------
    # #
    # Creates normalized images, where each row is divided by its own median
    # and each column is subtracted by its own average and then divided by the
    # standard deviation {1}. The rejected pixels (coming from this sigma
    # clipping selection) are also flagged on the second images, for an easy
    # detection. Or simply uses the previous data processing {other strings}.
    # #
    v = input(' - Sigma clipping... ')
    fname = 'Data/Temp/t03clip.fits'
    if v is '1':
        psigm = fc.fsig(pextr, 5.0, 'Data/Temp/bpaut.txt')
        fc.flis(psigm, 'FLAG').writeto(fname, overwrite=True)
        print('     Process complete.')
    else:
        psigm = fc.frea(fname)
        print('     Read-in complete.')

    # ---------------------
    # - BAD PIXEL LOADING -
    # ---------------------
    # #
    # Generates an organized list with all the coordinates of the bad pixels,
    # both to be masked (marked with 0) and interpolated (marked with 1).
    # #
    print(' - Bad pixels loading...')
    pbpix = fc.fbad(pextr.shape, 'Data/Temp/bpman.txt')
    print('     Process complete.')

    # ----------------------
    # - SPECTRA CORRECTION -
    # ----------------------
    # #
    # Corrects the original images, interpolating some of the bad pixels and
    # masking the ones that cannot be interpolated {1}. It also creates the
    # .fits file with the corrected spectra (two per chip as before). Or simply
    # uses the previous data processing {other strings}.
    # #
    v = input(' - Spectra correction... ')
    fname = 'Images/s02corr.fits'
    if v is '1':
        pcorr = fc.fcor(pextr, pbpix)
        fc.flis(pcorr, 'NORM').writeto(fname, overwrite=True)
        print('     Process complete.')
    else:
        pcorr = fc.frea(fname)
        print('     Read-in complete.')

    # ---------------------
    # - SPECTRA ALIGNMENT -
    # ---------------------
    # #
    # Aligns the spectra (previously cleaned out of all the bad pixels), by
    # minimizing their squared differences relatively to a certain fixed
    # spectrum - the one with the lowest signal-to-noise ratio {1}. It also
    # creates the .fits file with the aligned spectra (two per chip as before).
    # Or simply uses the previous data processing {other strings}.
    # #
    v = input(' - Spectra alignment... ')
    fname = 'Images/s03alig.fits'
    tname = 'Data/Temp/alig.txt'
    if v is '1':
        palig, salig = fc.fali(pcorr, perro, pbpix, 0.001, 2, tname)
        fc.flis(palig, 'NORM').writeto(fname, overwrite=True)
        print('     Process complete.')
    else:
        palig = fc.frea(fname)
        salig = np.swapaxes(np.array(np.genfromtxt(tname), dtype=None), 0, 1)
        print('     Read-in complete.')

    # -------------------------
    # -                       -
    # - LEAST SQUARES - MODEL -
    # -                       -
    # -------------------------

    # MODEL INJECTION
    # Extraction of the observations' time stamps.
    obs_times = fc.fhea(flist, ['MJD-OBS'], 0)

    # Definition of the important physical parameters.
    rdt = 1.0
    tag = str(rdt)
    period = 3.52474859
    sol, wl = 299792.458, 3.5
    t_zero = 2452826.628521-2400000.5

    # Calculation of the observed orbital phases.
    phase = ((obs_times[0]-t_zero)/period) % 1

    mname = 'Data/Models/spec_1-1-1-1-1-1.dat'
    model = np.flip(np.swapaxes(np.genfromtxt(mname, dtype=float), 0, 1), 1)

    sigma = 3/(2*(2*np.log(2))**0.5)
    xmdl, ymdl = model[0], ndimage.filters.gaussian_filter(model[1], sigma)
    rep = interpolate.splrep(xmdl, ymdl, s=0)

    kp, vsys, vbar = 135.0, -140.0, 0.0
    wldel = (wl/sol)*(kp*np.sin(2*np.pi*phase)+vsys+vbar)

    ptemp = []
    for i in range(int(len(palig)/2)):
        ptemp.append([])
        ptemp.append([])
        for j in range(len(palig[2*i+1])):
            ins = np.copy(pwave[i][j])
            pflux = interpolate.splev(ins + wldel[j], rep, der=0)
            sflux = blackbody_lambda((ins + wldel[j])*10**4, 6071)
            inje = 0.1*(pflux/sflux)*((1.41*7.1492e7)/(1.2*6.957e8))**2
            ins = (1+rdt*np.array(inje, dtype=float))*np.copy(palig[2*i][j])
            mdn = np.median(ins[ins != 0])
            ptemp[2*i].append(ins)
            ptemp[2*i+1].append(ins/mdn)
    palig = np.array(ptemp)
    fname = 'Data/Temp/P' + tag + '/t04inje' + tag + '.fits'
    fc.flis(palig, 'NORM').writeto(fname, overwrite=True)

    # TELLURIC FEATURES REMOVAL
    # Data to be collected from the header.
    elist = ['HIERARCH ESO TEL AIRM END',
             'HIERARCH ESO TEL AIRM START',
             'HIERARCH ESO TEL AMBI TEMP',
             'HIERARCH ESO TEL AMBI FWHM END',
             'HIERARCH ESO TEL AMBI FWHM START']
    head = fc.fhea(flist, elist, 0)
    # Airmass will be the mean of the two available values.
    am = (head[0]+head[1])*0.5
    # Temperature will be the only value available.
    tp = head[2]
    # Seeing will be the mean of the two available values, corrected from the
    # unexpected -1 values.
    si = []
    for i in range(len(head[3])):
        if head[3][i] != -1 and head[4][i] != -1:
            si.append((head[3][i]+head[4][i])*0.5)
        if head[3][i] == -1:
            si.append(head[4][i])
        if head[4][i] == -1:
            si.append(head[3][i])
    si = np.array(si)

    # Check which is the best fitting model (AM, TP, SI).
    lin = [[], [], []]
    log = [[], [], []]
    linlog = [[], [], []]
    loglin = [[], [], []]
    for i in range(int(len(palig)/2)):
        for j in range(len(palig[i][0])):
            flux = np.copy(palig[2*i, :, j])
            if not np.array_equal(flux, np.zeros(len(flux))):
                lin[0].append(stats.linregress(am, flux)[2]**2)
                log[0].append(stats.linregress(np.log10(am),
                                               np.log10(flux))[2]**2)
                linlog[0].append(stats.linregress(am, np.log10(flux))[2]**2)
                loglin[0].append(stats.linregress(np.log10(am), flux)[2]**2)
                lin[1].append(stats.linregress(tp, flux)[2]**2)
                log[1].append(stats.linregress(np.log10(tp),
                                               np.log10(flux))[2]**2)
                linlog[1].append(stats.linregress(tp, np.log10(flux))[2]**2)
                loglin[1].append(stats.linregress(np.log10(tp), flux)[2]**2)
                lin[2].append(stats.linregress(si, flux)[2]**2)
                log[2].append(stats.linregress(np.log10(si),
                                               np.log10(flux))[2]**2)
                linlog[2].append(stats.linregress(si, np.log10(flux))[2]**2)
                loglin[2].append(stats.linregress(np.log10(si), flux)[2]**2)

    print('\nAIRMASS CORRELATIONS (AM VS. FL):')
    print('LIN-LIN: '+str(np.mean(lin[0])))
    print('LOG-LOG: '+str(np.mean(log[0])))
    print('LIN-LOG: '+str(np.mean(linlog[0])))
    print('LOG-LIN: '+str(np.mean(loglin[0])))

    print('\nTEMPERATURE CORRELATIONS (TP VS. FL):')
    print('LIN-LIN: '+str(np.mean(lin[1])))
    print('LOG-LOG: '+str(np.mean(log[1])))
    print('LIN-LOG: '+str(np.mean(linlog[1])))
    print('LOG-LIN: '+str(np.mean(loglin[1])))

    print('\nSEEING CORRELATIONS (SI VS. FL):')
    print('LIN-LIN: '+str(np.mean(lin[2])))
    print('LOG-LOG: '+str(np.mean(log[2])))
    print('LIN-LOG: '+str(np.mean(linlog[2])))
    print('LOG-LIN: '+str(np.mean(loglin[2])))

    # # # AIRMASS REMOVAL BEFORE ALIGNMENT # # #
    parba = []
    for i in range(int(len(pcorr)/2)):
        parba.append([])
        parba.append([])
        for j in range(len(pcorr[2*i][0])):
            flux = np.copy(pcorr[2*i, :, j])
            if np.array_equal(flux, np.zeros(len(flux))):
                parba[2*i].append(flux)
                parba[2*i+1].append(flux)
            else:
                linr_linlog = stats.linregress(am, np.log10(flux))
                model = (10**linr_linlog[1])*(10**(linr_linlog[0]*am))
                parba[2*i].append(flux/model)
                parba[2*i+1].append(flux/model)
    parba = np.swapaxes(np.array(parba), 1, 2)
    for i in range(int(len(parba)/2)):
        for j in range(len(parba[2*i+1])):
            mdn = np.median(parba[2*i+1][j][parba[2*i+1][j] != 0])
            parba[2*i+1][j] /= mdn
    fname = 'Data/Temp/P' + tag + '/m01arba' + tag + '.fits'
    fc.flis(parba, 'NORM').writeto(fname, overwrite=True)

    # # # AIRMASS REMOVAL AFTER ALIGNMENT # # #
    paraa = []
    for i in range(int(len(palig)/2)):
        paraa.append([])
        paraa.append([])
        for j in range(len(palig[2*i][0])):
            flux = np.copy(palig[2*i, :, j])
            if np.array_equal(flux, np.zeros(len(flux))):
                paraa[2*i].append(flux)
                paraa[2*i+1].append(flux)
            else:
                linr_linlog = stats.linregress(am, np.log10(flux))
                model = (10**linr_linlog[1])*(10**(linr_linlog[0]*am))
                paraa[2*i].append(flux/model)
                paraa[2*i+1].append(flux/model)
    paraa = np.swapaxes(np.array(paraa), 1, 2)
    for i in range(int(len(paraa)/2)):
        for j in range(len(paraa[2*i+1])):
            mdn = np.median(paraa[2*i+1][j][paraa[2*i+1][j] != 0])
            paraa[2*i+1][j] /= mdn
    fname = 'Data/Temp/P' + tag + '/m02araa' + tag + '.fits'
    fc.flis(paraa, 'NORM').writeto(fname, overwrite=True)

    # # # TEMPERATURE REMOVAL BEFORE ALIGNMENT # # #
    ptrba = []
    for i in range(int(len(parba)/2)):
        ptrba.append([])
        ptrba.append([])
        for j in range(len(parba[2*i][0])):
            flux = np.copy(parba[2*i, :, j])
            if np.array_equal(flux, np.zeros(len(flux))):
                ptrba[2*i].append(flux)
                ptrba[2*i+1].append(flux)
            else:
                linr_loglin = stats.linregress(np.log10(tp), flux)
                model = linr_loglin[0]*np.log10(tp)+linr_loglin[1]
                ptrba[2*i].append(flux/model)
                ptrba[2*i+1].append(flux/model)
    ptrba = np.swapaxes(np.array(ptrba), 1, 2)
    for i in range(int(len(ptrba)/2)):
        for j in range(len(ptrba[2*i+1])):
            mdn = np.median(ptrba[2*i+1][j][ptrba[2*i+1][j] != 0])
            ptrba[2*i+1][j] /= mdn
    fname = 'Data/Temp/P' + tag + '/m03trba' + tag + '.fits'
    fc.flis(ptrba, 'NORM').writeto(fname, overwrite=True)

    # # # TEMPERATURE REMOVAL AFTER ALIGNMENT # # #
    ptraa = []
    for i in range(int(len(paraa)/2)):
        ptraa.append([])
        ptraa.append([])
        for j in range(len(paraa[2*i][0])):
            flux = np.copy(paraa[2*i, :, j])
            if np.array_equal(flux, np.zeros(len(flux))):
                ptraa[2*i].append(flux)
                ptraa[2*i+1].append(flux)
            else:
                linr_loglin = stats.linregress(np.log10(tp), flux)
                model = linr_loglin[0]*np.log10(tp)+linr_loglin[1]
                ptraa[2*i].append(flux/model)
                ptraa[2*i+1].append(flux/model)
    ptraa = np.swapaxes(np.array(ptraa), 1, 2)
    for i in range(int(len(ptraa)/2)):
        for j in range(len(ptraa[2*i+1])):
            mdn = np.median(ptraa[2*i+1][j][ptraa[2*i+1][j] != 0])
            ptraa[2*i+1][j] /= mdn
    fname = 'Data/Temp/P' + tag + '/m04traa' + tag + '.fits'
    fc.flis(ptraa, 'NORM').writeto(fname, overwrite=True)

    # # # SEEING REMOVAL BEFORE ALIGNMENT # # #
    psrba = []
    for i in range(int(len(ptrba)/2)):
        psrba.append([])
        psrba.append([])
        for j in range(len(ptrba[2*i][0])):
            flux = np.copy(ptrba[2*i, :, j])
            if np.array_equal(flux, np.zeros(len(flux))):
                psrba[2*i].append(flux)
                psrba[2*i+1].append(flux)
            else:
                linr_linlog = stats.linregress(si, np.log10(flux))
                model = (10**linr_linlog[1])*(10**(linr_linlog[0]*si))
                psrba[2*i].append(flux/model)
                psrba[2*i+1].append(flux/model)
    psrba = np.swapaxes(np.array(psrba), 1, 2)
    for i in range(int(len(psrba)/2)):
        for j in range(len(psrba[2*i+1])):
            mdn = np.median(psrba[2*i+1][j][psrba[2*i+1][j] != 0])
            psrba[2*i+1][j] /= mdn
    fname = 'Data/Temp/P' + tag + '/m05srba' + tag + '.fits'
    fc.flis(psrba, 'NORM').writeto(fname, overwrite=True)

    # # # SEEING REMOVAL AFTER ALIGNMENT # # #
    psraa = []
    for i in range(int(len(ptraa)/2)):
        psraa.append([])
        psraa.append([])
        for j in range(len(ptraa[2*i][0])):
            flux = np.copy(ptraa[2*i, :, j])
            if np.array_equal(flux, np.zeros(len(flux))):
                psraa[2*i].append(flux)
                psraa[2*i+1].append(flux)
            else:
                linr_linlog = stats.linregress(si, np.log10(flux))
                model = (10**linr_linlog[1])*(10**(linr_linlog[0]*si))
                psraa[2*i].append(flux/model)
                psraa[2*i+1].append(flux/model)
    psraa = np.swapaxes(np.array(psraa), 1, 2)
    for i in range(int(len(psraa)/2)):
        for j in range(len(psraa[2*i+1])):
            mdn = np.median(psraa[2*i+1][j][psraa[2*i+1][j] != 0])
            psraa[2*i+1][j] /= mdn
    fname = 'Data/Temp/P' + tag + '/m06sraa' + tag + '.fits'
    fc.flis(psraa, 'NORM').writeto(fname, overwrite=True)

    pclea = fc.frea('Data/Temp/P' + tag + '/m06sraa' + tag + '.fits')

    # ''
    # WATERFAL PLOTS #
    # # # ARGUMENTS # # #
    fc.flis(psraa, 'NORM').writeto('Images/s04clea.fits', overwrite=True)
    f = [[0, 1], [1, 1], [2, 1], [2, 0]]
    pgrou = [pclea, palig, pextr]

    # # # PARAMETERS # # #
    spec_num = len(f)
    rows, cols = len(pgrou[0][0]), len(pgrou[0][0][0])
    ratio = float(rows)/float(cols)
    bot, top, lef, rig = 0.15, 0.25, 0.05, 0.05
    quant = (1-(bot+top))/float(spec_num)
    fig_w = 18
    fig_h = ratio*float(fig_w*spec_num)*((1-(lef+rig))/(1-(bot+top)))

    for i in range(int(len(pgrou[0])/2)):
        # # # FIGURE # # #
        fig = pl.figure(figsize=(fig_w, fig_h))
        fig.suptitle(r'$Chip\,\,%d$' % (i+1), fontsize=24)

        # # # GENERAL AXIS # # #
        gb = fig.add_axes([lef, bot, 1-(lef+rig), 1-(bot+top)])
        gb.axis([0, cols-1, 0, rows-1])
        gb.xaxis.set_ticks(np.arange(15, cols, 100))
        gent = np.round(pwave[i][0][gb.get_xticks()], 3)
        gb.set_xticklabels(map(str, gent), fontsize=8)
        gb.set_xlabel(r'$Wavelength-\lambda\,\,(\mu m)$', fontsize=12)

        gt = gb.twiny()
        gt.axis([0, cols-1, 0, rows-1])
        gt.xaxis.set_ticks(np.arange(0, cols, 200))
        gt.set_xticklabels(map(str, gt.get_xticks()), fontsize=8)
        gt.set_xlabel(r'$Pixel\,\,number$', fontsize=12)

        gb.yaxis.set_ticks([])
        gb.set_ylabel(r'$Spectrum\,\,number$', fontsize=12, labelpad=20)

        # # # SUBPLOTS # # #
        sub = []
        for j in range(spec_num):
            a = pgrou[f[j][0]][2*i+f[j][1]]
            per, ipol = 2, 'nearest'
            p = np.percentile(a[a != 0], per)
            q = np.percentile(a[a != 0], 100-per)
            sub.append(fig.add_axes([lef, bot+j*quant, 1-(lef+rig), quant]))
            sub[j].imshow(a, interpolation=ipol, cmap='gray').set_clim(p, q)
            sub[j].xaxis.set_ticks([])
            sub[j].yaxis.set_ticks(np.arange(0, 33, 11, dtype=int))
            sub[j].set_yticklabels(map(str, sub[j].get_yticks()), fontsize=8)
            sub[j].invert_yaxis()

        iname = 'Images/c0' + str(i+1) + 'wfal.png'
        pl.savefig(iname, dpi=500)
        pl.close()
    # ''

    # LEAST SQUARES PROCESS
    pcomp = [[], []]
    for i in range(len(pwave[0])):
        pcomp[0].append([])
        pcomp[1].append([])
        for j in range(len(pwave)):
            arr = [pclea[2*j+1][i] != 0]
            ins = np.copy(pwave[j][i][arr])
            pcomp[0][i] = np.insert(pcomp[0][i], len(pcomp[0][i]), ins)
            ins = np.copy(pclea[2*j+1][i][arr])
            pcomp[1][i] = np.insert(pcomp[1][i], len(pcomp[1][i]), ins)
    pcomp = np.array(pcomp)

    pdiv, ampv = 1.816e-5, 400
    pamp = pdiv*ampv
    slist = np.arange(-pamp+pdiv, pamp, pdiv)
    vlist = slist*(sol/wl)

    npile = []
    for i in range(len(pcomp[0])):
        npile.append([])
        for s in slist:
            pflux = interpolate.splev(pcomp[0][i] + s, rep, der=0)
            npile[i].append(sum((pflux-pcomp[1][i])**2))
    for i in range(len(npile)):
        npile[i] /= np.median(npile[i])
    ml = np.mean(npile, axis=0)
    for i in range(len(npile[0])):
        for j in range(len(npile)):
            npile[j][i] -= ml[i]
    sl = np.std(npile, axis=0)
    for i in range(len(npile[0])):
        for j in range(len(npile)):
            if sl[i] != 0:
                npile[j][i] /= sl[i]
    npile = np.array(npile)

    fname = 'Data/Temp/P' + tag + '/t05mode' + tag + '.fits'
    temp = fits.HDUList(fits.PrimaryHDU())
    temp.append(fits.ImageHDU(data=npile, name='MOD1'))
    temp.writeto(fname, overwrite=True)

    fig, ax = pl.subplots(nrows=1, ncols=1, figsize=(18, 3))
    ax.imshow(npile, interpolation='nearest', cmap='gray', aspect='auto')
    ax.xaxis.set_ticks(np.arange(9, len(npile[0]), 30))
    gent = np.round(vlist[ax.get_xticks()], 0).astype(int)
    ax.set_xticklabels(map(str, gent), fontsize=12)
    ax.set_xlabel(r'$Radial\,\,velocity-RV_p$', fontsize=16)
    ax.yaxis.set_ticks(np.arange(0, 33, 5, dtype=int))
    ax.set_yticklabels(map(str, ax.get_yticks()), fontsize=12)
    ax.set_ylabel(r'$Spectrum\,\,number$', fontsize=16)
    ax.invert_yaxis()
    pl.tight_layout()
    pl.savefig('Data/Temp/P' + tag + '/mod' + tag + '.png', dpi=500)
    pl.close()

    # LUCKY MODEL RECOVERY
    wrg = ''  # 'wrg'
    wldel = (wl/sol)*(kp*np.sin(2*np.pi*phase)+vsys)

    cutv = int(max(np.absolute(wldel))/pdiv)
    pampn = (ampv-cutv)*pdiv
    snew = np.arange(-pampn+pdiv, pampn, pdiv)
    vnew = snew*(sol/wl)

    other = []
    for i in range(len(pcomp[0])):
        rep = interpolate.splrep(slist - wldel[i], npile[i], s=0)
        new = interpolate.splev(snew, rep, der=0)
        other.append(np.array(new))
    other = np.array(other)

    imptnt = -np.sum(other, axis=0)
    fig, ax = pl.subplots(nrows=1, ncols=1)
    ax.plot(range(len(imptnt)), imptnt, 'kx')
    ax.xaxis.set_ticks(np.arange(8, len(other[0]), 95))
    gent = np.round(vnew[ax.get_xticks()], 0).astype(int)
    ax.set_xticklabels(map(str, gent), fontsize=12)
    ax.set_xlabel(r'$System\,\,velocity-V_{sys}$', fontsize=16)
    ax.grid(True)
    pl.tight_layout()
    pl.savefig('Data/Temp/P' + tag + '/sum' + tag + wrg + '.png', dpi=500)
    pl.close()

    fname = 'Data/Temp/P' + tag + '/t06reco' + tag + wrg + '.fits'
    temp = fits.HDUList(fits.PrimaryHDU())
    temp.append(fits.ImageHDU(data=other, name='REC1'))
    temp.writeto(fname, overwrite=True)

    fig, ax = pl.subplots(nrows=1, ncols=1, figsize=(18, 3))
    ax.imshow(other, interpolation='nearest', cmap='gray', aspect='auto')
    ax.xaxis.set_ticks(np.arange(8, len(other[0]), 95))
    gent = np.round(vnew[ax.get_xticks()], 0).astype(int)
    ax.set_xticklabels(map(str, gent), fontsize=12)
    ax.set_xlabel(r'$System\,\,velocity-V_{sys}$', fontsize=16)
    ax.yaxis.set_ticks(np.arange(0, 33, 5, dtype=int))
    ax.set_yticklabels(map(str, ax.get_yticks()), fontsize=12)
    ax.set_ylabel(r'$Spectrum\,\,number$', fontsize=16)
    ax.invert_yaxis()
    pl.tight_layout()
    pl.savefig('Data/Temp/P' + tag + '/rec' + tag + wrg + '.png', dpi=500)
    pl.close()

    # MODEL RECOVERY
    kpdiv, kpampv = 1.5, 200
    kpamp = kpdiv*kpampv
    kplist = np.arange(-kpamp+kpdiv, kpamp, kpdiv)

    cutv = []
    for i in range(len(kplist)):
        wldel = (wl/sol)*(kplist[i]*np.sin(2*np.pi*phase))
        cutv.append(int(max(np.absolute(wldel))/pdiv))
    cutv = max(cutv)

    final = []
    for i in range(len(kplist)):
        wldel = (wl/sol)*(kplist[i]*np.sin(2*np.pi*phase))
        pampn = (ampv-cutv)*pdiv
        snew = np.arange(-pampn+pdiv, pampn, pdiv)
        vnew = snew*(sol/wl)

        now = []
        for j in range(len(pcomp[0])):
            rep = interpolate.splrep(slist - wldel[j], npile[j], s=0)
            arr = interpolate.splev(snew, rep, der=0)
            now.append(np.array(arr))
        final.append(np.sum(np.array(now), axis=0))
    final = np.array(final)/np.std(final)
    print(np.unravel_index(final.argmin(), final.shape), final.min(), rdt)

    fname = 'Data/Temp/P' + tag + '/t07fina' + tag + '.fits'
    temp = fits.HDUList(fits.PrimaryHDU())
    temp.append(fits.ImageHDU(data=final, name='FINAL1'))
    temp.writeto(fname, overwrite=True)

    fig, ax = pl.subplots(nrows=1, ncols=1, figsize=(8, 8))
    ax.imshow(final, interpolation='nearest', cmap='gray', origin='lower')
    ax.xaxis.set_ticks(np.arange(8, len(final[0]), 30))
    gent = np.round(vnew[ax.get_xticks()], 0).astype(int)
    ax.set_xticklabels(map(str, gent), fontsize=12)
    ax.set_xlabel(r'$System\,\,velocity-V_{sys}$', fontsize=16)
    ax.yaxis.set_ticks(np.arange(19, len(final), 30))
    gent = np.round(kplist[ax.get_yticks()], 0).astype(int)
    ax.set_yticklabels(map(str, gent), fontsize=12)
    ax.set_ylabel(r'$Radial\,\,velocity-K_p$', fontsize=16)
    pl.tight_layout()
    pl.savefig('Data/Temp/P' + tag + '/map' + tag + '.png')
    pl.close()
