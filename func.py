#!/usr/bin/env python

import copy as cp
import glob as gb
import numpy as np
import argparse as ap
import matplotlib.pyplot as pl
import matplotlib.ticker as tk
from astropy.io import fits
from scipy import interpolate
from scipy import stats


# ---------------------------
# - FILE'S NAMES EXTRACTION -
# ---------------------------
def fnam():
    '''
    Returns a list of file's names, given a specific argument.
     < temp: list of file's names
    '''

    # Receives the extensions of all the needed files.
    parser = ap.ArgumentParser()
    parser.add_argument('files', metavar='FILE', type=str,
                        nargs='+', help='File extension.')

    # Lists the corresponding file's names in a single array.
    temp = []
    for name in parser.parse_args().files:
        temp += cp.deepcopy(gb.glob(name))

    return np.copy(temp)


# ----------------------
# - SPECTRA EXTRACTION -
# ----------------------
def fima(flist, col):
    '''
    Returns an array of images (original and normalized), given the spectra.
     > flist: list of file's names
     > col: column's name
     < npile: array (or pile) of images
    '''

    # Calculates the number of chips used to capture the spectra.
    with fits.open(flist[0]) as mome:
        cnum = len(mome)-1

    # Creates an empty main array and, for each chip, creates two other arrays
    # within the main one - these will contain the stacked spectra (both the
    # original and the normalized).
    npile = []
    for i in range(cnum):
        npile.append([])
        npile.append([])
        # For each .fits file (corresponding to a single spectrum), stores its
        # spectrum, both the original and the one normalized to the median.
        for name in flist:
            with fits.open(name) as doc:
                got = np.copy(doc[i+1].data.field(col))
                npile[2*i].append(list(got))
                npile[2*i+1].append(list(got/np.median(got[got != 0])))

    return np.copy(npile)


# -----------------------
# - HDU LIST GENERATION -
# -----------------------
def flis(data, tag):
    '''
    Returns a HDU list, given a simple pile of images.
     > data: array (or pile) of images
     > tag: addition to the name of the possible second version of each image
     < temp: HDU list
    '''

    # Creates an empty list of HDU objects.
    temp = fits.HDUList(fits.PrimaryHDU())

    # For each chip (slice of the array of images), creates some HDU objects.
    if len(data) == 4:
        for i in range(len(data)):
            title = 'CHIP' + str(int(i+1))
            temp.append(cp.deepcopy(fits.ImageHDU(data=data[i], name=title)))
    if len(data) == 8:
        for i in range(len(data)):
            title = 'CHIP' + str(int(i*0.5+1))
            if i % 2 == 1:
                title += '_' + tag
            temp.append(cp.deepcopy(fits.ImageHDU(data=data[i], name=title)))

    return temp


# --------------------
# - IMAGE COLLECTION -
# --------------------
def frea(fname):
    '''
    Returns an array of images, given a previously created .fits file.
     > fname: name of the source .fits file
     < npile: array (or pile) of images
    '''

    # Creates an empty array and uses it to collect all the images previously
    # captured and piled up in the respective .fits file.
    npile = []
    with fits.open(fname) as doc:
        [npile.append(doc[i+1].data) for i in range(len(doc)-1)]

    return np.copy(npile)


# --------------------
# - ERROR COLLECTION -
# --------------------
def ferr(flist, col):
    '''
    Returns an array of flux errors, given the spectra.
     > flist: list of file's names
     > col: column's name
     < npile: array (or pile) of errors
    '''

    # Calculates the number of chips used to capture the spectra.
    with fits.open(flist[0]) as mome:
        cnum = len(mome)-1

    # Creates an empty array which will contain all the flux errors per chip.
    npile = []
    for i in range(cnum):
        npile.append([])
        # For each .fits file (or spectrum), stores its errors.
        for name in flist:
            with fits.open(name) as doc:
                npile[i].append(list(np.copy(doc[i+1].data.field(col))))

    return np.copy(npile)


# ------------------------
# - HISTOGRAM GENERATION -
# ------------------------
def fhis(data, axis, tag):
    '''
    Produces one histogram per chip, containing all the flux values.
     > data: array (or pile) of images
     > axis: array with the edges of both axis
     > tag: addition to the name of the images
     < temp: information about all the produced histograms
    '''

    # For each chip, one histogram with all its flux values is created and the
    # corresponding histogram information is stored.
    temp = []
    for i in range(int(len(data)/2)):
        iname = 'Images/' + tag + str(i+1) + 'hist.png'
        temp.append(pl.hist(data[i*2].flatten(), bins='auto'))
        pl.title('Chip ' + str(i+1) + ' - Histogram of the flux values')
        pl.xlabel('Flux values')
        pl.ylabel('Frequency')
        pl.axis(axis)
        pl.grid(True)
        pl.savefig(iname, dpi=500)
        pl.close()

    return np.copy(temp)


# ------------------
# - SIGMA CLIPPING -
# ------------------
def fsig(data, sigma, ndoc):
    '''
    Returns a sigma clipping of bad pixels, given the data and a sigma value.
     > data: array (or pile) of original images
     > sigma: sigma value of the sigma clipping
     > ndoc: name of the document that will contain the list of bad pixels
     < npile: array (or pile) of images with sigma clipping
    '''

    # Creates a new pile that will contain the sigma-clipped images (for each
    # chip, the first image corresponds to the data with normalized rows - to
    # the respective median - and columns subtracted by their respective mean,
    # and the second image contains all the rejected pixels).
    npile = []
    for i in range(int(len(data)/2)):
        npile.append(cp.deepcopy(list(data[i*2+1])))
        npile.append(cp.deepcopy(list(data[i*2+1])))
    npile = np.copy(npile)

    # Calculates arrays with the average per column of data in every chip.
    ml, sl = [], []
    [ml.append(np.mean(npile[i*2], axis=0)) for i in range(int(len(npile)/2))]

    # For each chip, every column is subtracted by its own average and divided
    # by its standard deviation as well; the rows were already normalized to
    # their respective medians.
    for i in range(int(len(npile)/2)):
        for j in range(len(npile[i*2])):
            for k in range(len(npile[i*2][j])):
                npile[i*2][j][k] -= ml[i][k]
        sl.append(np.std(npile[i*2], axis=0))
        for j in range(len(npile[i*2])):
            for k in range(len(npile[i*2][j])):
                if sl[i][k] != 0.0:
                    npile[i*2][j][k] /= sl[i][k]

    # Flags rejected pixels (and also detects columns which have null standard
    # deviation) with one, and all the other pixels with 0.5, so the rejected
    # ones can easily be detected on the output images.
    with open(ndoc, 'w') as doc:
        for i in range(int(len(npile)/2)):
            npile[i*2+1] = list(np.copy(npile[i*2]))
        for i in range(int(len(npile)/2)):
            for j in range(len(npile[i*2+1])):
                for k in range(len(npile[i*2+1][j])):
                    if sl[i][k] == 0.0:
                        npile[i*2+1][j][k] = sigma
                    if abs(npile[i*2+1][j][k]) >= sigma:
                        npile[i*2+1][j][k] = 1.0
                        doc.write(str(i+1)+'\t'+str(k+1)+'\t'+str(j+1)+'\n')
                    else:
                        npile[i*2+1][j][k] = 0.5

    return np.copy(npile)


# ---------------------
# - BAD PIXEL LOADING -
# ---------------------
def fbad(size, ndoc):
    '''
    Returns the list of pixels to be corrected, given the known bad regions.
     > size: dimensions of the array (or pile) of images
     > ndoc: name of the document containing the list of bad pixels
     < temp: list of all the pixels that must be corrected
    '''

    # Opens the list of all the bad pixels (these were previously identified
    # manually and written down in a specific file).
    mome = np.loadtxt(ndoc, dtype='int')

    # Creates an empty array to store both the list of pixels that will be
    # masked and the list of pixels that will be interpolated.
    temp = [[], []]
    for i in range(int(size[0]/2)):
        temp[0].append([])
        temp[1].append([])
        for j in range(size[1]):
            temp[0][i].append([])
            temp[1][i].append([])

    # For each line in the manually created list of bad pixel regions, the bad
    # pixels' coordinates are added to an organized structure that can then be
    # used throughout the rest of the software in an easier way. The returned
    # object has two lists (0 - to mask; 1 - to interpolate), and each one of
    # these contains four different lists (corresponding to the four chips);
    # within each chip list there is a list of bad pixels per spectrum.
    for l in mome:
        for j in np.arange(l[2]-1, l[3], 1):
            for k in np.arange(l[4]-1, l[5], 1):
                if l[1] == 0:
                    temp[0][l[0]-1][k] = np.append(temp[0][l[0]-1][k], j)
                if l[1] == 1:
                    temp[1][l[0]-1][k] = np.append(temp[1][l[0]-1][k], j)

    return np.copy(temp)


# ----------------------
# - SPECTRA CORRECTION -
# ----------------------
def fcor(data, bpix):
    '''
    Returns an array of corrected images, given the originals and bad pixels.
     > data: array (or pile) of original images
     > bpix: list of all the pixels that must be corrected
     < npile: array (or pile) of corrected images
    '''

    # Creates a new pile that will contain the corrected images, without bad
    # pixels or bad regions (as before, a normalized version of the spectra
    # will also be created for each chip).
    npile = []
    for i in range(int(len(data)/2)):
        npile.append(cp.deepcopy(list(data[i*2])))
        npile.append(cp.deepcopy(list(data[i*2])))
    npile = np.copy(npile)

    # Interpolates the bad pixels in every chip, using cubic splines without
    # smoothing based only on the values coming from the good pixels.
    for i in range(len(bpix[1])):
        for j in range(len(bpix[1][i])):
            # Then, for each original image and for each spectrum a new x array
            # is created by subtracting the bad pixels to the total amount of
            # pixels (0 to 1023) and cubic splines without smoothing are gene-
            # rated based on these correct values.
            xogn = np.arange(0, len(data[i][j]), 1)
            xnew = np.append(bpix[0][i][j], bpix[1][i][j])
            xold = np.delete(xogn, xnew)
            yold = np.delete(npile[i*2][j], xnew)
            rep = interpolate.splrep(xold, yold, s=0)
            # Finally, python uses this type of interpolation to calculate what
            # the fluxes should be in all the bad pixels.
            ynew = interpolate.splev(xnew, rep, der=0)
            for k in range(len(xnew)):
                npile[i*2][j][int(xnew[k])] = cp.deepcopy(ynew[k])

    # Masks, in every chip, the bad pixels that cannot be interpolated, by
    # setting them all to zero.
    for i in range(len(bpix[0])):
        for j in range(len(bpix[0][i])):
            for k in bpix[0][i][j]:
                npile[i*2][j][int(k)] = 0

    # Creates the normalized version of the corrected spectra.
    for i in range(int(len(npile)/2)):
        for j in range(len(npile[i*2])):
            got = np.copy(npile[i*2][j])
            npile[i*2+1][j] = got/np.median(got[got != 0])

    return np.copy(npile)


# ---------------------
# - SPECTRA ALIGNMENT -
# ---------------------
def fali(data, erro, bpix, pdiv, pamp, tname):
    '''
    Returns an array of aligned images, given the corrected ones.
     > data: array (or pile) of corrected images
     > erro: array (or pile) of errors
     > bpix: list of all the pixels that must be corrected
     > pdiv: minimum fraction of pixel by which spectra can be shifted
     > pamp: maximim amplitude of pixels by which spectra can be shifted
     > tname: name of the file that will contain the pixel shifts
     < npile: array (or pile) of aligned images
     < shif: array of pixel shifts
    '''

    # Creates two arrays: one containing all the chip numbers and another one
    # containing all the possible amounts of shifts applicable to the spectra.
    clist = range(int(len(data)/2))
    plist = np.arange(-float(pamp)+pdiv, float(pamp), pdiv)

    # Creates two empty arrays the will store all the new aligned spectra and
    # all the pixel shifts that will be applied to them. After that, the cycle
    # that will be applied to every chip - in order to align them - is started.
    npile, shif = [], []
    for chip in clist:
        # Extends the previously created arrays so all the chips can be stored
        # and creates an array that will contain all the spectra shift indexes.
        npile.append([])
        npile.append([])
        shif.append([])
        slist = range(len(data[2*chip]))

        # Creates an array to store the sum (pixel by pixel) of the squared
        # flux-to-error (signal-to-noise) ratios for every spectrum, calculates
        # these sums and stores the spectrum number that maximizes the sums.
        total = []
        for s in slist:
            oaux = np.arange(0, len(data[2*chip][s]), 1)
            baux = np.append(bpix[0][chip][s], bpix[1][chip][s])
            naux = np.delete(oaux, baux)
            # Calculates the sum (pixel by pixel) of the squared flux-to-error
            # (signal-to-noise) ratios, corrected by the amount of good pixels.
            aux = sum((data[2*chip][s][naux]/erro[chip][s][naux])**2.0)
            total.append(aux/float(len(naux)))
            # Plots the obtained signal-to-noise per spectrum.
            for i in range(len(shif)):
                iname = 'Images/c0' + str(chip+1) + 'sgtn.png'
                pl.plot(range(len(total)), np.array(total)/max(total), 'rx')
                pl.title('Chip %d - Normalized signal-to-noise' % (chip+1))
                pl.xlabel('Spectrum number')
                pl.ylabel('Signal-to-noise')
                pl.axis([0, 35, 0.0, 1.0])
                pl.grid(True)
                pl.savefig(iname, dpi=500)
                pl.close()
        # Fixes the spectrum which shows a higher signal-to-noise value.
        spec = total.index(max(total))
        leng = len(data[2*chip][spec])

        # Creates arrays containing the pixel indexes and fluxes of the fixed
        # spectrum. Only pixels that are further from masked regions by, at
        # least, the chosen amplitude of pixels are considered - this avoids
        # the masked regions to contribute to the alignment process.
        xbad = []
        for pix in bpix[0][chip][spec]:
            for i in np.arange(pix-pamp, pix+pamp+1, 1):
                if i not in xbad and i >= 0 and i < leng:
                    xbad.append(cp.deepcopy(i))
        xfix = np.delete(np.arange(0, leng, 1), xbad)
        yfix = np.copy(data[2*chip][spec][xfix])

        # Moves all the other spectra sideways, in order to get the shift value
        # which minimizes the squared differences to the fixed spectrum.
        for s in slist:
            # Collects all the flux values of the considered spectrum.
            vals = np.copy(data[2*chip][s])

            # If this happens to be the fixed spectrum, then its flux values
            # are just directly copied to the new aligned group of spectra.
            if s == spec:
                # Stores zero as the pixel shift applied to this spectrum and
                # stores both the fixed spectrum and its normalized version.
                shif[chip].append(0)
                npile[2*chip].append(list(vals))
                npile[2*chip+1].append(list(vals/np.median(vals[vals != 0])))

            # If this spectrum has to be aligned, then the sum of the squared
            # differences are calculated across all the possible shifts that
            # can be applied to the sepctrum - always relatively to the fixed
            # spectrum - and the shift that minimizes this sum is chosen.
            else:
                # Creates arrays containing the pixel indexes and fluxes of the
                # considered spectrum, in which the previously masked regions
                # are not included, so that a spline representation of the
                # spectrum can be generated - this allows the interpolation.
                xogn = np.arange(0, len(vals), 1)
                xbad = np.copy(bpix[0][chip][s])
                x = np.delete(xogn, xbad)
                y = np.delete(vals, xbad)
                rep = interpolate.splrep(x, y, s=0)

                # Creates an array to store the sum of the squared differences
                # for every possible shift considered, calculates these sums
                # and stores the pixel shift that minimizes the sums.
                total = []
                for pix in plist:
                    # Generates interpolated flux values for the shifted pixels
                    # based on the spline representation of the spectrum.
                    ynew = interpolate.splev(xfix + pix, rep, der=0)
                    # Calculates the sum of the squared differences relatively
                    # to the fixed spectrum.
                    total.append(sum((yfix-ynew)**2.0))
                shif[chip].append(cp.deepcopy(plist[total.index(min(total))]))

                # Generates the definite interpolated flux values, considering
                # already the final pixel positions that align the two spectra.
                yfin = interpolate.splev(xogn + shif[chip][s], rep, der=0)
                # Sets all the pixels within the masked regions to zero.
                for i in range(len(vals)):
                    if i in bpix[0][chip][s]:
                        yfin[i] = 0
                # Stores both the fixed spectrum and its normalized version.
                npile[2*chip].append(list(yfin))
                npile[2*chip+1].append(list(yfin/np.median(yfin[yfin != 0])))

    # Plots the obtained pixel shifts that were applied to the spectra.
    for i in range(len(shif)):
        iname = 'Images/c0' + str(i+1) + 'alig.png'
        xbas, ybas = range(len(shif[i])), shif[i]
        line = stats.linregress(xbas, ybas)
        pnts = line[0]*xbas+line[1]
        pl.plot(xbas, ybas, 'rx', xbas, pnts, 'k-.', markersize=3)
        pl.title('Chip %d - Spectra alignment' % (i+1))
        pl.xlabel('Spectrum number')
        pl.ylabel('Pixel shift')
        pl.axis([0, 35, -1.0, 1.0])
        pl.grid(True)
        pl.savefig(iname, dpi=500)
        pl.close()

    # Saves the obtained pixel shifts in a new file.
    with open(tname, 'w') as doc:
        for i in range(len(shif[0])):
            doc.write(str(i))
            for j in range(len(shif)):
                doc.write('\t%.3f' % shif[j][i])
            doc.write('\n')

    # Inserts a column (with the spectra numbers) in the pixel shifts array.
    shif = np.insert(shif, 0, np.array(range(len(shif[0]))), axis=0)

    return np.copy(npile), np.copy(shif)


# ---------------------
# - HEADER COLLECTION -
# ---------------------
def fhea(flist, elist):
    '''
    Returns an array with the lists of the desired values for all the spectra.
     > flist: list of file's names
     > elist: list of the desired entry's names or indexes
     < temp: array with the desired lists
    '''

    # Creates an empty array that will contain the lists of the desired values
    # for all the spectra.
    temp = []
    [temp.append([]) for i in range(len(elist))]

    # For each .fits file (corresponding to a single spectrum), stores all the
    # desired values from the header in the corresponding list.
    for name in flist:
        with fits.open(name) as doc:
            for i in range(len(elist)):
                temp[i].append(cp.deepcopy(doc[0].header[elist[i]]))

    return np.copy(temp)
