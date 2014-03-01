#! /usr/bin/env python
from __future__ import division
import numpy as np
import matplotlib.pylab as plt
import argparse

from base import plot as pf
from base import optics as o

import eye_analysis as ea 
import eye as eye

def EncircledIntensityPlot(meta, intensity):
    '''
    '''
    fig = plt.figure()
    fig.set_tight_layout(True)
    ax = fig.add_subplot(111)
    pf.AxisFormat(20)
    pf.TufteAxis(ax, ['left', 'bottom'], [5, 5], integer='on')

    for i in range(0, meta['iterations']):
        ax.plot(meta['xvals'] / meta['mm/deg'], intensity['intensity'][i, :].T, 
            label='{0} D'.format(meta['lens_accom'][i]))

    ax.legend(loc='lower right').set_title('lens accommodation')
    plt.ylim([0, 1.05])
    plt.ylabel('encircled intensity')
    plt.xlabel('degrees')
    plt.show()

    
def PSFplot(meta, intensity, symmetric=True):
    '''
    '''
    intensity = ea.getPSF(meta, intensity)
    variable = ea.findDependentVar(meta)

    fig = plt.figure()
    fig.set_tight_layout(True)
    ax = fig.add_subplot(111)
    pf.AxisFormat(20)
    pf.TufteAxis(ax, ['left', 'bottom'], [5, 5], integer='on')
    size = meta['samples']

    if symmetric:
        _x = np.linspace(-meta['retImg'] / meta['mm/deg'], 
            meta['retImg'] / meta['mm/deg'],
            (meta['samples'] * 2) + 1)

    for i in range(0, meta['iterations']):

        if variable is not None:
            ax.plot(_x * 60,
                intensity['PSFtotal'][i].T, 
                label = '{0}'.format(meta[variable][i]))
        else:
            ax.plot(_x * 60, intensity['PSFtotal'][i].T)

    if variable is not None:
        ax.legend().set_title(variable.replace('_',' ')
            .replace('%','(mm)')
            .replace('&','(y)')
            .replace('D', '(D)')
            .replace('*', '$\\degree$'))

    if symmetric:
        plt.xlim([-3, 3])
    else:
        plt.xlim([0, 2])
    plt.ylabel('point spread')
    plt.xlabel('arc min')
    plt.show()
    

def MTFplot(meta, intensity, density=False):
    '''
    '''
    intensity = ea.getMTF(meta, intensity)
    variable = ea.findDependentVar(meta)

    fig = plt.figure()
    fig.set_tight_layout(True)
    ax = fig.add_subplot(111)
    pf.AxisFormat(20)
    pf.TufteAxis(ax, ['left', 'bottom'], [5, 5], integer='on')

    cycles = (np.arange(0, intensity['MTF'].shape[1] )) / 2
    cpd = cycles / meta['retImg'] * meta['mm/deg'] 

    for i in range(0, meta['iterations']):

        if variable is not None:
            ax.plot(cpd,
                intensity['MTF'][i, :].T,
                label = '{0}'.format(meta[variable][i]))
        else:
            ax.plot(cpd, intensity['MTF'][i, :].T)

    # plot diffraction limited case:
    diffract, _x = o.diffraction(intensity['MTF'].shape[1], 
                    meta['pupil_size_%'][0], 16.6)
    if not density: 
        ax.plot(_x, diffract, 'k')
    if density:
        ax.semilogx(cpd, decibels(diffract), 'k')

    if variable is not None:
        ax.legend().set_title(variable.replace('_',' ')
            .replace('%','(mm)')
            .replace('&','(y)')
            .replace('D', '(D)')
            .replace('*', '$\\degree$'))

    plt.ylim([0, 1.0])
    plt.xlim([0, 60.0])
    plt.ylabel('modulation transfer')
    plt.xlabel('cycles / deg')
    plt.show()

    
def LSAplot(acc = [0, 6], age = [13, 20] ):
    '''This function is not correct.
    '''
    
    lsa_data = ea.loadLSA_data()

    fig = plt.figure()
    fig.set_tight_layout(True)
    ax = fig.add_subplot(111)
    pf.AxisFormat(20)
    pf.TufteAxis(ax, ['left', 'bottom'], [5,5], integer='on')
    
    for j in range(0,len(age)):
        for i in range(0,len(acc)):

            ea.findPlottingData(acc[j],'lens_accommodation_D', age[i], 'age_y')
            
            x1dat = data['defocus'][loc]
            y1dat = data['pupil_size_mm'][loc]
            
            ax.plot(x1dat, y1dat, linewidth=2.5,
                    label='{0} D, {1} y'.format(acc[i], age[j]));
            
    ax.legend(loc = 'lower right').set_title('Accommodation')

    plt.ylabel('pupil size')
    plt.xlabel('diopters')
    plt.show()


def PowerPlots(acc = [0,2,4,6,8], age = [10, 20]):
    '''
    '''
    lsa_data = ea.loadLSA_data()

    fig = plt.figure()
    fig.set_tight_layout(True)
    ax = fig.add_subplot(111)
    pf.AxisFormat(20)
    pf.TufteAxis(ax, ['left', 'bottom'], [5,5], integer='on')

    linestyle = ['ko-','ko--','ko-.','koD']
    for j in range(0, len(age)):
        x1dat, y1dat, topaxis = [],[],[]
        for i in range(0, len(acc)):
            
            ea.findPlottingData(acc[i],'lens_accommodation_D', 
                age[j], 'age', 5, 'pupil_size_mm')
            x1dat = np.append(x1dat, 
                data['lens_accommodation_D'][loc])
            y1dat = np.append(y1dat, 
                data['accomPower'][loc])
            topaxis = np.append(topaxis, 
                data['defocus'][loc])
        
        
        ax.plot(x1dat, y1dat, linestyle[j],
            linewidth=2.5, ms=10, 
            label = '{0} years'.format(age[j]))

    ax.legend(loc = 'upper left')
    
    plt.xlim([min(acc)-0.1, max(acc)+0.1])
    plt.ylabel('power (D)')
    plt.xlabel('accommodation (D)')
    plt.show()


def plotComparisonMTF(save_plots=False, legend=False):
    """
    Plot peripheral MTF with a comparison to Navarro et al 1993 or 
    Williams et al. 1996
    
    :param save_plots: decide whether to save plots (True) or not (False).
    :type save_plots: bool
    :param legend: turn legend on (True) or off (False). Default = True. \n
    :type legend: bool
    
    Currently supports 0, 10, 20, 40 degrees eccentricity.
    
    **This produces:**
    
    .. figure:: ../../Figures/MTFperiphery.png
       :height: 300px
       :width: 400px
       :align: center        
       
       **Fig 1:** A family of MTF curves from experimental data (dotted) 
       and schematic eye.        
    """   
    samples = 399
    cycles = np.arange(0, samples) / 2
    cpd = cycles / 0.1995 * 0.417492437017 

    #Ray trace params:
    wavelength = 632.8
    pupil = 3

    PAPER = False
    if PAPER:
        paper = 'Williams1996_clc' 
    else:
        paper = 'Navarro1993' 

    Fovea = o.MTF_analytical(cpd, 0, paper)
    TenDeg = o.MTF_analytical(cpd, 10, paper)
    TwentyDeg = o.MTF_analytical(cpd, 20, paper)
    
    fig = plt.figure(figsize=(8,6))
    fig.set_tight_layout(True)
    ax = fig.add_subplot(111)

    pf.AxisFormat(TickSize=15, FONTSIZE=26)
    pf.TufteAxis(ax, ['left', 'bottom'], [5,5])
    
    #Navarro et al 1993 analytical func:
    ax.plot(cpd, Fovea, 'm--')
    ax.plot(cpd, TenDeg, 'r--')
    ax.plot(cpd, TwentyDeg, 'g--')
    
    intensity = traceEye(1e8, 0, pupil, 0, wavelength)
    psf = o.PSF(intensity, samples)[1]
    mtf = o.MTF(psf)
    ax.plot(cpd, mtf, 'm-', label='fovea')

    intensity = traceEye(1e8, 10, pupil, 0, wavelength)
    psf = o.PSF(intensity, samples)[1]
    mtf = o.MTF(psf)
    ax.plot(cpd, mtf, 'r-', label='10 deg')  

    intensity = traceEye(1e8, 20, pupil, 0, wavelength)
    psf = o.PSF(intensity, samples)[1]
    mtf = o.MTF(psf)    
    ax.plot(cpd, mtf, 'g-', label='20 deg')
    
    ax.legend(loc='upper right')#,title='object dist, retinal location')

    plt.ylim([0, 1.0])
    plt.xlim([0, 60.0])
    
    plt.xlabel('spatial frequency (cycles / deg)')
    plt.ylabel('modulation transfer')
    
    if save_plots:
        fig.savefig('../img/MTFperiphery.png')

    plt.show()


def traceEye(object_distance=1e8, off_axis=0, pupil_size=3, diopters=0, wavelength=550):
    '''
    '''
    intensity = eye.py_eye(object_distance, off_axis, pupil_size, diopters, wavelength)

    return intensity


def main(args):
    '''
    '''
    intensity, meta = ea.loadData()

    if args.lsa or args.verbose:
        LSAplot(acc = [2,6], age = [10, 15])

    if args.power or args.verbose:
        PowerPlots()

    if args.intensity or args.verbose:
        EncircledIntensityPlot(meta, intensity)

    if args.psf or args.verbose:
        PSFplot(meta, intensity)

    if args.mtf or args.verbose:
        MTFplot(meta, intensity)

    if args.comp or args.verbose:
        plotComparisonMTF()


if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    
    parser.add_argument("-l", "--lsa", action="store_true", 
                        help="display LSA plot (not working)")
    parser.add_argument("-d", "--power", action="store_true",
                        help="display optical power plot")
    parser.add_argument("-i", "--intensity", action="store_true", 
                         help="display encircled intensity plot")
    parser.add_argument("-p", "--psf", action="store_true",
                        help="display PSF plot")
    parser.add_argument("-m", "--mtf", action="store_true", 
                         help="display MTF plot")
    parser.add_argument("-c", "--comp", action="store_true", 
                         help="display comparison plot")

    parser.add_argument("-v", "--verbose", action="store_true", 
                         help="display all plots")    
    args = parser.parse_args()

    main(args)
