#! /usr/bin/env python
from __future__ import division
import numpy as np
import matplotlib.pylab as plt
import sys

from base import optics as o
from base import data as dm


def loadData(downsample=1.0):
    '''
    '''
    dat = np.genfromtxt('dat/Intensity.csv', delimiter = ',', 
        dtype="f8,f8,f8,f8,f8,f8,f8,S20,f8,f8",
        names=True, skip_header=0)

    intensity = {}
    meta = {}

    meta['model'] = dat['model'][0]
    meta['iterations'] = int(dat['iterations'][0])
    meta['samples'] = int(dat.shape[0] / meta['iterations'])

    # Preallocate memory:
    intensity['rawintensity'] = np.zeros((meta['iterations'],
                                    meta['samples'])) 
    intensity['intensity'] = np.zeros((meta['iterations'], 
                                    meta['samples']))

    meta['lens_accom'] = np.zeros(meta['iterations'])
    meta['pupil_size_%'] = np.zeros(meta['iterations'])
    meta['off_axis_*'] = np.zeros(meta['iterations'])
    meta['object_distance'] = np.zeros(meta['iterations'])
    meta['age_&'] = np.zeros(meta['iterations'])
    meta['eye_length_%'] = np.zeros(meta['iterations'])

    for i in range(0, meta['iterations']):
        ind1 = i * meta['samples']
        ind2 = (i + 1) * meta['samples']

        intensity['rawintensity'][i,:] = dat['encircled_intensity'][ind1:ind2]
        intensity['intensity'][i,:] = (dat['encircled_intensity'][ind1:ind2]
                             / np.max(intensity['rawintensity'][i,:]))

        meta['lens_accom'][i] = dat['lens_focus_D'][ind1]
        meta['pupil_size_%'][i] = dat['pupil_size_mm'][ind1]
        meta['off_axis_*'][i] = dat['offaxis_deg'][ind1]
        meta['object_distance'][i] = dat['obj_distance_mm'][ind1]
        meta['age_&'][i] = dat['age_y'][ind1]
        meta['eye_length_%'][i] = dat['axial_len_mm'][ind1]

    # set up parameters of the eye:
    xvals = dat['radius_mm'][0: meta['samples']]

    # retinal image, here is radius of image area
    meta['retImg'] = np.max(xvals) # 1/2 size of image in mm 
    radians = 2 * np.arctan(meta['retImg'] / (meta['eye_length_%'][0]))
    meta['deg'] = dm.rad2deg(radians)
    meta['mm/deg'] = meta['retImg'] * 2 / meta['deg']
    meta['xvals'] = xvals

    variable = findDependentVar(meta)
    printParams(meta)

    return intensity, meta


def findDependentVar(meta):
    '''
    '''
    variable = None

    for key in meta.iterkeys():
        try :
            if len(meta[key]) > 0:
                if key != "model":
                    if meta[key][0] != meta[key][1]:
                        variable = key
        except (IndexError, TypeError):
            pass
    return variable


def printParams(meta):
    '''
    '''
    print ' '
    for key in meta.iterkeys():
        if key != 'xvals':
            out = (key.replace('_',' ')
                .replace('%','(mm)')
                .replace('&','(y)')
                .replace('D', '(D)')
                .replace('*', '$\\degree$'))
            if len(out) < 3:
                print out, "\t\t\t\t||\t", meta[key]
            elif len(out) < 7:
                print out, "\t\t\t||\t", meta[key]
            elif len(out) < 15:
                print out, "\t\t||\t", meta[key]
            else:
                print out, " \t||\t", meta[key]
    print ' '


def findPlottingData(data, arg1, type1, arg2, 
                type2, arg3=None, type3=None):
    '''Called by LSA and power plots
    '''

    ind1 = data[type1] == arg1
    ind2 = data[type2] == arg2
    
    if arg3 == None:
        loc = ind1 & ind2 == True
    else:
        ind3 = data[type3] == arg3
        loc = ind1 & ind2 & ind3 == True

    return loc


def loadLSA_data():
    '''
    '''
    dat = np.genfromtxt('../dat/EyeLSA.csv', delimiter = ',')
    
    lsa_data = {
                'lens_accommodation_D': dat[:, 0], 
                'pupil_size_mm': dat[:, 1], 
                'age': dat[:, 2], 
                'accomPower': dat[:, 3], 
                'relaxPower': dat[:, 4],
                'defocus': dat[:, 5]
                }
    return lsa_data


def getMTF(meta, intensity):
    '''
    '''
    if 'PSF' not in intensity:
        intensity = getPSF(meta, intensity)

    intensity['MTF'] = np.zeros((meta['iterations'], 
                                    meta['samples'] ))
    for i in range(0, meta['iterations']):
        # normalize MTF
        intensity['MTF'][i, :] = o.MTF(intensity['PSFtotal'][i, :])

    return intensity


def getPSF(meta, intensity):
    '''
    '''
    intensity['PSF'] = np.zeros((meta['iterations'],
                                meta['samples']))

    intensity['PSFtotal'] = np.zeros((meta['iterations'],
                                 (meta['samples'] * 2) + 1))

    for i in range(0, meta['iterations']):
        intensity['PSF'][i, :], intensity['PSFtotal'][i, :] = (
            o.PSF(intensity['rawintensity'][i, :], meta['samples'])
        )
    return intensity
