#! /usr/bin/env python
from __future__ import division
import numpy as np
#import matplotlib
#matplotlib.use('SVG')
import matplotlib.pylab as plt
import sys
import argparse

import PlottingFun as pf


# [ ] fix mtfs: DC component?, chromatic?
# [ ] move PSF and MTF into analysis.cpp


class EyePlot():
    '''
    '''

    def __init__(self):

        self._meta = {}
        self.Intensity = {}
        self._loadData()
    
    def _loadData(self, downsample=1.0):
        '''
        '''
        dat = np.genfromtxt('dat/Intensity.csv', delimiter = ',', 
            dtype="f8,f8,f8,f8,f8,f8,f8,S20,f8,f8",
            names=True, skip_header=0)

        self._meta['model'] = dat['model'][0]
        self._meta['iterations'] = int(dat['iterations'][0])
        self._meta['samples'] = int(dat.shape[0] / self._meta['iterations'])

        # Preallocate memory:
        self.Intensity['rawintensity'] = np.zeros((self._meta['iterations'],
                                        self._meta['samples'])) 
        self.Intensity['intensity'] = np.zeros((self._meta['iterations'], 
                                        self._meta['samples']))

        self._meta['lens_accom'] = np.zeros(self._meta['iterations'])
        self._meta['pupil_size_%'] = np.zeros(self._meta['iterations'])
        self._meta['off_axis_*'] = np.zeros(self._meta['iterations'])
        self._meta['object_distance'] = np.zeros(self._meta['iterations'])
        self._meta['age_&'] = np.zeros(self._meta['iterations'])
        self._meta['eye_length_%'] = np.zeros(self._meta['iterations'])

        for i in range(0, self._meta['iterations']):
            ind1 = i * self._meta['samples']
            ind2 = (i + 1) * self._meta['samples']

            self.Intensity['rawintensity'][i,:] = dat['encircled_intensity'][ind1:ind2]
            self.Intensity['intensity'][i,:] = (dat['encircled_intensity'][ind1:ind2]
                                 / np.max(self.Intensity['rawintensity'][i,:]))

            self._meta['lens_accom'][i] = dat['lens_focus_D'][ind1]
            self._meta['pupil_size_%'][i] = dat['pupil_size_mm'][ind1]
            self._meta['off_axis_*'][i] = dat['offaxis_deg'][ind1]
            self._meta['object_distance'][i] = dat['obj_distance_mm'][ind1]
            self._meta['age_&'][i] = dat['age_y'][ind1]
            self._meta['eye_length_%'][i] = dat['axial_len_mm'][ind1]

        # set up parameters of the eye:
        self.xvals = dat['radius_mm'][0: self._meta['samples']]

        self._meta['retImg'] = np.max(self.xvals) # size of image in mm 
        radians = 2 * np.arctan(self._meta['retImg'] / (2 * 16.6))
        self._meta['deg'] = rad2deg(radians)
        self._meta['mm/deg'] = self._meta['retImg'] / self._meta['deg']

        self._findDependentVar()
        self._printParams()

    def _loadLSA_data(self):
        '''
        '''
        dat = np.genfromtxt('../dat/EyeLSA.csv', delimiter = ',')
        
        self.lsa_data = {
                    'lens_accommodation_D': dat[:, 0], 
                    'pupil_size_mm': dat[:, 1], 
                    'age': dat[:, 2], 
                    'accomPower': dat[:, 3], 
                    'relaxPower': dat[:, 4],
                    'defocus': dat[:, 5]
                    }

    def _genPSF(self):
        '''
        '''
        self.Intensity['PSF'] = np.zeros((self._meta['iterations'],
                                    self._meta['samples']))
        self.Intensity['PSFtotal'] = np.zeros((self._meta['iterations'],
                                     (self._meta['samples'] * 2) + 1))

        for i in range(0, self._meta['iterations']):

            self.Intensity['PSF'][i, :], self.Intensity['PSFtotal'][i, :] = (
                genPSF(self.Intensity['rawintensity'][i, :], self.xvals)
            )

    def _genMTF(self):
        '''
        '''
        if 'PSF' not in self.Intensity:
            self._genPSF()

        self.Intensity['MTF'] = np.zeros((self._meta['iterations'], 
                                        self._meta['samples'] + 1))
        for i in range(0, self._meta['iterations']):
            # normalize MTF
            self.Intensity['MTF'][i, :] = genMTF(self.Intensity['PSFtotal'][i, :])

    def _findDependentVar(self):
        '''
        '''
        self.variable = None
        for key in self._meta.iterkeys():
            try :
                if len(self._meta[key]) > 0:
                    if key != "model":
                        if self._meta[key][0] != self._meta[key][1]:
                            self.variable = key
            except (IndexError, TypeError):
                pass
        if self.variable == None:
            raise IOError("No dependent variable found")

    def _printParams(self):
        '''
        '''
        print ' '
        for key in self._meta.iterkeys():
            out = (key.replace('_',' ')
                .replace('%','(mm)')
                .replace('&','(y)')
                .replace('D', '(D)')
                .replace('*', '$\\degree$'))
            if len(out) < 3:
                print out, "\t\t\t\t||\t", self._meta[key]
            elif len(out) < 7:
                print out, "\t\t\t||\t", self._meta[key]
            elif len(out) < 15:
                print out, "\t\t||\t", self._meta[key]
            else:
                print out, " \t||\t", self._meta[key]
        print ' '

    def _findPlottingData(self, arg1, type1, arg2, 
                    type2, arg3=None, type3=None):
        '''Called by LSA and power plots
        '''

        ind1 = self.data[type1] == arg1
        ind2 = self.data[type2] == arg2
        
        if arg3 == None:
            self.loc = ind1 & ind2 == True
        else:
            ind3 = self.data[type3] == arg3
            self.loc = ind1 & ind2 & ind3 == True
    
    def EncircledIntensityPlot(self):
        '''
        '''

        fig = plt.figure()
        ax = fig.add_subplot(111)
        pf.AxisFormat(20)
        pf.TufteAxis(ax, ['left', 'bottom'], [5, 5], integer='on')

        print (self._meta['samples'], len(self.xvals), 
            self.Intensity['intensity'][0,:].shape)
        for i in range(0, self._meta['samples']):

            ax.plot(self.xvals / self._meta['mm/deg'],
                self.Intensity['intensity'][i, :].T, 
                label = '{0} D'.format(self._meta['lens_accom'][i]))

        ax.legend(loc = 'lower right').set_title('lens accommodation')
        plt.ylim([0, 1.05])
        plt.ylabel('encircled intensity')
        plt.xlabel('degrees')
        plt.tight_layout()
        plt.show()
        
    def PSFplot(self, symmetric=True):
        '''
        '''
        self._genPSF()

        fig = plt.figure()
        ax = fig.add_subplot(111)
        pf.AxisFormat(20)
        pf.TufteAxis(ax, ['left', 'bottom'], [5, 5], integer='on')
        size = self._meta['samples']
        if symmetric:
            _x = np.linspace(-self._meta['retImg'], self._meta['retImg'],
                (self._meta['samples'] * 2) + 1)

        for i in range(0, self._meta['iterations']):

            if symmetric:
                ax.plot(_x / self._meta['mm/deg'] * 60,
                    self.Intensity['PSFtotal'][i].T, 
                    label = '{0}'.format(self._meta[self.variable][i]))
            else:
                ax.plot(self.xvals / self._meta['mm/deg'] * 60,
                    self.Intensity['PSF'][i].T, 
                    label = '{0}'.format(self._meta[self.variable][i]))

        ax.legend().set_title(self.variable.replace('_',' ')
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
        plt.tight_layout()
        plt.show()
        
    def MTFplot(self, density=False):
        '''
        '''
        self._genMTF()

        fig = plt.figure()
        ax = fig.add_subplot(111)
        pf.AxisFormat(20)
        pf.TufteAxis(ax, ['left', 'bottom'], [5, 5], integer='on')

        cycles = (np.arange(0, self.Intensity['MTF'].shape[1] )) / 2
        cpd = cycles / self._meta['retImg'] * self._meta['mm/deg'] 

        for i in range(0,self._meta['iterations']):

            if not density:
                ax.plot(cpd,
                    self.Intensity['MTF'][i, :].T,
                    label = '{0}'.format(self._meta[self.variable][i]))
            if density:
                ax.semilogx(cpd, 
                    decibels(self.Intensity['MTF'][i, :].T),
                    label = '{0}'.format(self._meta[self.variable][i]))

        # plot diffraction limited case:
        diffract, _x = diffraction(self._meta['mm/deg'], 
                        self.Intensity['MTF'].shape[1], 
                        self._meta['pupil_size_%'][0],
                        16.6)#self._meta['eye_length_%'][0])

        if not density: 
            ax.plot(_x * (self._meta['samples'] / 2) / self._meta['deg'],
                diffract, 'k')
        if density:
            ax.semilogx(cpd, decibels(diffract), 'k')

        ax.legend().set_title(self.variable.replace('_', ' ')
            .replace('%', '(mm)')
            .replace('&',' (y)')
            .replace('D', '(D)')
            .replace('*', '$\\degree$'))

        plt.ylim([0, 1.0])
        plt.xlim([0, 60.0])
        plt.ylabel('modulation transfer')
        plt.xlabel('cycles / deg')
        plt.tight_layout()
        plt.show()
        
    def LSAplot(self, acc = [0, 6], age = [13, 20] ):
        '''This function is not correct.
        '''
        
        self._loadLSA_data()

        fig = plt.figure()

        ax = fig.add_subplot(111)
        pf.AxisFormat(20)
        pf.TufteAxis(ax, ['left', 'bottom'], [5,5], integer='on')
        
        for j in range(0,len(age)):
            for i in range(0,len(acc)):

                self._findPlottingData(acc[j],'lens_accommodation_D', age[i], 'age_y')
                
                self.x1dat = self.data['defocus'][self.loc]
                self.y1dat = self.data['pupil_size_mm'][self.loc]
                
                ax.plot(self.x1dat,self.y1dat, linewidth=2.5,
                        label='{0} D, {1} y'.format(acc[i], age[j]));
                
        ax.legend(loc = 'lower right').set_title('Accommodation')

        plt.ylabel('pupil size')
        plt.xlabel('diopters')
        plt.tight_layout()
        plt.show()

    def PowerPlots(self, acc = [0,2,4,6,8], age = [10, 20]):
        '''
        '''
        self._loadLSA_data()

        fig = plt.figure()
        
        ax = fig.add_subplot(111)
        pf.AxisFormat(20)
        pf.TufteAxis(ax, ['left', 'bottom'], [5,5], integer='on')
        linestyle = ['ko-','ko--','ko-.','koD']
        for j in range(0, len(age)):
            self.x1dat, self.y1dat, self.topaxis = [],[],[]
            for i in range(0, len(acc)):
                
                self._findPlottingData(acc[i],'lens_accommodation_D', 
                    age[j], 'age', 5, 'pupil_size_mm')
                self.x1dat = np.append(self.x1dat, 
                    self.data['lens_accommodation_D'][self.loc])
                self.y1dat = np.append(self.y1dat, 
                    self.data['accomPower'][self.loc])
                self.topaxis = np.append(self.topaxis, 
                    self.data['defocus'][self.loc])
            
            
            ax.plot(self.x1dat,self.y1dat, linestyle[j],
                linewidth=2.5, ms=10, 
                label = '{0} years'.format(age[j]))

        ax.legend(loc = 'upper left')
        
        plt.xlim([min(acc)-0.1, max(acc)+0.1])
        plt.ylabel('power (D)')
        plt.xlabel('accommodation (D)')
        plt.tight_layout()
        plt.show()

def genPSF(intensity, xvals):
    '''

    '''
    samples = len(xvals)
    PSF = np.zeros(samples)
    PSFtotal = np.zeros((samples * 2) + 1)

    # we have an integral, therefore take the deriv to get rays / bin
    deriv = np.zeros((samples))
    deriv[0] = intensity[0]
    deriv[1:] = intensity[1:] - intensity[0:-1]

    for i in range(0, samples - 1):

        # account for increasing size of area
        radius0 = xvals[i]
        radius1 = xvals[i + 1]

        # subtract inner and outer circle area to get sliver of interest
        area = (np.pi * radius1 ** 2.0) - (np.pi * radius0 ** 2.0)

        # deriv = amount in each circle; then divide by area
        PSF[i] = deriv[i]  / area 

    # normalize so that each PSF has same integral of 1.

    PSF = PSF / np.sum(PSF)

    PSFtotal[1:samples + 1] = PSF[::-1]
    PSFtotal[samples:-1] = PSF

    return PSF, PSFtotal


def genMTF(PSFtotal):
    '''
    '''
    samples = len(PSFtotal)
    MTF = np.zeros(samples / 2)

    # make sure PSF is normalized
    normPSF = PSFtotal / np.sum(PSFtotal)

    # do the FFT, take only right half
    temp = np.abs(np.fft.fftshift(np.fft.fft(normPSF)))
    temp = temp[np.floor(samples / 2):]

    # make sure we only get real part
    MTF = np.real(temp)

    return MTF


def diffraction(deg, samples, pupil_size_mm, focal_len, ref_index=1.336, wavelength=550.0):
    '''See Appendix B of "Light, the Retinal Image and Photoreceptors"
    Packer & Williams.

    '''
    lam = wavelength / 1000 # convert mm into meters.
    NA = NumericalAperature(ref_index, D=pupil_size_mm, focal_len=focal_len)

    s_0 = NA / lam # convert to radians
    s =  np.linspace(0, s_0, samples)
    print "NA: ", NA, "s_0", s_0

    dif = (2.0 / np.pi) * (np.arccos(s / s_0) - 
            (s / s_0) * np.sqrt(1.0 - (s / s_0) ** 2.0))

    return dif, s

def NumericalAperature(n, theta=None, D=None, focal_len=None):
    '''
    Find the numerical aperature of a system

    :param n: refractive index
    :param theta: angle of marginal rays

    According to the formula

    $$NA = n \\sin(\\theta)$$

    or 

    $$NA = n \\sin(\\arctan(\\frac{D}{2f}))$$
       
    '''
    if D is None and focal_len is None and theta is not None:
        out = n * np.sin(theta)
    elif theta is None and D is not None and focal_len is not None:
        out = n * np.sin(np.arctan(D / (2 * focal_len)))
    else:
        raise IOError("check parameters.")

    return out

def nextpow2(n):
    '''
    '''
    m_f = np.log2(n)
    m_i = np.ceil(m_f)

    return 2 ** m_i


def rad2deg(radians):
    '''Convert radians to degrees.
    '''
    return radians * 180.0 / np.pi


def decibels(x):
    '''
    '''
    return 10.0 * np.log10(x)


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
    parser.add_argument("-v", "--verbose", action="store_true", 
                         help="display all plots")    
    args = parser.parse_args()

    out = EyePlot()

    if args.lsa or args.verbose:
        out.LSAplot(acc = [2,6], age = [10, 15])
    if args.power or args.verbose:
        out.PowerPlots()
    if args.intensity or args.verbose:
        out.EncircledIntensityPlot()
    if args.psf or args.verbose:
        out.PSFplot()
    if args.mtf or args.verbose:
        out.MTFplot()
