#!/usr/bin/python
from __future__ import division
import numpy as np
import matplotlib.pylab as plt
import sys
import argparse

import PlottingFun as pf


# [ ] fix mtfs: DC component, smoothing, prob density, power of 2?
# [ ] figure out psfs       
# [ ] pretty print parameters
# [ ] dictionary of parameters, use to figure out legend of plots.
# [ ] move PSF and MTF into analysis.cpp

class EyePlot():
    '''
    '''

    def __init__(self):

        self._loadData()
    
    def _loadData(self, downsample=1.0):
        '''
        '''
        dat = np.genfromtxt('dat/Intensity.csv', delimiter = ',', 
            dtype="f8,f8,f8,f8,f8,f8,f8,S20,f8",
            names=True, skip_header=0)

        self._meta = {}
        self._meta['model'] = dat['model'][0]
        self._meta['iterations'] = 4 # dat['iterations'][0]
        self._meta['samples'] = int(dat.shape[0] / self._meta['iterations'])

        # Preallocate memory:
        self.Intensity = {}
        self.Intensity['rawintensity'] = np.zeros((4, self._meta['samples'])) 
        self.Intensity['intensity'] = np.zeros((4, self._meta['samples']))

        self._meta['lens_accom'] = np.zeros(self._meta['iterations'])
        self._meta['pupil_size'] = np.zeros(self._meta['iterations'])
        self._meta['offaxis'] = np.zeros(self._meta['iterations'])
        self._meta['obj_dist'] = np.zeros(self._meta['iterations'])
        self._meta['age'] = np.zeros(self._meta['iterations'])

        for i in range(0, self._meta['iterations']):
            ind1 = i * self._meta['samples']
            ind2 = (i + 1) * self._meta['samples']

            self.Intensity['rawintensity'][i,:] = dat['encircled_intensity'][ind1:ind2]
            self.Intensity['intensity'][i,:] = (dat['encircled_intensity'][ind1:ind2]
                                 / np.max(self.Intensity['rawintensity'][i,:]))

            self._meta['lens_accom'][i] = dat['lens_focus_D'][ind1]
            self._meta['pupil_size'][i] = dat['pupil_size_mm'][ind1]
            self._meta['offaxis'][i] = dat['offaxis_deg'][ind1]
            self._meta['obj_dist'][i] = dat['obj_distance_mm'][ind1]
            self._meta['age'][i] = dat['age_y'][ind1]

        # set up parameters of the eye:
        self.xvals = np.zeros((1, self._meta['samples']))
        self.xvals = dat['radius_mm'][0: self._meta['samples']]

        self._meta['retImg'] = np.max(self.xvals) # size of regino in mm 
        radians = np.tan(self._meta['retImg'] / 24)
        self._meta['mm/deg'] = rad2deg(radians)
        self._meta['eye_len'] = 24 # diameter in mm

    def _printParams(self):
        '''
        '''
        print "eye length: ", self.eye_length
        print "degrees off axis: ", self.deg_off_axis

    def _loadLSA_data(self):
        '''
        '''
        
        dat = np.genfromtxt('../dat/EyeLSA.csv', delimiter = ',')
        
        self.lsa_data = {
                    'lens_accomm': dat[:, 0], 
                    'pupil_size': dat[:, 1], 
                    'age': dat[:, 2], 
                    'accomPower': dat[:, 3], 
                    'relaxPower': dat[:, 4],
                    'defocus': dat[:, 5]
                    }

    def _genPSF(self):
        '''
        '''
        self.Intensity['PSF'] = np.zeros((4, self._meta['samples']))
        self.Intensity['PSFtotal'] = np.zeros((4, (self._meta['samples'] * 2)))
        
        # we have a integral, therefore take the deriv to get rays / bin
        deriv = np.zeros((4, self._meta['samples']))
        deriv[:, 0] = self.Intensity['rawintensity'][:, 0]
        deriv[:, 1:] = self.Intensity['rawintensity'][:,1:] - self.Intensity[
                                                    'rawintensity'][:, 0:-1]
                                    
        for i in range(0, self._meta['samples'] - 1):
            # account for increasing size of area
            radius0 = self.xvals[i]
            radius1 = self.xvals[i + 1]
            # subtract inner and outer circle area to get sliver of interest
            area = (np.pi * radius1 ** 2.0) - (np.pi * radius0 ** 2.0)
            # deriv = amount in each circle; then divide by area
            self.Intensity['PSF'][:, i] = deriv[:, i]  / area 

        # normalize so that each PSF has same integral of 1.
        for i in range(0, 4):
            self.Intensity['PSF'][i, :] = (self.Intensity['PSF'][i, :] 
                        / np.sum(self.Intensity['PSF'][i, :]))
        
        self.Intensity['PSFtotal'][:, 1:self._meta['samples'] + 1] = self.Intensity[
                                                            'PSF'][:, ::-1]
        self.Intensity['PSFtotal'][:, self._meta['samples'] + 1:] = self.Intensity[
                                                            'PSF'][:, 1:]

    def _genMTF(self):
        '''
        '''
        if 'PSF' not in self.Intensity:
            self._genPSF()

        self.Intensity['MTF'] = np.zeros((4, np.floor(self._meta['samples'] / 2.0)))
        for i in range(0, 4):

            # find next power of 2 to make more efficient.
            print self.Intensity['MTF'].shape, np.floor(self._meta['samples'] / 2)
            self.Intensity['MTF'][i, :] = np.abs(np.fft.fftshift(
                                np.fft.fft(self.Intensity['PSF'][i,:])))[
                                np.ceil(self._meta['samples'] / 2.0):]

            self.Intensity['MTF'][i, :] = (self.Intensity['MTF'][i] / 
                                        np.max(self.Intensity['MTF'][i]))

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
        for i in range(0, 4):

            ax.plot(self.xvals / self._meta['mm/deg'],
                self.Intensity['intensity'][i, :].T, 
                label = '{0} D'.format(self._meta['lens_accom'][i]))

        ax.legend(loc = 'lower right').set_title('lens accommodation')
        plt.ylim([0, 1.05])
        plt.ylabel('encircled intensity')
        plt.xlabel('degrees')
        plt.tight_layout()
        plt.show()
        
    def PSFplot(self):
        '''
        '''
        self._genPSF()

        fig = plt.figure()
        ax = fig.add_subplot(111)
        pf.AxisFormat(20)
        pf.TufteAxis(ax, ['left', 'bottom'], [5, 5], integer='on')
        size = self._meta['samples']

        for i in range(0, self._meta['iterations']):

            ax.plot(self.xvals / self._meta['mm/deg'],
                self.Intensity['PSF'][i].T, 
                label = '{0} D'.format(self._meta['lens_accom'][i]))

        ax.legend(loc = 'upper right').set_title('lens accommodation')
        #plt.ylim([-0.05, 1.01])
        plt.ylabel('point spread')
        plt.xlabel('degrees')
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

        cycles = (np.arange(1, np.ceil((self._meta['samples'] + 1) / 2))) / 2
        cpd = cycles * self._meta['mm/deg'] #/ self.deg

        for i in range(0,self._meta['iterations']):

            if not density:
                ax.plot(cpd,
                    self.Intensity['MTF'][i, :].T,
                    label = '{0} D'.format(self._meta['lens_accom'][i]))
            if density:
                ax.semilogx(cpd, 
                    decibels(self.Intensity['MTF'][i, :].T),
                    label = '{0} D'.format(self._meta['lens_accom'][i]))

        # plot diffraction limited case:
        diffract = diffraction(self._meta['mm/deg'], self._meta['samples'], self._meta['eye_len'],
                        self._meta['pupil_size'][0])
        diffract = diffract #* self._meta['mm/deg'] ** -1
        if not density: 
            ax.plot(cpd, diffract, 'k')
        if density:
            ax.semilogx(cpd, decibels(diffract), 'k')

        ax.legend(loc = 'upper right').set_title('lens accommodation')

        #plt.ylim([0, 1.0])
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

                self._findPlottingData(acc[j],'lens_accomm', age[i], 'age')
                
                self.x1dat = self.data['defocus'][self.loc]
                self.y1dat = self.data['pupil_size'][self.loc]
                
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
                
                self._findPlottingData(acc[i],'lens_accomm', 
                    age[j], 'age', 5, 'pupil_size')
                self.x1dat = np.append(self.x1dat, 
                    self.data['lens_accomm'][self.loc])
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


def diffraction(deg, samples, pupil_size, focal_len, ref_index=1.55, wavelength=550.0):
    '''See Appendix B of "Light, the Retinal Image and Photoreceptors"
    Packer & Williams.

    '''
    lam = wavelength / 1000 # convert mm into meters.
    NA = NumericalAperature(ref_index, D=pupil_size, focal_len=focal_len)

    s_0 = NA / lam # convert to radians
    s =  np.linspace(0, s_0, np.floor(samples / 2))
    #print "pupil_size: ", pupil_size, "NA: ", NA, "s_0", s_0

    dif = (2.0 / np.pi) * (np.arccos(s / s_0) - 
            (s / s_0) * np.sqrt(1.0 - (s / s_0) ** 2.0))

    return dif

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
        out = n * np.sin(np.arctan(D / 2 * focal_len))
    else:
        raise IOError("check parameters.")

    return out


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
