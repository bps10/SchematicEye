#!/usr/bin/python
from __future__ import division
import numpy as np
import matplotlib.pylab as plt
import sys

import PlottingFun as pf
       

# [ ] for self.xvals, start using self.intensity['radius'], 
#       but correct for mm_per_deg 
# [ ] figure out s_0 in diffraction limited case
# [ ] use fabric to send command line arguments - or use make to ensure up to
#       date?
# [ ] fix mtfs: DC component, smoothing, prob density, power of 2?
# [ ] figure out psfs       

class EyePlot():

    def __init__(self):

        self.loadData()
        
        self.eye_length = 24 # diameter in mm
        self.deg = 1 # in deg
        mm_per_deg = (self.eye_length  * np.pi) / 360.0
        self.xvals = np.linspace(0, 1, self.samples) / mm_per_deg
        #* (self.samples / self.samples) #* self.Intensity['downsample']
        
    ## Data Manipulation Methods ##
    ###############################
    
    def loadData(self, downsample=1.0):
        '''
        '''
        dat = np.genfromtxt('dat/Intensity.csv', delimiter = ',',
            skip_header=1)
        self.Intensity = {}
        self.Intensity['downsample'] = 1.0 # default = 1.0 or none.
        dat = dat[::int(self.Intensity['downsample'])]

        self.samples = int( dat.shape[0] / 4.0 )

        # Preallocate memory:
        self.Intensity['rawintensity'] = np.zeros((4, self.samples)) 
        self.Intensity['intensity'] = np.zeros((4, self.samples))
        self.Intensity['radius'] = np.zeros((4, self.samples))
        self.Intensity['lens_accom'] = np.zeros((4, self.samples))
        self.Intensity['pupil_size'] = np.zeros((4, self.samples))
        for i in range(0,4):

            totals = np.max(dat[i * self.samples:(i + 1) * self.samples, 0])
            self.Intensity['rawintensity'][i,:] = dat[
                            i * self.samples:(i + 1) * self.samples, 0]
            self.Intensity['intensity'][i,:] = dat[
                        i * self.samples:(i + 1) * self.samples, 0] / totals
            self.Intensity['radius'][i,:] = dat[
                        i * self.samples:(i + 1) * self.samples, 1]
            self.Intensity['lens_accom'][i,:] = dat[
                        i * self.samples:(i + 1) * self.samples, 2]
            self.Intensity['pupil_size'][i,:] = dat[
                        i * self.samples:(i + 1) * self.samples, 3]

    def loadLSA_data(self):
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
        self.Intensity['PSF'] = np.zeros((4, self.samples))
        self.Intensity['PSFtotal'] = np.zeros((4, (self.samples * 2)))
           
        deriv = np.zeros((4, self.samples))
        deriv[:, 0]     = self.Intensity['intensity'][:,0]
        deriv[:, 1:] = self.Intensity['intensity'][:,1:] - self.Intensity[
                                                    'intensity'][:, 0:-1]
                                    
                                    
        self.Intensity['PSF'][:, 0] = 1.0
        for i in range(1,self.samples-1):
            self.Intensity['PSF'][:, i] = (self.Intensity['PSF'][:, i - 1] - 
                                            deriv[:, i]) ** 1.0
        
        self.Intensity['PSFtotal'][:, 1:self.samples + 1] = self.Intensity[
                                                            'PSF'][:, ::-1]
        self.Intensity['PSFtotal'][:, self.samples + 1:] = self.Intensity[
                                                            'PSF'][:, 1:]

    def _genMTF(self):
        '''
        '''
        if 'PSF' not in self.Intensity:
            self._genPSF()

        self.Intensity['MTF'] = np.zeros((4, int(self.samples / 2.0) ))
        for i in range(0, 4):

            # find next power of 2 to make more efficient.

            self.Intensity['MTF'][i, ] = np.abs(np.fft.fftshift(
            np.fft.fft(self.Intensity['PSF'][i,:])))[
                                int(self.samples / 2.0):]

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
    
    ## Plotting Methods ##
    ######################
        
    def EncircledIntensityPlot(self):
        '''
        '''

        fig = plt.figure()
        ax = fig.add_subplot(111)
        pf.AxisFormat(20)
        pf.TufteAxis(ax, ['left', 'bottom'], [5, 5], integer='on')

        print (self.samples, len(self.xvals), 
            self.Intensity['intensity'][0,:].shape)
        for i in range(0, 4):

            ax.plot(self.xvals,
                self.Intensity['intensity'][i, :].T, 
                label = '{0} D'.format(self.Intensity['lens_accom'][i, 0]))

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
        size = self.Intensity['PSFtotal'].shape[1] / 2.0

        for i in range(0, self.Intensity['PSFtotal'].shape[0]):

            ax.plot(self.xvals,
                self.Intensity['PSFtotal'][i, size:size + self.samples].T, 
                label = '{0} D'.format(self.Intensity['lens_accom'][i, 0]))

        ax.legend(loc = 'upper right').set_title('lens accommodation')
        plt.ylim([-0.05, 1.01])
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

        cycles = (np.arange(1, (self.samples + 1) / 2)) / 2
        cpd = cycles / self.deg

        for i in range(0,self.Intensity['MTF'].shape[0]):

            if not density:
                ax.semilogy(cpd,
                    self.Intensity['MTF'][i, :].T,
                    label = '{0} D'.format(self.Intensity['lens_accom'][i,0]))
            if density:
                ax.semilogx(cpd, 
                    decibels(self.Intensity['MTF'][i, :].T),
                    label = '{0} D'.format(self.Intensity['lens_accom'][i,0]))

        # plot diffraction limited case:
        diffract = diffraction(cpd, 3) #self.Intensity['pupil_size'][0,0])
        if not density: 
            ax.semilogy(cpd, diffract, 'k')
        if density:
            ax.semilogx(cpd, decibels(diffract), 'k')

        ax.legend(loc = 'upper right').set_title('lens accommodation')

        #plt.ylim([0, 1.0])
        plt.ylabel('modulation transfer')
        plt.xlabel('cycles / deg')
        plt.tight_layout()
        plt.show()
        
    def LSAplot(self, acc = [0, 6], age = [13, 20] ):
        '''This function is not correct.
        '''
        
        self.loadLSA_data()

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
        self.loadLSA_data()

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


def diffraction(spatial_freq, pupil_size, wavelength=550.0):
    '''See Appendix B of "Light, the Retinal Image and Photoreceptors"
    Packer & Williams.

    '''
    lam = wavelength / 1000 # convert mm into meters.
    NA = NumericalAperature(1.5, D=pupil_size, focal_len=24.0)
    s = np.linspace(0, 1, len(spatial_freq))

    s_0 = NA / lam # * (180 / np.pi)# convert to radians
    print "pupil_size: ", pupil_size, "NA: ", NA, "s_0", s_0
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


def decibels(x):
    '''
    '''
    return 10.0 * np.log10(x)


if __name__ == "__main__":

    # main function can change options
        
    out = EyePlot()
    #out.LSAplot(acc = [2,6], age = [10, 15])
    #out.PowerPlots()
    out.EncircledIntensityPlot()
    out.PSFplot()
    out.MTFplot()
