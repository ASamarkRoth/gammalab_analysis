#!/usr/bin/env python3
import csv                               ## for reading in our data files
import logging                           ## for orderly print output
import numpy as np                       ## for handling of data
import matplotlib.pyplot as plt          ## for plotting
from scipy.optimize import curve_fit     ## to fit functions to the data
from scipy import signal                 ## to find maxima in the data
from scipy import stats                  ## for linear regressions
import sys                               ## useful system calls (used to exit cleanly)

class Spectrum:
    """A class to hold our spectrum measurement data and meta data (such as duration)"""
    def __init__(self, filename):
        self.filename = filename
        self.x = np.array(np.zeros(1))    ## creates a new empty array; will later store our x values
        self.y = np.array(np.zeros(1))    ## creates a new empty array; will later store our y values
        self.name = filename   ## a more descriptive name, can be used e.g. in legends
        self.duration = 0
    def subtract(self, m):
        self.y = self.y - m.y
        ## these are spectra: cannot have counts below 0, so remove these here and set them to 0 instead:
        self.y[self.y < 0] = 0;
    def scale(self, scale):
        self.y *= scale
    def calibrate(self, slope, intercept):
        self.x = self.x*slope + intercept
        
def read_mca_data_file(filename):
    """Reads in a data file (csv format) stored by the Maestro MCA software and returns a 'Spectrum' object. Tested with Maestro Version 6.05 """
    log = logging.getLogger('betalab_analysis') ## set up logging
    m = Spectrum(filename) ## create a new Spectrum measurement object; this is what we return in the end
    log.info("Reading data from file '" + filename + "'")
    try:
        with open(filename, newline='') as f:
            reader = csv.reader(f) ## use the python csv module to parse the file
            interval = []          ## start/stop channel numbers used to assign correct x values to the data points
            ## first parse the "header" of the data file (until the '$DATA:' line) containing all the meta data
            for row in reader:
                if row[0] == '$MEAS_TIM:':
                    ## this item gives the duration of the measurement
                    log.debug("Parsing MEAS_TIM header info")
                    row = next(reader)
                    duration = [int(s) for s in row[0].split(' ')]
                    m.duration = duration[1] ## two parts: real time/live time; take the second
                if row[0] == '$DATA:':
                    ## this is the last part of the header and contains the start/stop channel numbers
                    log.debug("Parsing DATA header info")
                    row = next(reader)
                    interval = [int(s) for s in row[0].split(' ')]
                    ## "DATA" is the last item: stop with the header processing
                    break
            ## TODO: make sure that the file does not end before we have parsed the header!
            log.debug("Done with header parsing")
            nbins = int(interval[1]-interval[0])+1
            m.y = np.array(np.zeros(nbins))
            ## continue, now reading data
            for idx, row in enumerate(reader):
                if idx >= nbins:
                    break
                m.y[idx] = int(row[0])
            m.x = np.arange(interval[0], interval[1]+1,1)
            log.debug("Loaded all data from file")
    except IOError:
        log.error("Could not find the file '"+str(filename)+"'")
        sys.exit(-1)
    return m


def gaussfcn(x, *p):
    """ gauss function to be used for fits to the data"""
    A, mu, sigma = p
    return A*np.exp(-(x-mu)**2/(2.*sigma**2))

class Gauss:
    """A class to hold coefficients for Gaussian distributions"""
    def __init__(self, A, mu, sigma):
        self.A = A
        self.mu = mu
        self.sigma = sigma
    def value(self, x):
        return gaussfcn(x, self.A, self.mu, self.sigma)

def line(x, *p):
    """ straight line function to be used for fits to the data"""
    a, b = p
    return a*x+b

def fit_gaussians_to_measurement(m):
    """ fits all gaussians in a spectrum measurement and returns a list of coefficients"""
    ## list to store the paramters of the fitted gaussians in
    gaussians = []
    
    ## find peaks in m.y with range of widths given by an array (range from X1 to X2 in steps of X3)
    peakind = signal.find_peaks_cwt(m.y, np.arange(10,80,5)) 

    for p in peakind:
        log.info("Found peak in the data at position x = " + str(m.x[p]) + "; fitting it with a gaussian")
        ## p0 is the initial guess for the fitting coefficients (A, mu and sigma above)
        p0 = [1., m.x[p], 1.] ## mu is given by one of the found peaks positions
        ## scypi gives a warning if the fit does not work; we want to know about those, so we set them up to be caught here:
        import warnings
        from scipy.optimize import OptimizeWarning
        with warnings.catch_warnings():
            warnings.simplefilter("error", OptimizeWarning)
            ## perform the gaussian fit to the data:
            try:
                ## use the scipy curve_fit routine (uses non-linear least squares to perform the fit)
                ## see http://docs.scipy.org/doc/scipy-0.16.1/reference/generated/scipy.optimize.curve_fit.html
                coeff, var_matrix = curve_fit(gaussfcn, m.x[p-10:p+10], m.y[p-10:p+10], p0=p0) ## fit using "slices" of the arrays with +/- 10 around peak
                ## create a Gauss object with the fitted coefficients for better code readability
                g = Gauss(*coeff)
            except (RuntimeError, OptimizeWarning, TypeError):
                ## the minimization did not work out... log it and continue to next peak
                log.info("  - gaussian fit failed!")
                continue
        ## filter the results -- or we get a lot of "spurious" peaks
        xdynrange = m.x.shape[0] ## dynamic range in x
        if g.sigma > 0.1*xdynrange: ## check width of gaussian in percent of the dynamic range
            log.info("  - sigma out of bounds: " + str(g.sigma))
            continue
        if g.mu > xdynrange-0.1*xdynrange or g.mu < 6: ## check center position: should not be at the limits of measurement range
            log.info("  - mu out of bounds: " + str(g.mu))
            continue
        if g.A < 0.25*np.average(m.y): ## check if the peak hight is at least 25% over the average data value
            log.info("  - A out of bounds: " + str(g.A))
            continue            
        log.info("  - fit result: A = " + str(g.A) + ", mu = " + str(g.mu) + ", sigma = " + str(g.sigma) + ". ")
        ## store the results
        gaussians.append(g)
    return gaussians

##                  _
##  _ __ ___   __ _(_)_ __    _ __  _ __ ___   __ _ _ __ __ _ _ __ ___
## | '_ ` _ \ / _` | | '_ \  | '_ \| '__/ _ \ / _` | '__/ _` | '_ ` _ \
## | | | | | | (_| | | | | | | |_) | | | (_) | (_| | | | (_| | | | | | |
## |_| |_| |_|\__,_|_|_| |_| | .__/|_|  \___/ \__, |_|  \__,_|_| |_| |_|
##                           |_|              |___/
if __name__ == '__main__':
    ## set up some print-out routines (logging)
    FORMAT = '%(asctime)s %(name)s:line %(lineno)-4d %(levelname)-8s %(message)s'
    logging.basicConfig(format=FORMAT)
    log = logging.getLogger('betalab_analysis') ## set up logging
    log.setLevel("INFO")

    ## ____            ___   ___       __ __   __    ___   ___
    ##/ ___| _ __     / _ \ / _ \     / / \ \ / /   / _ \ / _ \
    ##\___ \| '__|___| (_) | | | |   / /   \ V /___| (_) | | | |
    ## ___) | | |_____\__, | |_| |  / /     | |_____\__, | |_| |
    ##|____/|_|         /_/ \___/  /_/      |_|       /_/ \___/

    ## setup the (first) plot
    plt.xlabel('channel number')
    plt.ylabel('counts')
    plt.title("Sr-90/Y-90 MCA spectrum")
    plt.xlim([0, 512])                                    ## to set the x-axis range ([xmin, xmax]), use ylim() to set y axis limits
    plt.annotate('Now all you need is data! :)',          # text to put there
                xy=(0.35, 0.5),                           # coordinates for text
                xycoords='figure fraction')               # set relative coordinates
    ##plt.yscale('log')                                     ## set y axis to log scale
    plt.grid(True)                                        ## enable a grid to guide the eye

    ## Delete this to continue!
    plt.show()           ## <-- shows the plot
    log.info("Stopping analysis here... modify code to continue! ")    
    sys.exit() ## quit for now...
    ## ... delete until here :)
    
    ## read in the sr90 measurement file
    sr90 = read_mca_data_file('data/sr90.Spe')
    if not sr90:
        ## looks like we couldn't open the file, so just exit here
        sys.exit()
    ## plot the sr90 raw measurement
    plt.plot(sr90.x, sr90.y, marker='o', label="Sr-90/Y-90 raw data")  ## marker='o' parameter: plot with markers

    ## Side quest: What about backgrounds? One should look at this...
    ## -- load a file with a (long!) background measurement
    ## => read in the background measurement file here using 'read_mca_data_file'
    ## -- make sure that the background is correctly normalized to the Sr-90 data by taking measurement times into account
    ## => use the "duration" property of our spectra to calculate the factor and use "scale" method to scale the _background_
    ## -- subtract the background from the measurement data
    ## => use the 'subtract' method of the Spectrum class 
    ## -- plot both the background and the data w/o background into the same plot
    ## => how much background is there? where is it located? will it affect the measurement precision?

    plt.legend()     ## generate the legend (with the "label" information from the plots)

    ## Delete this to continue!
    plt.show()           ## <-- shows the plot
    log.info("Stopping analysis here... modify code to continue! ")    
    sys.exit() ## quit for now...
    ## ... delete until here :)

    
    ##                _ __________
    ##  ___ ___      / |___ /___  |
    ## / __/ __|_____| | |_ \  / /
    ##| (__\__ \_____| |___) |/ /
    ## \___|___/     |_|____//_/
    ##    

    ## read in the Cs137 file measured without Al plate
    cs137 = read_mca_data_file('data/cs137.Spe')

    if not cs137:
        ## no data file could be loaded..
        sys.exit()
    
    ## plot into a new figure
    plt.figure()
    plt.grid(True)
    plt.xlabel('channel number')
    plt.ylabel('counts')
    plt.title("Cs-137")

    ## plot the cs-137 data
    plt.plot(cs137.x, cs137.y, marker='o',label="Cs-137")

    ## Side quest: Is that peak in the data really electrons? And what about the gamma background?
    ## -- load a file with a measurement of only the gamma background (how to do that?)
    ## => read in the background measurement file here using 'read_mca_data_file'
    ## -- make sure that the background is correctly normalized to the CS-137 data by taking measurement times into account
    ## => use the "duration" property of our spectra to calculate the factor and use "scale" method to scale the _background_
    ## -- subtract the background from the CS-137 measurement data
    ## => use the 'subtract' method of the Spectrum class 
    ## -- plot both the background and the CS-137 data w/o background into the same plot
    ## => how much background is there? where is it located? will it affect the fit precision?

    ## fit all gaussians in our measurement
    fits = fit_gaussians_to_measurement(cs137)
    ## loop over fit results
    for g in fits:
        ## plot the gaussian fit
        plt.plot(cs137.x, g.value(cs137.x), label="Gauss fit, $\sigma$="+str(g.sigma))

    ## now we have some data for our energy calibration:
    ## peak 0: internal conversion peak, 0.630 MeV

    if len(fits)>=1: ## check if we have at least one peak
        ## now create two arrays, one holding the channel number of the fitted peak
        ## the other the corresponding energy
        ecalib_cs137_channel = np.array( [ fits[0].mu ] )   ## [Ch# peak1, Ch# peak2, .... ]
        ecalib_cs137_energy  = np.array( [ 0.630 ] )        ## [E peak1, E peak2, ... ] in MeV
        ## copy the sigma (the width of the gaussian) from the fit results
        ecalib_cs137_sigma   = np.array( [ fits[0].sigma ] )## [Sigma] values for each peak
    else:
        ## could not reliably determine peak... need manual intervention!
        log.error("Unexpected number of Gaussians found in spectrum! Please modify fit parameters!")
        plt.show()     ## <-- shows the plot
        sys.exit()     ## quit 
                
    ## generate the legend (with the "label" information from the plots)
    plt.legend()
    ## to set the x-axis range ([xmin, xmax]), use ylim() to set y axis limits
    plt.xlim([0, 512])

    ## Delete this to continue!
    plt.show()           ## <-- shows the plot
    log.info("Stopping analysis here... modify code to continue! ")    
    sys.exit() ## quit for now...
    ## ... delete until here :)

    ## ____  _      ____   ___ _____
    ##| __ )(_)    |___ \ / _ \___  |
    ##|  _ \| |_____ __) | | | | / /
    ##| |_) | |_____/ __/| |_| |/ /
    ##|____/|_|    |_____|\___//_/
    ##
    ## ... nothing here yet ... time to copy&paste from Cs-137!
    
        
    ##                                               _ _ _               _   _                 
    ##  ___ _ __   ___ _ __ __ _ _   _      ___ __ _| (_) |__  _ __ __ _| |_(_) ___  _ __      
    ## / _ \ '_ \ / _ \ '__/ _` | | | |    / __/ _` | | | '_ \| '__/ _` | __| |/ _ \| '_ \     
    ##|  __/ | | |  __/ | | (_| | |_| |   | (_| (_| | | | |_) | | | (_| | |_| | (_) | | | |    
    ## \___|_| |_|\___|_|  \__, |\__, |    \___\__,_|_|_|_.__/|_|  \__,_|\__|_|\___/|_| |_|    
    ##                     |___/ |___/
    ## plot into a new figure
    plt.figure()
    plt.grid(True)
    plt.xlabel('channel number')
    plt.ylabel('energy [MeV]')
    plt.title("Energy Calibration")

    ## PLOT the data from Cs-137 including uncertainties (based on width of Gaussian fit)
    plt.errorbar(x=ecalib_cs137_channel, y=ecalib_cs137_energy, xerr=ecalib_cs137_sigma, marker='o',label="Cs-137")

    ## might want to COMBINE data arrays from different calibration sources for the fit:
    ## use
    ## new_array = np.concatenate( array1, array2 )
    ## to do so. Then change the data set in the fit command.
        
    ## linear regression of the data
    ## http://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.linregress.html
    slope = 1.
    intercept = 0.
    # .... something is missing here....
    
    ## ALTERNATIVE METHODS TO FIT:
    ## use "curve_fit" which allows to take uncertainties into account!
    ## http://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.curve_fit.html

    log.info("Determined calibration constants from linear regression: E [MeV] = "+str(slope)+"*N_ch + " + str(intercept))
    x = np.arange(1,512)
    plt.plot(x,slope*x+intercept,label="linear regression")
    plt.legend()
    
    ## apply energy calibration
    log.info("Applying calibration constants")
    sr90.calibrate(slope,intercept)
    
    ## plot into a new figure
    plt.figure()
    plt.grid(True)
    plt.xlabel('energy [MeV]')
    plt.ylabel('counts')
    plt.title("Sr-90/Y-90 energy spectrum")
    plt.plot(sr90.x, sr90.y, marker='o')

    ## Delete this to continue!
    plt.show()           ## <-- shows the plot
    log.info("Stopping analysis here... modify code to continue! ")    
    sys.exit() ## quit for now...
    ## ... delete until here :)

    
    ## _____                   _       _  __          _        ____  _       _
    ##|  ___|__ _ __ _ __ ___ (_)     | |/ /   _ _ __(_) ___  |  _ \| | ___ | |_
    ##| |_ / _ \ '__| '_ ` _ \| |_____| ' / | | | '__| |/ _ \ | |_) | |/ _ \| __|
    ##|  _|  __/ |  | | | | | | |_____| . \ |_| | |  | |  __/ |  __/| | (_) | |_
    ##|_|  \___|_|  |_| |_| |_|_|     |_|\_\__,_|_|  |_|\___| |_|   |_|\___/ \__|    
    ## fermi-kurie calculations:
    mec2 = 0.510998910 ## MeV
    pc = np.sqrt((sr90.x + mec2)**2 - mec2**2)
    A = pc/mec2
    f = 1.3604*A*A + 0.1973*A + 0.0439
    Ee = (sr90.x + mec2)
    QminTe = np.sqrt((sr90.y*pc)/(Ee*f))

    ## plot into a new figure
    plt.figure()
    plt.grid(True)
    plt.xlabel('Te [MeV]')
    plt.ylabel('Q-Te')
    plt.title("Sr-90 Fermi-Kurie")
    plt.plot(sr90.x, QminTe, marker='o', label="data")
    
    ## linear regression of the FM plot see
    ## http://docs.scipy.org/doc/scipy-0.16.1/reference/generated/scipy.stats.linregress.html
    ## the fit does not really work on the edges of the FM plot, so we
    ## take the linear region of lower_limit<E [MeV]<upper_limit (to
    ## be determined from the plot)
    lower_limit, upper_limit = 1,2 ## initialize
    try:
        # search for the bins that match our criteria
        lower_limit = np.where(sr90.x>0.2)[0][0] ## first elements indicate first bin matching our criteria
        upper_limit = np.where(sr90.x>1.5)[0][0]
    except IndexError:
        log.error("Could not find any bins to fit! Please check the limit settings!")

    slope, intercept, r_value, p_value, std_err = stats.linregress(sr90.x[lower_limit:upper_limit], QminTe[lower_limit:upper_limit])
    x = np.arange(0,2.5,0.05) ## generate x axis for fit result (start, stop, stepsize)
    plt.plot(x,slope*x+intercept,label="linear regression")

    plt.legend()

    ## now the Q value is determined by where the linear regression intersects with the x axis (Q-Te = 0)
    Q = -intercept/slope

    ## print results
    log.info("Determined linear regression to Fermi-Kurie plot: Q-Te = "+str(slope)+"*Te + " + str(intercept))
    log.info("===> Q value: Q = "+str(Q)+" MeV ")
    ## label plot with Q value
    plt.annotate('Q = '+"{:.3f}".format(Q)+' MeV', # text to put there
                xy=(Q, 0),                         # where to point to
                xytext=(0., 60),                   # offset for text
                textcoords='offset points',
                arrowprops=dict(arrowstyle="->", color = 'red'))
    
    ## final step:
    plt.show()           ## <-- shows the plot
