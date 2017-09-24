#!/usr/bin/env python3
import numpy as np                       ## for handling of data
from scipy.optimize import curve_fit     ## to fit functions to the data
from scipy import signal                 ## to find maxima in the data
from scipy import stats                  ## for linear regressions
import sys                               ## useful system calls (used to exit cleanly)


def gaussfcn(x, *p):
    """ gauss function to be used for fits to the data"""
    A, mu, sigma = p
    return A*np.exp(-(x-mu)**2/(2.*sigma**2))

class Gauss:
    """A class to hold coefficients for Gaussian distributions"""
    def __init__(self, A, mu, sigma, covar_matrix):
        self.A = A
        self.mu = mu
        self.sigma = sigma
        self.covar_matrix = covar_matrix
    def value(self, x):
        return gaussfcn(x, self.A, self.mu, self.sigma)
    def as_string(self, ndigits=4):
        return str("A: {}, mu: {}, sigma: {}".format(round(self.A, ndigits),
                                                     round(self.mu, ndigits),
                                                     round(self.sigma, ndigits)))

def line(x, *p):
    """ straight line function to be used for fits to the data"""
    a, b = p
    return a*x+b

def fit_gaussian_at_idx(x, y, idx, npoints=10):
    """ takes a spectrum measurement and an index for a position in the data (x/y) where
    a Gaussian fit should be performed. Takes into account N surrounding data points given by 'npoints' argument """
    ## p0 is the initial guess for the fitting coefficients (A, mu and sigma above)
    p0 = [y[idx], x[idx], 1.] ## mu is given by one of the found peaks positions
    ## scypi gives a warning if the fit does not work; we want to know about those, so we set them up to be caught here:
    import warnings
    from scipy.optimize import OptimizeWarning
    with warnings.catch_warnings():
        warnings.simplefilter("error", OptimizeWarning)
        ## perform the gaussian fit to the data:
        try:
            ## use the scipy curve_fit routine (uses non-linear least squares to perform the fit)
            ## see http://docs.scipy.org/doc/scipy-0.16.1/reference/generated/scipy.optimize.curve_fit.html
            ## fit using "slices" of the arrays with +/- 10 around peak (or less if out of bounds)
            uppvar = min(len(x)-idx, npoints)
            lowvar = min(idx, npoints)
            coeff, var_matrix = curve_fit(gaussfcn,
                                          x[idx-lowvar:idx+uppvar],
                                          y[idx-lowvar:idx+uppvar],
                                          p0=p0)
            ## create a Gauss object with the fitted coefficients for better code readability
            g = Gauss(*coeff, var_matrix)
            return g
        except (RuntimeError, OptimizeWarning, TypeError):
            ## the minimization did not work out... log it and continue to next peak
            print("  - gaussian fit failed!")

def fit_gaussian_at_pos(x, y, pos, npoints=10):
    """ fits x,y values at given x position with a Gaussian. """
    g = fit_gaussian_at_idx(x, y, idx=np.where(x>=pos)[0][0], npoints=npoints)
    if g is None:
        print("Fit at x = {}: failed! :(".format(pos))
    return g

def fit_all_gaussians(x, y, npoints=10, widths = np.arange(10,80,5), loglevel="WARNING"):
    """ fits all gaussians in a spectrum measurement and returns a list of coefficients. 
    The range of widths considered for fit is given by an array (e.g. 'np.arange(X1,X2,X3)': 
    range from X1 to X2 in steps of X3)."""

    ## list to store the parameters of the fitted gaussians in
    gaussians = []

    ## find peaks in y with range of widths given by an array (range from X1 to X2 in steps of X3)
    peakind = signal.find_peaks_cwt(y, widths) 
    print("Found {} peak(s) in the data, fitting them with gaussians:".format(len(peakind)))
    for p in peakind:
        g = fit_gaussian_at_idx(x, y, p, npoints=npoints)
        if not g:
            continue
        ## filter the results -- or we get a lot of "spurious" peaks
        xdynrange = x.shape[0] ## dynamic range in x
        if g.sigma > 0.1*xdynrange: ## check width of gaussian in percent of the dynamic range
            print("  - sigma out of bounds: " + str(g.sigma))
            continue
        if g.mu > xdynrange-0.1*xdynrange or g.mu < 6: ## check center position: should not be at the limits of measurement range
            print("  - mu out of bounds: " + str(g.mu))
            continue
        if g.A < 0.25*np.average(y): ## check if the peak hight is at least 25% over the average data value
            print("  - A out of bounds: " + str(g.A))
            continue
        print("  - fit result: A = " + str(g.A) + ", mu = " + str(g.mu) + ", sigma = " + str(g.sigma) + ". ")
        ## store the results
        gaussians.append(g)
    return gaussians
