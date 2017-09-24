#!/usr/bin/env python3
import csv                               ## for reading in our data files
import logging                           ## for orderly print output
import numpy as np                       ## for handling of data
import os                                ## path manipulations

class Spectrum:
    """A class to hold our spectrum measurement data and meta data (such as duration)"""
    def __init__(self, filename):
        self.filename = filename
        self.x = np.array(np.zeros(1))    ## creates a new empty array; will later store our x values
        self.y = np.array(np.zeros(1))    ## creates a new empty array; will later store our y values
        self.name = os.path.splitext(     ## a more descriptive name, can be used e.g. in legends
                os.path.basename(filename))[0] ## init with file name without extension
        self.duration = 0
    def subtract(self, m):
        self.y = self.y - m.y
        ## these are spectra: cannot have counts below 0, so remove these here and set them to 0 instead:
        self.y[self.y < 0] = 0;
    def scale(self, scale):
        self.y *= scale
    def calibrate(self, slope, intercept):
        self.x = self.x*slope + intercept

def load_spectrum(filename):
    """Reads in a data file (csv format) stored by the Maestro MCA software and returns a 'Spectrum' object. Tested with Maestro Version 6.05 """
    log = logging.getLogger('load_spectrum') ## set up logging
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
            nchannel = int(interval[1]-interval[0])+1
            m.y = np.array(np.zeros(nchannel))
            ## continue, now reading data
            for idx, row in enumerate(reader):
                if idx >= nchannel:
                    break
                m.y[idx] = int(row[0])
            m.x = np.arange(interval[0], interval[1]+1,1)
            log.debug("Loaded all data from file")
    except IOError:
        log.error("Could not find the file '"+str(filename)+"'")
        return None
    return m

#This function reads the calibrated background spectrum that is to be analysed in the gamma lab.
def load_calibrated_spectrum(filename):
    log = logging.getLogger('gammalab_analysis') ## set up logging
    m = Spectrum(filename) ## create a new Spectrum measurement object; this is what we return in the end
    m.x = np.zeros(8192)
    m.y = np.zeros(8192)
    log.info("Reading calibrated data from file '" + filename + "'")
    try:
        with open(filename) as f: #the with keyword handles the opening (__enter__ method) and closing (__exit__ method) of the file automatically
            reader = csv.reader(f)
            for idx, row in enumerate(reader):
                channel, energy = row[0].split() 
                m.y[idx] = int(channel)
                m.x[idx] = float(energy)
    except IOError:
        log.error("Could not find the file '"+str(filename)+"'")
        sys.exit(-1)
    return m
