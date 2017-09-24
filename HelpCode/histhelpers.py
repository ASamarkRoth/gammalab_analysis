""" Helper routines to import histograms stored in csv files as exported by pyroexport from a ROOT file
example use:
    h = import_hist_from_csv('/Users/hperrey/src/pyroplot/hc1_4_0.csv')
    import matplotlib.pyplot as plt
    plt.plot(h.getBinCenters(), h.data)
    plt.show()

"""


import csv
import numpy as np
import logging
import sys
import os

class AquaHist:
    """A class to hold AquaDAQ measurement data and meta data"""
    def __init__(self, filename):
        self.filename = filename
        self.data = np.array(np.zeros(1))    ## creates a new empty histogram; will later store our data
        self.bins = np.array(np.zeros(1))
        self.name = os.path.splitext(     ## a more descriptive name, can be used e.g. in legends
                os.path.basename(filename))[0] ## init with file name without extension
        self.title = ""
        self.xAxisTitle = ""
        self.yAxisTitle = ""

def getBinCenters(bins):
    """ calculate center values for given bins """
    return np.array([np.mean([bins[i],bins[i+1]]) for i in range(0, len(bins)-1)])
def getRandDist(bins, data):
    """ calculates a random distribution for each bin based on bin content, thus allowing rebinning of histogram """
    return np.array([
        # get random number in the range lowedge < x < highedge for the current bin
        (bins[binno + 1] - bins[binno])*np.random.random_sample() + bins[binno]
        # loop over all bins
        for binno in range(0, len(bins)-1)
        # for each entry in bin; convert float to int; does not handle negative entries!
        for count in range (1, max(0,int(round(data[binno]))))
    ])

def getEvenDist(bins, data):
    """ calculates a random distribution for each bin based on bin content, thus allowing rebinning of histogram """
    return np.array([
        # get fraction number in the range lowedge < x < highedge for the current bin
        (bins[binno + 1] - bins[binno])*(count/data[binno]) + bins[binno]
        # loop over all bins
        for binno in range(0, len(bins)-1)
        # for each entry in bin; convert float to int; does not handle negative entries!
        for count in range (1, max(0,int(round(data[binno]))))
    ])



def import_hist_from_csv(file):
    logging.basicConfig(format='%(asctime)s %(message)s')
    log = logging.getLogger(__name__)
    # now open the input file:
    log.debug("Opening input file {}".format(file))
    with open(file, 'rt') as csvfile:
        reader = csv.reader(csvfile)
        ## first line should contain "pyroexport" string or this is likely an unsupported file type
        row = next(reader)
        if not 'PYROEXPORT' in row[0]:
            log.error("File '"+str(file)+"' is not a valid pyroexport histogram file!")
            sys.exit(-1)

        # set up the AquaHist object to be returned later
        h = AquaHist(file)
        row = next(reader) # parse next line
        h.name = os.path.basename(row[1]) # the file name of the original ROOT file
        row = next(reader) # parse next line
        h.title = row[1]
        row = next(reader)
        h.xAxisTitle = row[1]
        row = next(reader)
        h.yAxisTitle = row[1]

        # read number of bins and initialize the arrays to store the data
        row = next(reader)
        nbins = int(row[1])
        h.data = np.array(np.zeros(nbins))
        h.bins = np.array(np.zeros(nbins+1))
        row = next(reader) # skip one line
        row = next(reader) # skip column headers

        ## continue, now reading remaining data
        for idx, row in enumerate(reader):
            if (idx) >= nbins:
                break
            h.data[idx] = float(row[1])
            h.bins[idx] = float(row[0])
        h.bins[nbins] = float(row[0]) # complete array with final bin edge

        log.debug("Loaded all data from file")
        return h

if __name__ == '__main__':
    print("Loading data...")
    h = import_hist_from_csv('hc1_4_0.csv')
    import matplotlib.pyplot as plt
    plt.plot(getBinCenters(h.bins), h.data, label="data")
    from timeit import default_timer as timer
    print("Testing 'rebin' module on data...")
    start = timer()
    import rebin
    def shift(tof1, tof2, amount):
        # rebin the (shifted) data to match the first spectrum
        tof2.data = rebin.rebin(tof2.bins + amount, tof2.data, tof1.bins)
        # adjust the binning to match spectrum one's
        tof2.bins = tof1.bins
    h2 = h
    shift(h, h2, 20)
    plt.plot(getBinCenters(h2.bins), h2.data, label="data")
#    newdata = rebin.rebin(h.bins, h.data, h.bins)
#    plt.plot(getBinCenters(h.bins), newdata, label="rebin")
    end = timer()
    print(".. which took {} seconds to run. ".format(end - start))
#    print("Testing internal random redistribution method on data...")
#    start = timer()
#    plt.hist(getRandDist(h.bins, h.data), bins=h.bins, label="random")
#    end = timer()
#    print(".. which took {} seconds to run. ".format(end - start))
#    print("Testing internal uniform redistribution method on data...")
#    start = timer()
#    plt.hist(getEvenDist(h.bins, h.data), bins=h.bins, label="even")
#    end = timer()
#    print(".. which took {} seconds to run. ".format(end - start))
    plt.legend()
    plt.show()
