## \package MD.analysis.histogram 
# \brief This module is used for creating a histogram of data

import numpy as np
from scipy import interpolate

## \brief return normalized histogram which has been smootthed 
# These values are normalizd to unity
# \returns hist_s  y values of hist
# \returns xs  x values for hist 
# \returns max_hist  max value in histogram
#
# \param distance: distance values between points
# \param bins: number of bins to use
# \param end: cutoff for the end of the histogram
def histogram(distance,bins,max_hist=0,end=1):
    hist,bins=np.histogram(distance,bins)
    bin2=[]
    hist2=[]
    for i in range(len(bins)-end):
        bin2.append(bins[i+1])
        hist2.append(hist[i])
    bins=np.array(bin2)
    hist=np.array(hist2)
    #normalize the histogram so that the highest point is 1
    if float(max(hist))>max_hist:
        max_hist=float(max(hist))
    hists=[]
    for i in range(hist.shape[0]):
        hists.append(hist[i]/max_hist)
    hist=np.array(hists)
    s = interpolate.InterpolatedUnivariateSpline(bins, hist)
    xs = np.arange(bins[0],bins[-1],0.1)
    hist_s=s(xs)
    return hist_s,xs,max_hist

## \brief return histogram which has been normalized 
#
# \returns hist_s  y values of hist
# \returns xs  x values for hist 
# \returns max_hist max value in histogram
#
# \param distance: distance values between points
# \param bins: number of bins to use
# \param end: cutoff for the end of the histogram
def histogram_normal(distance,bins,max_hist=0,end=1):
    if type(bins) is int:
        hist,bins=np.histogram(distance,bins)
    else:
        hist,bins=np.histogram(distance,bins=bins)
    bin2=[]
    hist2=[]
    for i in range(len(bins)-end):
        bin2.append((bins[i+1]+bins[i])/2)
        hist2.append(hist[i])
    bins=np.array(bin2)
    hist=np.array(hist2)
    if float(max(hist))>max_hist:
        max_hist=float(max(hist))
        print max_hist
    hists=[]
    max_hist = float(sum(hist))
    for i in range(hist.shape[0]):
        hists.append(hist[i]/max_hist)
    hist=np.array(hists)
    return hist,bins,max(hist)


## \brief return histogram which is ready for plotting 
#
# \returns hist_s  y values of hist
# \returns xs  x values for hist 
# \returns max_hist max value in histogram
#
# \param distance: distance values between points
# \param bins: number of bins to use
# \param end: cutoff for the end of the histogram
def histogram_reg(distance,bins,max_hist=0,end=1):
    hist,bins=np.histogram(distance,bins)
    bin2=[]
    hist2=[]
    for i in range(len(hist)-end):
        bin2.append(bins[i+1])
        hist2.append(hist[i])
    bins=np.array(bin2)
    hist=np.array(hist2)
    return hist,bins
