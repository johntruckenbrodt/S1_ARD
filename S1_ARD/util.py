import math
import time

import numpy as np
from numpy.polynomial.polynomial import polyfit

from osgeo import gdal

import matplotlib.pyplot as plt
from matplotlib.offsetbox import AnchoredText

from scipy import stats
from scipy.stats import gaussian_kde
from sklearn.metrics import mean_squared_error, r2_score

from spatialist import Raster

# Function to generate one to one plots for each land cover class
# from specified images. (Utilises function one.)
def one2oneSARplots(in_sar1_img, in_sar2_img, band, sar1_ref, sar2_ref, lc, lc_ints, lc_string):
    # Read data.
    sar1 = readimg2array(in_sar1_img, band, 1)
    sar2 = readimg2array(in_sar2_img, band, 1)
    # Iterate land covers.
    for p, i in enumerate(lc_ints):
        start_time = time.time()
        # Calculate correlation statistics (RMSE & R^2).
        s1_c = sar1[lc == i]
        s2_c = sar2[lc == i]
        
        sar3 = s1_c[(np.isfinite(s1_c)) & (np.isfinite(s2_c))]
        sar4 = s2_c[(np.isfinite(s1_c)) & (np.isfinite(s2_c))]
        print(np.shape(sar3))
        print(np.shape(sar4))
        lowest = min([min(sar3), min(sar4)])
        highest = max([max(sar3), max(sar4)])
        o2o_min = [lowest, highest]
        o2o_max = [lowest, highest]
        
        num = len(sar3)
        
        rmse = math.sqrt(mean_squared_error(sar3, sar4))
        print(rmse)
        r2 = r2_score(sar3, sar4)
        print(r2)
        
        slope, intercept, r_value, p_value, std_err = stats.linregress(sar3, sar4)
        line = slope * s1_c + intercept
        print(r_value ** 2, p_value, std_err, intercept, slope)
        
        # Plot SAR vars sliced to pixels of land cover.
        plt.scatter(sar1[lc == i], sar2[lc == i], s=0.1)
        plt.xlabel('%s' % (sar1_ref))
        plt.ylabel('%s' % (sar2_ref))
        plt.title('%s \n %s vs %s \n N=%s' % (lc_string[p], sar1_ref, sar2_ref, num))
        plt.xlim(lowest, highest)
        plt.ylim(lowest, highest)
        # Add one to one line.
        plt.plot(o2o_min, o2o_max, color="black")
        
        # Add best fit line.
        #   plt.plot(np.unique(s1_c), np.poly1d(np.polyfit(s1_c, s2_c, 1))(np.unique(s1_c)))
        plt.plot(s1_c, line, color="red")
        
        plt.grid()
        plt.show()
        print("--- %s seconds ---" % (time.time() - start_time))
    sar1 = None
    sar2 = None
    sar3 = None
    sar4 = None


def one2oneSARplotsratio(in_sar1_img, in_sar2_img, vv_band, vh_band, sar1_ref, sar2_ref, lc, lc_ints, lc_string):
    # Read data.
    sar1 = readimg2array(in_sar1_img, vv_band, 2) / readimg2array(in_sar1_img, vh_band, 2)
    sar2 = readimg2array(in_sar2_img, vv_band, 2) / readimg2array(in_sar2_img, vh_band, 2)
    sar1 = 10 * np.log10(sar1)
    sar2 = 10 * np.log10(sar2)
    
    for p, i in enumerate(lc_ints):
        start_time = time.time()
        
        s1_c = sar1[lc == i]
        s2_c = sar2[lc == i]
        sar3 = s1_c[(np.isfinite(s1_c)) & (np.isfinite(s2_c))]
        sar4 = s2_c[(np.isfinite(s1_c)) & (np.isfinite(s2_c))]
        print(np.shape(sar3))
        print(np.shape(sar4))
        lowest = min([min(sar3), min(sar4)])
        highest = max([max(sar3), max(sar4)])
        o2o_min = [lowest, highest]
        o2o_max = [lowest, highest]
        
        num = len(sar3)
        
        rmse = math.sqrt(mean_squared_error(sar3, sar4))
        print(rmse)
        r2 = r2_score(sar3, sar4)
        print(r2)
        
        slope, intercept, r_value, p_value, std_err = stats.linregress(sar3, sar4)
        line = slope * s1_c + intercept
        print(r_value ** 2, p_value, std_err, intercept, slope)
        
        # Plot SAR vars sliced to pixels of land cover.
        plt.scatter(sar1[lc == i], sar2[lc == i], s=0.1)
        plt.xlabel('%s' % (sar1_ref))
        plt.ylabel('%s' % (sar2_ref))
        plt.title('%s \n %s vs %s \n N=%s' % (lc_string[p], sar1_ref, sar2_ref, num))
        plt.xlim(lowest, highest)
        plt.ylim(lowest, highest)
        # Add one to one line.
        plt.plot(o2o_min, o2o_max, color="black")
        
        # Add best fit line.
        #   plt.plot(np.unique(s1_c), np.poly1d(np.polyfit(s1_c, s2_c, 1))(np.unique(s1_c)))
        plt.plot(s1_c, line, color="red")
        
        plt.grid()
        plt.show()
        print("--- %s seconds ---" % (time.time() - start_time))
    sar1 = None
    sar2 = None
    sar3 = None
    sar4 = None


# Function to read first band of image and convert to dB
# is SAR input is specified.
def readimg2array(in_img, band, sar):
    start_time = time.time()
    ds = gdal.Open(in_img)
    band = ds.GetRasterBand(band)
    cols = ds.RasterXSize
    rows = ds.RasterYSize
    print(cols, rows)
    print("--- %s seconds ---" % (time.time() - start_time))
    if sar == 1:
        x = band.ReadAsArray()
        x[(x == -32768) | (x == 0.0)] = 'NaN'
        return 10 * np.log10(x)
    elif sar == 2:
        x = band.ReadAsArray()
        x[(x == -32768) | (x == 0.0)] = 'NaN'
        return x
    else:
        x = band.ReadAsArray()
        return x


def simplify_lc(in_lc):
    in_lc = np.where((in_lc <= 11), 1, in_lc)
    in_lc = np.where((in_lc >= 12) & (in_lc <= 22), 12, in_lc)
    in_lc = np.where((in_lc >= 23) & (in_lc <= 34), 23, in_lc)
    in_lc = np.where((in_lc >= 35) & (in_lc <= 39), 35, in_lc)
    in_lc = np.where((in_lc >= 40) & (in_lc <= 44), 40, in_lc)
    return in_lc


def sar_vs_inc(file_sar, file_inc, nsamples, nodata=-99, db_convert=False, title='', xlabel='', ylabel='',
               regfun=False, ymin=None, ymax=None):
    with Raster(file_sar) as ras:
        # print(ras)
        sar = ras.matrix()
    
    with Raster(file_inc) as inc:
        # print(inc)
        inc = inc.matrix()
    
    inc = inc * 180 / math.pi
    
    sar[sar == nodata] = np.nan
    
    mask = ~np.isnan(sar)
    
    step = int(math.floor(sar.size / nsamples))
    sar_sub = sar[mask][0::step]
    inc_sub = inc[mask][0::step]
    
    if db_convert:
        sar_sub = 10 * np.log10(sar_sub)
    
    # # Calculate the point density
    xy = np.vstack([sar_sub, inc_sub])
    z = gaussian_kde(xy)(xy)
    
    # # Sort the points by density, so that the densest points are plotted last
    idx = z.argsort()
    x, y, z = inc_sub[idx], sar_sub[idx], z[idx]
    
    # fig, ax = plt.subplots()
    plt.scatter(x, y, c=z, s=20, edgecolor='')
    
    # add linear regression equation
    if regfun:
        b, m = polyfit(x, y, 1)
        text_box = AnchoredText('y = {:.2f} + {:.2f} * x'.format(b, m), frameon=True, loc='lower left', pad=0.1)
        plt.setp(text_box.patch, facecolor='white', alpha=0.5)
        plt.gca().add_artist(text_box)
    
    plt.title(title)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)

    # set axes range
    bottom, top = plt.ylim()
    plt.ylim(ymin if ymin is not None else bottom, ymax if ymax is not None else ymax)
