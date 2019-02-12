import math
import time
import numpy as np
from numpy.polynomial.polynomial import polyfit

from scipy import stats
from scipy.stats import gaussian_kde
from sklearn.metrics import mean_squared_error, r2_score

from astropy.convolution import convolve, CustomKernel

import matplotlib.pyplot as plt
from matplotlib.offsetbox import AnchoredText

from osgeo import gdal


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


def sar_vs_inc(sar, inc, nsamples=1000, nodata=-99, db_convert=False, title='', xlabel='', ylabel='',
               regfun=False, ymin=None, ymax=None, mask=None):
    inc = np.rad2deg(inc)
    
    sar[sar == nodata] = np.nan
    
    if mask is not None:
        sar[~mask] = np.nan
    
    nanmask = (~np.isnan(sar)) & (~np.isnan(inc))

    sample_ids = sampler(nanmask, nsamples)
    
    sar_sub = sar.flatten()[sample_ids]
    inc_sub = inc.flatten()[sample_ids]
    
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


def dem_aspect(img):
    boundary = 'extend'
    kernel = CustomKernel(np.array([[-1, 0, 1], [-2, 0, 2], [-1, 0, 1]]) / 8.)
    xchangerate_array = -convolve(img, kernel, normalize_kernel=False, boundary=boundary,
                                  nan_treatment='fill', fill_value=np.nan)
    
    kernel = CustomKernel(np.array([[-1, -2, -1], [0, 0, 0], [1, 2, 1]]) / 8.)
    ychangerate_array = -convolve(img, kernel, normalize_kernel=False, boundary=boundary,
                                  nan_treatment='fill', fill_value=np.nan)
    
    aspect = np.rad2deg(np.arctan2(ychangerate_array, -xchangerate_array))
    aspect_value = np.copy(aspect)
    
    # numpy raises a warning if np.nan values are set to False in a mask
    with np.errstate(invalid='ignore'):
        mask = aspect < 0.0
        aspect_value[mask] = 90.0 - aspect[mask]
        
        mask = aspect > 90.0
        aspect_value[mask] = 360.0 - aspect[mask] + 90.0
        
        mask = (aspect >= 0.0) & (aspect < 90.0)
        aspect_value[mask] = 90.0 - aspect[mask]
    return aspect_value


def dem_distribution(slope, aspect, head_angle, inc_angle, look_dir='right', nsamples=1000):
    
    nanmask = (~np.isnan(slope)) & (~np.isnan(aspect))
    sample_ids = sampler(nanmask, nsamples)
    
    slope = slope.flatten()[sample_ids]
    aspect = aspect.flatten()[sample_ids]
    
    # # Calculate the point density
    xy = np.vstack([slope, aspect])
    z = gaussian_kde(xy)(xy)
    
    # # Sort the points by density, so that the densest points are plotted last
    idx = z.argsort()
    x, y, z = aspect[idx], slope[idx], z[idx]
    
    # create visible angle map mesh grid
    xx, yy = np.meshgrid(np.radians(np.arange(0, 361)),
                         np.arange(0, 91),
                         sparse=False)
    vis_map = visible_sar_angle_map(head_angle, inc_angle, look_dir)
    
    # get current axis if a figure exists
    ax = plt.gca() if plt.get_fignums() else None
    
    # create a new figure if there is no axis or the projection of the current axis is not polar
    if ax is None or not ax.__class__.__name__ == 'PolarAxesSubplot':
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='polar')
    ax.set_theta_zero_location('N')
    c = ax.scatter(x, y, c=z, s=10, alpha=0.75)
    ax.contour(xx, yy, vis_map, colors='red')
    
    ax.set_ylim(0, 95)


def dem_slope(img, x_cell_size, y_cell_size):
    boundary = 'extend'
    kernel = CustomKernel(np.array([[-1, 0, 1], [-2, 0, 2], [-1, 0, 1]]) / (8. * x_cell_size))
    xchangerate_array = convolve(img, kernel, normalize_kernel=False, boundary=boundary,
                                 nan_treatment='fill', fill_value=np.nan)
    
    kernel = CustomKernel(np.array([[-1, -2, -1], [0, 0, 0], [1, 2, 1]]) / (8. * y_cell_size))
    ychangerate_array = convolve(img, kernel, normalize_kernel=False, boundary=boundary,
                                 nan_treatment='fill', fill_value=np.nan)
    
    slope_radians = np.arctan(np.sqrt(np.square(xchangerate_array) + np.square(ychangerate_array)))
    slope_degrees = np.rad2deg(slope_radians)
    return slope_degrees


def sampler(nanmask, nsamples=None, seed=42):
    indices = np.where(nanmask.flatten())[0]
    samplesize = min(indices.size, nsamples) if nsamples is not None else indices.size
    np.random.seed(seed)
    sample_ids = np.random.choice(a=indices, size=samplesize, replace=False)
    return sample_ids


def visible_sar_angle_map(head_angle, inc_angle, look_dir='right'):
    # convert deg to rad and geographic heading to mathematical angle
    head_ang = math.radians(90. - head_angle)
    inc_ang = math.radians(inc_angle)
    
    # flight direction vector
    vel_vec = np.array([math.cos(head_ang), math.sin(head_ang), 0])
    
    # viewing plane (look vector)
    if look_dir == 'right':
        look_vec = np.array([math.sin(inc_ang) * math.sin(head_ang),
                             -math.sin(inc_ang) * math.cos(head_ang),
                             -math.cos(inc_ang)])
        viewplane = np.cross(look_vec, vel_vec)
    
    elif look_dir == 'left':
        look_vec = np.array([-math.sin(inc_ang) * math.sin(head_ang),
                             math.sin(inc_ang) * math.cos(head_ang),
                             -math.cos(inc_ang)])
        viewplane = np.cross(vel_vec, look_vec)
    else:
        raise ValueError("look_dir must either be 'left' or 'right'")
    
    # map in aspect and slope coordinates
    angle, radius = np.meshgrid(np.radians(np.arange(0, 361)),
                                np.radians(np.arange(0, 91)),
                                sparse=False)
    
    # surface normal vector
    normal = np.array([np.sin(radius) * np.cos(angle),
                       np.sin(radius) * np.sin(angle),
                       np.cos(radius)])
    
    # condition that facet is visible by radar
    view_prod = normal[0, :, :] * look_vec[0] + \
                normal[1, :, :] * look_vec[1] + \
                normal[2, :, :] * look_vec[2]
    
    # condition that facet is tilted positive in viewplane
    mont_prod = normal[0, :, :] * viewplane[0] + \
                normal[1, :, :] * viewplane[1] + \
                normal[2, :, :] * viewplane[2]
    
    # visible map_array
    visible_map_array = (view_prod < 0) & (mont_prod > 0)
    
    return visible_map_array
