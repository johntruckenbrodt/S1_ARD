import os
import math
import time
import shutil
import numpy as np
import pathos.multiprocessing as mp
from numpy.polynomial.polynomial import polyfit

from scipy import stats
from scipy.stats import gaussian_kde, variation
from sklearn.metrics import mean_squared_error, r2_score

from astropy.convolution import convolve, CustomKernel

import matplotlib.pyplot as plt
from matplotlib.offsetbox import AnchoredText

from osgeo import gdal, ogr
from osgeo.gdalconst import GA_Update

from spatialist import haversine, Raster, Vector, crsConvert, gdalwarp, gdal_translate


def scatter(x, y, xlab='', ylab='', title='', nsamples=1000, mask=None, measures=None,regline=False,
            o2o=False, denscol=False, grid=False, xlim=None, ylim=None):
    """
    general function for creating scatter plots
    
    Parameters
    ----------
    x: numpy.ndarray
        dataset I
    y: numpy.ndarray
        dataset II
    xlab: str
        the x-axis label
    ylab: str
        the y-axis label
    title: str
        the plot title
    nsamples: int
        the number of data samples to plot
    mask: numpy.ndarray
        an optional array for masking the datasets
    measures: list
        additional measures to be printed in a text box; current options:
         - `eq`: the linear regression equation
         - `rmse`
         - `r2`
         - `n`: the number of samples
         - `cv_x`, `cv_y`: the coefficient of variation of either x or y
    regline: bool
        draw a linear regression line?
    o2o: bool
        draw a data one to one line?
    denscol: bool
        color the points by Gaussian density?
    grid: bool
        add a mesh grid to the plot?
    xlim: tuple
        the x-axis limits
    ylim: tuple
        the y-axis limits

    Returns
    -------

    """
    nanmask = (np.isfinite(x)) & (np.isfinite(y))
    if mask is None:
        mask = nanmask
    else:
        mask = mask & nanmask
    
    sample_ids = sampler(mask, nsamples)
    
    x = x.flatten()[sample_ids]
    y = y.flatten()[sample_ids]
    
    measures = [] if measures is None else measures
    fields = []
    if regline or 'eq' in measures:
        b, m = polyfit(x, y, 1)
    if 'eq' in measures:
        fields.append('y = {:.2f} {} {:.2f} * x'.format(b, '+' if m > 0 else '-',  abs(m)))
    if 'rmse' in measures:
        rmse = round(math.sqrt(mean_squared_error(x, y)), 2)
        fields.append('RMSE = {:.2f}'.format(rmse))
    if 'r2' in measures:
        r2 = round(r2_score(x, y), 2)
        fields.append('$R^2$ = {:.2f}'.format(r2))
    if 'n' in measures:
        fields.append('n = {}'.format(len(x)))
    if 'cv_x' in measures:
        cv = round(float(variation(x)), 4)
        fields.append('CV(x) = {}'.format(cv))
    if 'cv_y' in measures:
        # cv = variation(y)
        cv = round(float(variation(y)), 4)
        fields.append('CV(y) = {}'.format(cv))
    text = '\n'.join(fields)
    
    if denscol:
        # # Calculate the point density
        xy = np.vstack([x, y])
        z = gaussian_kde(xy)(xy)
        
        # # Sort the points by density, so that the densest points are plotted last
        idx = z.argsort()
        x, y, z = x[idx], y[idx], z[idx]
    else:
        z = None
    plt.scatter(x, y, c=z, s=20, edgecolor='', zorder=3)
    
    # set axes range
    if xlim is not None:
        plt.xlim(*xlim)
    if ylim is not None:
        plt.ylim(*ylim)
    
    # determine x-axis limits and compute border offset for lines
    xmin, xmax = plt.xlim()
    xo = (xmax - xmin) / 100 * 2
    
    if len(text) > 0:
        text_box = AnchoredText(text, frameon=True, loc='lower right')
        plt.setp(text_box.patch, facecolor='white')  # , alpha=0.5
        plt.gca().add_artist(text_box)
    if regline:
        ffit = np.poly1d((m, b))
        x_new = np.linspace(xmin + xo, xmax - xo, num=2)
        plt.plot(x_new, ffit(x_new), color='red', zorder=2)
    if o2o:
        ffit = np.poly1d((1, 1))
        x_new = np.linspace(xmin + xo, xmax - xo, num=2)
        plt.plot(x_new, ffit(x_new), color='black', zorder=1)
    plt.title(title)
    plt.xlabel(xlab)
    plt.ylabel(ylab)
    if grid:
        plt.grid()


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


# Function to read first band of image and convert to dB is SAR input is specified.
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


def dem_distribution(slope, aspect, head_angle, inc_angle, look_dir='right', nsamples=1000, title=''):
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
    # ax.set_theta_zero_location('N')
    ax.set_theta_direction(-1)
    c = ax.scatter(x, y, c=z, s=10, alpha=0.75)
    ax.contour(xx, yy, vis_map, colors='red')
    ax.set_title(title)
    ax.set_ylim(0, 95)


def dem_slope(img, xres_m, yres_m):
    """
    compute the slope of a DEM
    
    Parameters
    ----------
    img: numpy.ndarray
        the input DEM
    xres_m: int or float
        the x resolution of the DEM in same units as the height values
    yres_m: int or float
        the x resolution of the DEM in same units as the height values

    Returns
    -------

    """
    boundary = 'extend'
    kernel = CustomKernel(np.array([[-1, 0, 1], [-2, 0, 2], [-1, 0, 1]]) / (8. * xres_m))
    xchangerate_array = convolve(img, kernel, normalize_kernel=False, boundary=boundary,
                                 nan_treatment='fill', fill_value=np.nan)
    
    kernel = CustomKernel(np.array([[-1, -2, -1], [0, 0, 0], [1, 2, 1]]) / (8. * yres_m))
    ychangerate_array = convolve(img, kernel, normalize_kernel=False, boundary=boundary,
                                 nan_treatment='fill', fill_value=np.nan)
    
    slope_radians = np.arctan(np.sqrt(np.square(xchangerate_array) + np.square(ychangerate_array)))
    slope_degrees = np.rad2deg(slope_radians)
    return slope_degrees


def dem_degree2meter(demfile):
    with Raster(demfile) as ras:
        res_lon, res_lat = ras.res
        lat = (ras.geo['ymin'] + ras.geo['ymax']) / 2
        lon = (ras.geo['xmin'] + ras.geo['xmax']) / 2
    post_north = haversine(lat, lon, lat + res_lat, lon)
    post_east = haversine(lat, lon, lat, lon + res_lon)
    return post_east, post_north


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


def wkt2shp(wkt, srs, outname):
    geom = ogr.CreateGeometryFromWkt(wkt)
    geom.FlattenTo2D()
    
    srs = crsConvert(srs, 'osr')
    
    layername = os.path.splitext(os.path.basename(outname))[0]
    
    with Vector(driver='Memory') as bbox:
        bbox.addlayer(layername, srs, geom.GetGeometryType())
        bbox.addfield('area', ogr.OFTReal)
        bbox.addfeature(geom, fields={'area': geom.Area()})
        bbox.write(outname, format='ESRI Shapefile')
    geom = None


def parallel_apply_along_axis(func1d, axis, arr, cores=4, *args, **kwargs):
    """
    Like numpy.apply_along_axis(), but takes advantage of multiple cores.
    https://stackoverflow.com/questions/45526700/easy-parallelization-of-numpy-apply-along-axis
    
    Parameters
    ----------
    func1d: function
        the function to be applied
    axis: int
        the axis along which to apply func1d
    arr: the input array
    cores: the number of parallel cores
    args: any
        Additional arguments to func1d.
    kwargs: any
        Additional named arguments to func1d.

    Returns
    -------
    numpy.ndarray
    """
    # Effective axis where apply_along_axis() will be applied by each
    # worker (any non-zero axis number would work, so as to allow the use
    # of `np.array_split()`, which is only done on axis 0):
    effective_axis = 1 if axis == 0 else axis
    if effective_axis != axis:
        arr = arr.swapaxes(axis, effective_axis)
    
    def unpack(arguments):
        func1d, axis, arr, args, kwargs = arguments
        return np.apply_along_axis(func1d, axis, arr, *args, **kwargs)
    
    chunks = [(func1d, effective_axis, sub_arr, args, kwargs)
              for sub_arr in np.array_split(arr, mp.cpu_count())]
    
    pool = mp.Pool(cores)
    individual_results = pool.map(unpack, chunks)
    # Freeing the workers:
    pool.close()
    pool.join()
    
    return np.concatenate(individual_results)


def inc_stack(small, gamma, snap, outdir, prefix=''):
    outnames_base = ['small', 'gamma', 'snap']
    outnames = [os.path.join(outdir, prefix + x) + '.tif' for x in outnames_base]
    
    if all([os.path.isfile(x) for x in outnames]):
        return outnames
    
    # set SMALL product nodata GeoTiff value
    with Raster(small)[0:100, 0:100] as ras:
        if ras.nodata is None:
            print('setting nodata value for SMALL product')
            mat = ras.matrix()
            nodata = float(mat[0, 0])
            ras2 = gdal.Open(small, GA_Update)
            ras2.GetRasterBand(1).SetNoDataValue(nodata)
            ras2 = None
    
    tmpdir = os.path.join(outdir, 'tmp')
    if not os.path.isdir(tmpdir):
        os.makedirs(tmpdir)
    
    small_edit = os.path.join(tmpdir, os.path.basename(small).replace('.tif', '_edit.tif'))
    if not os.path.isfile(small_edit):
        print('reducing resolution of SMALL product')
        gdal_translate(small, small_edit,
                       options={'xRes': 90, 'yRes': 90, 'resampleAlg': 'average',
                                'format': 'GTiff'})
    
    # subtract 90 degrees from SMALL product
    small_out = outnames[0]
    if not os.path.isfile(small_out):
        print('subtracting 90 degrees from SMALL product')
        with Raster(small_edit) as ras:
            mat = ras.matrix() - 90
            ras.assign(mat, 0)
            print('creating {}'.format(small_out))
            ras.write(small_out, format='GTiff', nodata=-99)
    
    # set GAMMA product nodata value
    with Raster(gamma) as ras:
        if ras.nodata != 0:
            print('setting nodata value of GAMMA product')
            ras2 = gdal.Open(gamma, GA_Update)
            ras2.GetRasterBand(1).SetNoDataValue(0)
            ras2 = None
    
    # convert GAMMA product from radians to degrees
    gamma_deg = os.path.join(tmpdir, os.path.basename(gamma).replace('.tif', '_deg.tif'))
    if not os.path.isfile(gamma_deg):
        print('converting GAMMA product from radians to degrees')
        with Raster(gamma) as ras:
            mat = np.rad2deg(ras.matrix())
            ras.assign(mat, 0)
            ras.write(gamma_deg, format='GTiff')
    gamma = gamma_deg
    
    # use extent of SMALL product as reference
    ext = Raster(small_out).bbox().extent
    
    # create new directory for the stacked files
    if not os.path.isdir(outdir):
        os.makedirs(outdir)
    
    # warp the products to their common extent
    warp_opts = {'options': ['-q'], 'format': 'GTiff', 'multithread': True,
                 'outputBounds': (ext['xmin'], ext['ymin'], ext['xmax'], ext['ymax']),
                 'dstNodata': -99, 'xRes': 90, 'yRes': 90, 'resampleAlg': 'bilinear',
                 'dstSRS': 'EPSG:32632'}
    
    for i, item in enumerate([gamma, snap]):
        outfile = outnames[i + 1]
        if not os.path.isfile(outfile):
            print('creating {}'.format(outfile))
            gdalwarp(src=item, dst=outfile, options=warp_opts)
    shutil.rmtree(tmpdir)
    return outnames


def commonextent(*args):
    """
    compute the common extent of multiple extent dictionaries
    
    Parameters
    ----------
    args: dict
        an extent dictionary as returned by :meth:`spatialist.Vector.extent`

    Returns
    -------
    dict:
        the common extent
    """
    ext_new = {}
    for ext in args:
        if len(ext_new.keys()) == 0:
            ext_new = ext
        else:
            for key in ['xmin', 'ymin']:
                if ext[key] > ext_new[key]:
                    ext_new[key] = ext[key]
            for key in ['xmax', 'ymax']:
                if ext[key] < ext_new[key]:
                    ext_new[key] = ext[key]
    return ext_new
