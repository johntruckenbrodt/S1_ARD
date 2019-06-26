import os
import re
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
from matplotlib.lines import Line2D
from matplotlib.offsetbox import AnchoredText
from matplotlib import patches

from osgeo import gdal, ogr
from osgeo.gdalconst import GA_Update

from spatialist import haversine, Raster, Vector, crsConvert, gdalwarp, gdal_translate, intersect
from spatialist.ancillary import finder


def scatter(x, y, z=None, xlab='', ylab='', title='', nsamples=1000, mask=None,
            measures=None, regline=False, o2o=False, denscol=False, grid=False,
            xlim=None, ylim=None, sort_z=False, legend=False,
            regline_label='regression', o2o_label='1-to-1'):
    """
    general function for creating scatter plots.
    
    Parameters
    ----------
    x: numpy.ndarray
        dataset I
    y: numpy.ndarray
        dataset II
    z: numpy.ndarray
        dataset III for coloring the data points; overrides parameter `denscol`
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
         - `cv_x`, `cv_y`: the coefficient of variation of either `x` or `y`
         - `mean_x`, `mean_y`: the mean value of either `x` or `y`
    regline: bool
        draw a linear regression line?
    o2o: bool
        draw a data one-to-one line?
    denscol: bool
        color the points by Gaussian density?; overridden by parameter `z`
    grid: bool
        add a mesh grid to the plot?
    xlim: tuple
        the x-axis limits
    ylim: tuple
        the y-axis limits
    sort_z: bool
        if `z` is not None, sort its values so that points with high `z`
        values are plotted last?
    legend: bool
        add a legend for the regression line and one-to-one line if they exist?
    regline_label: str
        the legend label for the regression line
    o2o_label: str
        the legend label for the one-to-one line

    Returns
    -------

    """
    nanmask = np.isfinite(x) & np.isfinite(y)
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
        sign = '+' if m > 0 else '-'
        fields.append('y = {:.2f} {} {:.2f} * x'.format(b, sign, abs(m)))
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
    if 'mean_x' in measures:
        mean = np.mean(x)
        fields.append(r'$\mu$(x) = {:.2f}'.format(mean))
    if 'mean_y' in measures:
        mean = np.mean(y)
        fields.append(r'$\mu$(y) = {:.2f}'.format(mean))
    
    text = '\n'.join(fields)
    
    if z is not None:
        z = z.flatten()[sample_ids]
        denscol = False
        
        # # Sort the points by z, so that points with high z values are plotted last
        if sort_z:
            idx = z.argsort()
            x, y, z = x[idx], y[idx], z[idx]
    
    if denscol:
        # # Calculate the point density
        xy = np.vstack([x, y])
        z = gaussian_kde(xy)(xy)
        
        # # Sort the points by density, so that the densest points are plotted last
        idx = z.argsort()
        x, y, z = x[idx], y[idx], z[idx]
    
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
        plt.plot(x_new, ffit(x_new), color='red', zorder=2, label=regline_label)
    if o2o:
        ax = plt.gca()
        line = Line2D([0, 1], [0, 1], color='black', zorder=1, label=o2o_label)
        transform = ax.transAxes
        line.set_transform(transform)
        ax.add_line(line)
    plt.title(title)
    plt.xlabel(xlab)
    plt.ylabel(ylab)
    if legend:
        ax = plt.gca()
        handles, labels = ax.get_legend_handles_labels()
        if len(handles) > 0:
            plt.legend(handles, labels, loc='upper left')
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
    """
    compute the aspect of a DEM.
    
    Parameters
    ----------
    img: numpy.ndarray
        the DEM array

    Returns
    -------
    numpy.ndarray
        the computed aspect array
    """
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


def dem_distribution(slope, aspect, head_angle, inc_angle, look_dir='right',
                     nsamples=1000, title='', mask=None):
    """
    create a polar slope-aspect DEM plot superimposed with the area visible to a SAR sensor.
    
    Parameters
    ----------
    slope: numpy.ndarray
    aspect: numpy.ndarray
    head_angle: float
        the SAR sensor heading
    inc_angle: float
        the SAR sensor's incident angle
    look_dir: str
        the SAR sensor look direction; either `left` or `right`
    nsamples: int
        the number of samples to select from the `slope` and `aspect` arrays
        using function :func:`~S1_ARD.util.sampler`
    title: str
        the plot's title
    mask: numpy.ndarray
        an additional binary array to mask the slope and aspect values

    Returns
    -------
    
    See Also
    --------
    visible_sar_angle_map
    """
    nanmask = np.isfinite(slope) & np.isfinite(aspect)
    if mask is None:
        mask = nanmask
    else:
        mask = mask & nanmask
    
    sample_ids = sampler(mask, nsamples)
    
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
    compute the slope of a DEM.
    
    Parameters
    ----------
    img: numpy.ndarray
        the input DEM
    xres_m: int or float
        the x resolution of the DEM in same units as the height values
    yres_m: int or float
        the y resolution of the DEM in same units as the height values

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
    """
    compute the spatial resolution in meters for a DEM with WGS84 degree coordinates.
    
    Parameters
    ----------
    demfile: str
        the DEM file

    Returns
    -------
    tuple
        (posting_east, posting_north)
    
    See Also
    --------
    spatialist.auxil.haversine
    """
    with Raster(demfile) as ras:
        res_lon, res_lat = ras.res
        lat = (ras.geo['ymin'] + ras.geo['ymax']) / 2
        lon = (ras.geo['xmin'] + ras.geo['xmax']) / 2
    post_north = haversine(lat, lon, lat + res_lat, lon)
    post_east = haversine(lat, lon, lat, lon + res_lon)
    return post_east, post_north


def sampler(nanmask, nsamples=None, seed=42):
    """
    central function to select random samples from arrays.
    
    Parameters
    ----------
    nanmask: numpy.ndarray
        a mask to limit the sample selection
    nsamples: int
        the number of samples to select
    seed: int
        seed used to initialize the pseudo-random number generator

    Returns
    -------
    numpy.ndarray
        the generated random samples
    
    See Also
    --------
    numpy.random.seed
    numpy.random.choice
    """
    indices = np.where(nanmask.flatten())[0]
    samplesize = min(indices.size, nsamples) if nsamples is not None else indices.size
    np.random.seed(seed)
    sample_ids = np.random.choice(a=indices, size=samplesize, replace=False)
    return sample_ids


def visible_sar_angle_map(head_angle, inc_angle, look_dir='right'):
    """
    create a SAR sensor slope-aspect visibility mask;
    used by :func:`~S1_ARD.util.dem_distribution`.
    
    Parameters
    ----------
    head_angle: float
        the SAR sensor heading
    inc_angle: float
        the SAR sensor's incident angle
    look_dir: str
        the SAR sensor look direction; either `left` or `right`

    Returns
    -------
    numpy.ndarray
        the binary map with aspect-slope coordinates
    """
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
    """
    convert a well-known text string geometry to a shapefile.
    
    Parameters
    ----------
    wkt: str
        the well-known text description
    srs: int, str
        the spatial reference system; see :func:`spatialist.auxil.crsConvert` for options.
    outname: str
        the name of the shapefile to write

    Returns
    -------

    """
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
    Like :func:`numpy.apply_along_axis()`, but takes advantage of multiple cores.
    Adapted from `here <https://stackoverflow.com/questions/45526700/
    easy-parallelization-of-numpy-apply-along-axis>`_.
    
    Parameters
    ----------
    func1d: function
        the function to be applied
    axis: int
        the axis along which to apply `func1d`
    arr: numpy.ndarray
        the input array
    cores: int
        the number of parallel cores
    args: any
        Additional arguments to `func1d`.
    kwargs: any
        Additional named arguments to `func1d`.

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
    compute the common extent of multiple extent dictionaries.
    
    Parameters
    ----------
    args: dict
        an extent dictionary, see e.g. :attr:`spatialist.vector.Vector.extent`

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


def clc_legend(filename):
    """
    read clc meta data from the dedicated CSV file available
    `here <https://www.eea.europa.eu/data-and-maps/data/
    corine-land-cover-3/corine-land-cover-classes-and/clc_legend.csv>`_.
    
    Parameters
    ----------
    filename: str
        the CSV file to be read

    Returns
    -------
    dict
        the CSV values in a dictionary
    """
    pattern = r'([0-9]{1,2}),' \
              r'([0-9]{3}),' \
              r'([a-zA-Z ]*),' \
              r'(["]*[a-zA-Z ,-/]*["]*),' \
              r'(["]*[a-zA-Z ,-/]*["]*),' \
              r'([0-9-]{11})\n'
    
    with open(filename, 'r') as clc:
        keys = clc.readline().strip().split(',')
        out = dict(zip(keys, [[] for x in keys]))
        for line in clc:
            match = re.search(pattern, line)
            if match:
                items = match.groups()
                for i in range(0, len(keys)):
                    key = keys[i]
                    if key == 'RGB':
                        value = tuple([int(x) / 255. for x in items[i].split('-')])
                    else:
                        value = items[i].replace('"', '')
                    out[key].append(value)
    return out


def clc_prep(clc, reference, outname):
    """
    resample and crop the corine product to the resolution and extent of a reference image.
    
    Parameters
    ----------
    clc: str
        the name of the CLC input file
    reference: str
        the name of the reference file
    outname: str
        the named of the output image

    Returns
    -------

    """
    with Raster(reference).bbox() as box_ras:
        with Raster(clc).bbox() as box_clc:
            if intersect(box_ras, box_clc) is None:
                print('no intersect')
                return
    
    if not os.path.isfile(outname):
        with Raster(reference) as ras:
            ref_crs = ras.projection
            xres, yres = ras.res
            with ras.bbox() as box:
                ref_ext = box.extent
        
        outputBounds = (ref_ext['xmin'], ref_ext['ymin'], ref_ext['xmax'], ref_ext['ymax'])
        
        gdalwarp_opt = {'format': 'GTiff', 'outputBounds': outputBounds, 'multithread': True,
                        'xRes': xres, 'yRes': yres, 'dstSRS': ref_crs, 'resampleAlg': 'mode'}
        
        gdalwarp(src=clc, dst=outname, options=gdalwarp_opt)
    else:
        print('outfile already exists')


def clc_prepare(reference, outdir, source):
    """
    create a CLC subset resampled to a reference image.
    
    Parameters
    ----------
    reference: str
        the reference file with the target CRS and extent
    outdir: str
        the directory to write the new file to;
        new files are named clc{index}.tif, e.g. clc1.tif.
    source: str
        the original product to be subsetted

    Returns
    -------
    str
        the name of the file written to `outdir`
    """
    with Raster(reference) as ras:
        xRes, yRes = ras.res
        epsg = ras.epsg
        ext = ras.extent
    
    #########################################################################
    warp_opts = {'options': ['-q'], 'format': 'GTiff', 'multithread': True,
                 'dstNodata': -99, 'resampleAlg': 'mode'}
    
    if not os.path.isdir(outdir):
        os.makedirs(outdir)
    
    clc_subs = finder(outdir, ['clc[0-9].tif'], regex=True)
    
    match = False
    if len(clc_subs) > 0:
        for j, sub in enumerate(clc_subs):
            with Raster(sub) as ras:
                if ras.extent == ext:
                    clc_sub = sub
                    match = True
    if not match:
        clc_sub = os.path.join(outdir, 'clc{}.tif'.format(len(clc_subs)))
        print('creating', clc_sub)
        warp_opts['dstSRS'] = 'EPSG:{}'.format(epsg)
        warp_opts['xRes'] = xRes
        warp_opts['yRes'] = yRes
        warp_opts['outputBounds'] = (ext['xmin'], ext['ymin'],
                                     ext['xmax'], ext['ymax'])
        gdalwarp(src=source, dst=clc_sub, options=warp_opts)
    return clc_sub


def uzh_prepare(reference, outdir, source):
    """
    create an UZH incident angle subset resampled to a reference image.
    
    Parameters
    ----------
    reference: str
        the reference file with the target extent
    outdir: str
        the directory to write the new file to;
        new files are named uzh_{epsg}_{index}.tif, e.g. uzh_4326_1.tif.
    source: str
        the original product to be subsetted

    Returns
    -------
    numpy.ndarray
        the content of the file written to `outdir`
    """
    with Raster(reference) as ras:
        xRes, yRes = ras.res
        epsg = ras.epsg
        ext = ras.extent
    
    warp_opts = {'options': ['-q'], 'format': 'GTiff', 'multithread': True,
                 'dstNodata': -99, 'resampleAlg': 'bilinear'}
    
    if not os.path.isdir(outdir):
        os.makedirs(outdir)
    
    # find existing files
    uzh_subs = finder(outdir, ['uzh_[0-9]{4,5}_[0-9].tif'], regex=True)
    
    # check if any of the existing files matches the extent of the reference
    match = False
    if len(uzh_subs) > 0:
        for j, sub in enumerate(uzh_subs):
            with Raster(sub) as ras:
                if ras.extent == ext:
                    uzh_sub = sub
                    match = True
    if not match:
        with Raster(source) as ras:
            if ras.epsg != epsg:
                raise RuntimeError('CRS mismatch')
        
        basename = 'uzh_{}_{}.tif'.format(epsg, len(uzh_subs))
        uzh_sub = os.path.join(outdir, basename)
        print('creating', uzh_sub)
        warp_opts['dstSRS'] = 'EPSG:{}'.format(epsg)
        warp_opts['xRes'] = xRes
        warp_opts['yRes'] = yRes
        warp_opts['outputBounds'] = (ext['xmin'], ext['ymin'],
                                     ext['xmax'], ext['ymax'])
        gdalwarp(src=source, dst=uzh_sub, options=warp_opts)
    return uzh_sub


def dev_max(arr):
    """
    compute the maximum deviation from the median of all array
    values and the corresponding ID.
    
    Parameters
    ----------
    arr: numpy.ndarray
        the 1D array
    Returns
    -------
    tuple
        (maximum deviation, ID)
    """
    if len(arr[~np.isnan(arr)]) > 0:
        med = np.nanmedian(arr)
        absdev = abs(arr - med)
        maxdev = np.nanmax(absdev)
        if maxdev > 0:
            maxdev_id = np.nanargmax(absdev)
        else:
            maxdev_id = np.nan
        return maxdev, maxdev_id
    else:
        return np.nan, np.nan


def extent2patch(extent, edgecolor='r'):
    """
    create a matplotlib rectangle patch from an extent dictionary.
    
    Parameters
    ----------
    extent: dict
        an extent dictionary, see e.g. :attr:`spatialist.vector.Vector.extent`
    edgecolor: str
        the edge color of the path

    Returns
    -------
    matplotlib.patches.Rectangle
    """
    xdiff = extent['xmax'] - extent['xmin']
    ydiff = extent['ymax'] - extent['ymin']
    rect = patches.Rectangle((extent['xmin'], extent['ymin']),
                             xdiff, ydiff, linewidth=2,
                             edgecolor=edgecolor, facecolor='none')
    return rect
