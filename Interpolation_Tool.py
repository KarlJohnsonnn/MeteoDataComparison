import numpy as np

from scipy import interpolate

def Interpolate_2D(x1, y1, z1, x2, y2, method):
    len_x1 = len(x1)
    len_x2 = len(x2)
    len_y1 = len(y1)
    len_y2 = len(y2)

    if method in ['linear', 'quadratic', 'cubic', 'quintic']:
        fcn = interpolate.interp2d(x1, y1, z1, kind=method)

        interp_z = fcn(x2, y2)
        interp_z = np.ma.masked_equal(interp_z, 0.0)
        interp_z = np.ma.masked_less_equal(interp_z, -100.0)
        interp_z = np.ma.masked_invalid(interp_z)


    elif method in ['NearestNeighbour']:

        coord = np.empty((len_x1 * len_y1, 2))
        values = np.empty((len_x1 * len_y1, 1))
        cnt = 0
        for i in range(len_x1):
            for j in range(len_y1):
                coord[cnt, 0] = x1[i]
                coord[cnt, 1] = y1[j]
                values[cnt] = z1[j, i]
                cnt += 1
                print('in copy loop  ', i, j, x1[i], y1[j], z1[j, i])

        values = np.ma.masked_invalid(values)
        fcn = interpolate.NearestNDInterpolator(coord, values)

        interp_z = np.empty((len_x2 * len_y2))
        cnt = 0
        for xi in x2:
            for yi in y2:
                interp_z[cnt] = fcn(xi, yi)
                cnt += 1
                print('in fcn loop  ', cnt, xi, yi)

        interp_z = np.ma.masked_equal(interp_z, 0.0)
        interp_z = np.ma.masked_less_equal(interp_z, -100.0)
        interp_z = np.ma.masked_invalid(interp_z)
        interp_z = np.reshape(interp_z, (len_x2, len_y2)).T

    return interp_z

def Interpolate_2D_neu(ds1, z_inter, ds2, method):
    len_x1 = ds1.n_time;    len_x2 = ds2.n_time
    len_y1 = ds1.n_height;  len_y2 = ds2.n_height

    x1 = ds1.time;    x2 = ds2.time
    y1 = ds1.height;  y2 = ds2.height

    z1 = z_inter

    if method in ['linear', 'quadratic', 'cubic', 'quintic']:
        fcn = interpolate.interp2d(x1, y1, z1, kind=method)

        interp_z = fcn(x2, y2)
        interp_z = np.ma.masked_equal(interp_z, 0.0)
        interp_z = np.ma.masked_less_equal(interp_z, -100.0)
        interp_z = np.ma.masked_invalid(interp_z)


    elif method in ['NearestNeighbour']:

        coord = np.empty((len_x1 * len_y1, 2))
        values = np.empty((len_x1 * len_y1, 1))
        cnt = 0
        for i in range(len_x1):
            for j in range(len_y1):
                coord[cnt, 0] = x1[i]
                coord[cnt, 1] = y1[j]
                values[cnt] = z1[j, i]
                cnt += 1
                print('in copy loop  ', i, j, x1[i], y1[j], z1[j, i])

        values = np.ma.masked_invalid(values)
        fcn = interpolate.NearestNDInterpolator(coord, values)

        interp_z = np.empty((len_x2 * len_y2))
        cnt = 0
        for xi in x2:
            for yi in y2:
                interp_z[cnt] = fcn(xi, yi)
                cnt += 1
                print('in fcn loop  ', cnt, xi, yi)

        interp_z = np.ma.masked_equal(interp_z, 0.0)
        interp_z = np.ma.masked_less_equal(interp_z, -100.0)
        interp_z = np.ma.masked_invalid(interp_z)
        interp_z = np.reshape(interp_z, (len_x2, len_y2)).T

    return interp_z

def interpolate_data(x, y, xnew, method):
    # create a callable function from the actual data
    fnew = interpolate.interp1d(x, y, kind=method)

    # calculate the interpolation dataset
    ynew = fnew(xnew)

    # mask values ( minInterpol=minAcutalData, etc. )
    ynew = np.ma.masked_greater_equal(ynew, y.max())
    ynew = np.ma.masked_less_equal(ynew, y.min())

    return ynew