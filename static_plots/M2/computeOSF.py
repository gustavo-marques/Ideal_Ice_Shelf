import numpy
import matplotlib.pyplot as plt
import os


def computeOSF(v, vMask, dx, zInterface, zInterfaceOut, plot=False, yPlot=None,
               zPlot=None, tIndex=None):
    """
    Compute the overturning streamfunction.

    Arguments:
    v -- velocity in x (north) direction with shape (nz, ny-1, nx)

    vMask -- mask of where u is valid with the same shape as v

    dx -- the size of each cell in the x (east-west) with the same shape as v

    zInterface -- the location of interfaces between vertical layers at the
                  same horizontal locations as v, with shape (nz+1, ny-1, nx)

    zInterfaceOut -- the location of interfaces between vertical layers on the
                     output (z-level) grid, with shape nzOut+1

    Keword arguments:
    plot -- if True, plots a section of the velocity before and after
            interpolation to the output grid as well as the zonally integrated
            meridional transport and the overturning streamfunction

    yPLot -- locations of cell centers in y, used to plot the results, with
             shape nx

    zPLot -- locations of layer interfaces as horizontal cell centers, used to
             plot the results, with shape (nz+1, ny, nx)

    tIndex -- the time index, used to determine a file name for plot images

    Returns:
        osfOut -- a masked array containing the overturning streamfunciton on
                  the ISOMIP+ grid with shape (nzOut, ny)
    """

    nz = v.shape[0]
    ny = v.shape[1]+1
    nx = v.shape[2]

    nzOut = len(zInterfaceOut)-1

    vTransportOut = numpy.zeros((nzOut, ny-1, nx))
    vTransportMask = numpy.zeros((nzOut, ny-1, nx), int)
    # now, we need to conservatively interpolate uTransport to a z-level grid
    for zIndexOut in range(nzOut):
        for zIndexIn in range(nz):
            # find the overlap (if any) between this layer on the input and
            # output grids
            zTop = numpy.minimum(zInterface[zIndexIn, :, :],
                                 zInterfaceOut[zIndexOut])
            zBot = numpy.maximum(zInterface[zIndexIn+1, :, :],
                                 zInterfaceOut[zIndexOut+1])

            mask = numpy.logical_and(zTop > zBot, vMask[zIndexIn, :, :])
            area = mask*(zTop-zBot)*dx[zIndexIn, :, :]
            vTransportOut[zIndexOut, :, :] += area*v[zIndexIn, :, :]
            vTransportMask[zIndexOut, :, :] += numpy.array(mask, int)

    # the meridional transport is the sum of the transport in individual cells
    # across the x axis (axis=2)
    vTransportMerid = numpy.sum(vTransportOut, axis=2)
    # vTransportMask and uTransportMeridMask actually contain the number of
    # grid cells on the input grid that contributed to each point in
    # vTransportOut and vTransportMerid, respectively
    vTransportMeridMask = numpy.sum(vTransportMask, axis=2)

    # the overturning streamfunction is the cumulative vertical sum of the
    # meridional transport
    osf = numpy.zeros((nzOut+1, ny-1))
    osfMask = numpy.zeros(osf.shape, bool)
    osfMask[1:, :] = vTransportMeridMask > 0

    osf[1:, :] = numpy.cumsum(vTransportMerid, axis=0)

    # ISOMIP+ output grid is at cell centers, whereas osf is at u points
    # horizontally and w points vertically, so it needs to be interpolated.

    # first average in z
    osf_v = osf[0:-1, :]*osfMask[0:-1, :] + osf[1:, :]*osfMask[1:, :]
    weight_v = 1.0*osfMask[0:-1, :] + 1.0*osfMask[1:, :]

    # then eaverage in x keeping in mind that ISOMIP+ has one extra cell in
    # each direction than the ROMS output
    osfOut = numpy.zeros((nzOut, ny))
    weight = numpy.zeros((nzOut, ny))
    osfOut[:, 0:-1] += osf_v
    osfOut[:, 1:] += osf_v
    weight[:, 0:-1] += weight_v
    weight[:, 1:] += weight_v

    mask = weight > 0.
    osfOut[mask] /= weight[mask]
    osfOut = numpy.ma.masked_array(osfOut, mask=numpy.logical_not(mask))

    if plot:

        if not os.path.exists('plots'):
            os.makedirs('plots')

        vOut = numpy.ma.masked_all((nzOut, ny-1, nx))
        for zIndexOut in range(nzOut):
            zTop = numpy.minimum(zInterface[0, :, :],
                                 zInterfaceOut[zIndexOut])
            zBot = numpy.maximum(zInterface[-1, :, :],
                                 zInterfaceOut[zIndexOut+1])
            areaOut = dx[0, :, :]*(zTop-zBot)
            mask = areaOut > 0.

            vSlice = numpy.zeros((ny-1, nx))
            vSlice[mask] = vTransportOut[zIndexOut, :, :][mask]/areaOut[mask]
            vOut[zIndexOut, :, :] = vSlice

        xIndex = nx/2

        YIn = yPlot.reshape((1, ny))
        ZIn = zPlot[:, :, xIndex]

        plt.close('all')

        plt.figure(1)
        mask = numpy.logical_not(vMask[:, :, xIndex])
        plt.pcolor(1e-3*YIn, ZIn, numpy.ma.masked_array(v[:,:, xIndex],
                                                        mask=mask))
        plt.colorbar()
        plt.title('v (m/s) on the input grid')
        plt.xlabel('y (km)')
        plt.ylabel('z (m)')
        plt.axis('tight')
        plt.savefig('plots/vIn_{:03d}.png'.format(tIndex))

        (YOut, ZOut) = numpy.meshgrid(yPlot, zInterfaceOut, indexing='xy')

        plt.figure(2)
        mask = vOut[:, :, xIndex] == 0.
        plt.pcolor(1e-3*YOut, ZOut, numpy.ma.masked_array(vOut[:, : , xIndex],
                                                          mask=mask))
        plt.colorbar()
        plt.title('v (m/s) on the output grid')
        plt.xlabel('y (km)')
        plt.ylabel('z (m)')
        plt.axis('tight')
        plt.savefig('plots/vOut_{:03d}.png'.format(tIndex))

        plt.figure(3)
        mask = vTransportMeridMask == 0
        plt.pcolor(1e-3*YOut, ZOut, 1e-6*numpy.ma.masked_array(vTransportMerid,
                                                               mask=mask))
        plt.colorbar()
        plt.title('meridional transport (Sv) on the output grid')
        plt.xlabel('y (km)')
        plt.ylabel('z (m)')
        plt.axis('tight')
        plt.savefig('plots/meridTransOut_{:03d}.png'.format(tIndex))

        plt.figure(4)
        plt.pcolor(1e-3*YOut, ZOut, vTransportMeridMask)
        plt.colorbar()
        plt.title('input grid points contributing to meridional tranport')
        plt.xlabel('y (km)')
        plt.ylabel('z (m)')
        plt.axis('tight')
        plt.savefig('plots/meridTransCount_{:03d}.png'.format(tIndex))

        zOut = 0.5*(zInterfaceOut[0:-1] + zInterfaceOut[1:])
        (YOut, ZOut) = numpy.meshgrid(yPlot, zOut, indexing='xy')

        plt.figure(5)
        mask = numpy.logical_not(osfMask[1:-1, :])
        plt.pcolor(1e-3*YOut, ZOut, 1e-6*numpy.ma.masked_array(osf[1:-1, :],
                                                               mask=mask))
        plt.colorbar()
        plt.title('meridional overturning streamfunction (Sv) on the output '
                  'grid')
        plt.xlabel('y (km)')
        plt.ylabel('z (m)')
        plt.axis('tight')
        plt.savefig('plots/osf_{:03d}.png'.format(tIndex))

        dy = yPlot[1]-yPlot[0]
        y0 = yPlot[0]-0.5*dy
        (YOut, ZOut) = numpy.meshgrid(y0 + dy*(numpy.arange(ny+1)),
                                      zInterfaceOut,
                                      indexing='xy')

        plt.figure(6)
        plt.pcolor(1e-3*YOut, ZOut, 1e-6*osfOut)
        plt.colorbar()
        plt.title('meridional overturning streamfunction (Sv) ISOMIP+ grid')
        plt.xlabel('y (km)')
        plt.ylabel('z (m)')
        plt.axis('tight')
        plt.savefig('plots/osfInterp_{:03d}.png'.format(tIndex))

    return osfOut
