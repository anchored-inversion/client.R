#include <math.h>
#include <R.h>

/*
 * 2D travel time calculation for wave propagation.
 * Implemented based on the finite-difference method described in
 *   Vidale, 1988, Finite-difference calculation of travel times, in
 *   Bulletin of the Seismological Society of American, vol 78, no 6.
 * and
 *   Vidale, 1990, Finite-difference calculation of traveltimes in three
 *   dimensions, in Geophysics, vol 55, no 5.
 */

int imax(int x, int y)
{
    return((x > y) ? x : y);
}


int imin(int x, int y)
{
    return((x < y) ? x : y);
}


void travel_time_2d_expand(
    const double * const slowness_inner,
    const double * const slowness_outer,
    const double * const ttime_inner,
    double *ttime_outer,
    int *index,
    const int n,
    const double h,
    int *nsingular,
    int *nneg)
{
    double t0, t1, t2, t3, s;
    int idx;
    double ftemp;

    *nsingular = 0;
    *nneg = 0;

    for (int i = 0; i < n; i++)
    {
        ttime_outer[i] = ttime_inner[i];
        index[i] = i;
    }

    rsort_with_index(ttime_outer, index, n);
        // Now 'index' contains the indices
        // of points in the order of their visit,
        // from low 'ttime_inner' to high 'ttime_inner'.
        // The array 'ttime_inner' is not changed.

    for (int i = 0; i < n; i++)
        ttime_outer[i] = -1.;


    for (int i = 0; i < n; i++)
    {
        idx = index[i];

        if (idx == 0)
        {
            if (ttime_outer[idx+1] < 0.)
            {
                (*nsingular) ++;

                t0 = ttime_inner[idx];
                t1 = ttime_inner[idx+1];
                s = (slowness_inner[idx] + slowness_inner[idx+1]
                        + slowness_outer[idx]) / 3.;
                ftemp = h*s*h*s - (t1-t0) * (t1-t0);
                if (ftemp > 0)
                    ttime_outer[idx] = t0 + sqrt(ftemp);
                else
                {
                    ttime_outer[idx] = t0;
                    (*nneg) ++;
                }
            } else
            {
                t0 = ttime_inner[idx+1];
                t1 = ttime_inner[idx];
                t2 = ttime_outer[idx+1];
                s = (slowness_inner[idx+1] + slowness_inner[idx]
                    + slowness_outer[idx+1] + slowness_outer[idx]) / 4.;
                ftemp = 2. * (h*s) * (h*s) - (t2-t1) * (t2-t1);
                if (ftemp > 0)
                    ttime_outer[idx] = t0 + sqrt(ftemp);
                else
                {
                    ttime_outer[idx] = t0;
                    (*nneg) ++;
                }
            }
        } else if (idx == n - 1)
        {
            if (ttime_outer[idx-1] < 0.)
            {
                (*nsingular) ++;

                t0 = ttime_inner[idx-1];
                t1 = ttime_inner[idx];
                s = (slowness_inner[idx-1] + slowness_inner[idx]
                        + slowness_outer[idx]) / 3.;
                ftemp = h*s*h*s - (t1-t0) * (t1-t0);
                if (ftemp > 0)
                    ttime_outer[idx] = t0 + sqrt(ftemp);
                else
                {
                    ttime_outer[idx] = t0;
                    (*nneg) ++;
                }
            } else
            {
                t0 = ttime_inner[idx-1];
                t1 = ttime_inner[idx];
                t2 = ttime_outer[idx-1];
                s = (slowness_inner[idx-1] + slowness_inner[idx] +
                    slowness_outer[idx-1] + slowness_outer[idx]) / 4.;
                ftemp = 2. * (h*s) * (h*s) - (t2-t1) * (t2-t1);
                if (ftemp > 0)
                    ttime_outer[idx] = t0 + sqrt(ftemp);
                else
                {
                    ttime_outer[idx] = t0;
                    (*nneg) ++;
                }
            }
        } else
        {
            if (ttime_outer[idx-1] < 0. && ttime_outer[idx+1] < 0.)
            {
                t0 = ttime_inner[idx];
                t1 = ttime_inner[idx-1];
                t2 = ttime_inner[idx+1];
                s = (slowness_inner[idx] + slowness_inner[idx-1] +
                    slowness_inner[idx+1] + slowness_outer[idx]) / 4.;
                ftemp = h*s*h*s - .25 * (t2 - t1) * (t2 - t1);
                if (ftemp > 0)
                    ttime_outer[idx] = t0 + sqrt(ftemp);
                else
                {
                    ttime_outer[idx] = t0;
                    (*nneg) ++;
                }
            } else
            {
                if (ttime_outer[idx-1] >= 0.)
                {
                    t0 = ttime_inner[idx-1];
                    t1 = ttime_inner[idx];
                    t2 = ttime_outer[idx-1];
                    s = (slowness_inner[idx-1] + slowness_inner[idx] +
                        slowness_outer[idx-1] + slowness_outer[idx]) / 4.;
                    ftemp = 2. * (h*s) * (h*s) - (t2-t1) * (t2-t1);
                    if (ftemp > 0)
                        ttime_outer[idx] = t0 + sqrt(ftemp);
                    else
                    {
                        ttime_outer[idx] = t0;
                        (*nneg) ++;
                    }
                }

                if (ttime_outer[idx+1] >= 0.)
                {
                    t0 = ttime_inner[idx+1];
                    t1 = ttime_inner[idx];
                    t2 = ttime_outer[idx+1];
                    s = (slowness_inner[idx+1] + slowness_inner[idx] +
                        slowness_outer[idx+1] + slowness_outer[idx]) / 4.;
                    ftemp = 2. * (h*s) * (h*s) - (t2-t1) * (t2-t1);
                    if (ftemp > 0)
                        t3 = t0 + sqrt(ftemp);
                    else
                        t3 = t0;
                    if (ttime_outer[idx] < 0.)
                        ttime_outer[idx] = t3;
                    else
                    {
                        ttime_outer[idx] = fmin(ttime_outer[idx], t3);
                        (*nneg) ++;
                    }
                }
            }
        }
    }
}



void c_travel_time_2d(
    int *nnx, int *nny, double *hh,
    double *source,
    double *slowness,
    int *innerradius,
    double *ttime,
    int *nsingular,
    int *nneg)
{
/* *
 * The model modain is divided into _square_ grid boxes.
 * 'Grid points' are box centers.
 * Origin of the coord system is (0, 0) at the lower-left grid point.
 * 'X' goes right, 'Y' goes up.
 * nnx, nny:    Dimensions of the 2D grid.
 *              'nnx' ('nny') is number of grids along 'X' ('Y').
 * hh:          Common grid size in both 'X' and 'Y' directions.
 *              (meters)
 * source:      (x, y) coordinates of wave source.
 *              Must be within the grid domain;
 *              but need not coincide with a grid point.
 * slowness:    Array of slowness in each grid box (or at each grid point).
 *              Slowness = 1 / velocity. (seconds/meter)
 *              5-8 km/second are realistic values for seismic wave speed.
 *              Stored as if one scans an image one row at a time
 *              from left to right, starting at the bottom row
 *              (i.e. low Y to high Y).
 *              When passing in a R matrix, 1st col of the matrix
 *              contains first horizontal scan.
 *              In other words, row index of the R matrix corresponds to
 *              'x', col index corresponds to 'y'.
 *              Size of array is (nnx * nny).
 * innerradius: For the several smallest circles surrounding the source,
 *              the travel time will be calculated by directly
 *              integrating slowness, assuming a straight ray path.
 *              The most outside circel of all these has radius
 *              'innerradius'. The source itself has radius 0.
 * ttime:       Resultant array of travel time (seconds)
 *              from 'source' to each grid point.
 *              Layout of the elements corresponds to that of
 *              'slowness', that is, column major to R convention.
 *              The element layout is as 'slowness'.
 * */

    int nx = *nnx, ny = *nny;
    double h = *hh;
    double x0 = source[0];
    double y0 = source[1];
        // Source coord.
    int ix0 = (int) (x0 / h + 0.5);
    int iy0 = (int) (y0 / h + 0.5);
        // The grid box index of 'source'.

    double ftemp;
    double eps = h * 1.e-6;
    int n1, n2;

    *nsingular = 0;
    *nneg = 0;

    // Process the inner circles:
    // in these circles, travel time is calculated by assuming
    // straight ray paths and simply integrating the slowness
    // on the ray paths.


    if (*innerradius < 1)
        *innerradius = 1;

    for (int ix = imax(0, ix0 - *innerradius); ix <= imin(nx-1, ix0 + *innerradius); ix++)
    {
        for (int iy = imax(0, iy0 - *innerradius); iy <= imin(ny-1, iy0 + *innerradius); iy++)
        {
            double timetraveled = 0.;

            double dx = ix*h - source[0];
            double dy = iy*h - source[1];
            int ixnow = ix0;
            int iynow = iy0;

            if (ixnow == ix && iynow == iy)
            {
                timetraveled = sqrt(dx * dx + dy * dy) * slowness[nx * iy + ix];
                ttime[nx * iy + ix] = timetraveled;
                continue;
            }


            if (fabs(dx) < eps)
            {
                if (dy < 0)
                {
                    timetraveled += (y0 - (iy0 - .5) * h) * slowness[nx * iy0 + ix0];
                    for (int k = iy0 - 1; k > iy; k--)
                        timetraveled += h * slowness[nx * k + ix0];
                } else
                {
                    timetraveled += ((iy0 + .5) * h - y0) * slowness[nx * iy0 + ix0];
                    for (int k = iy0 + 1; k < iy; k++)
                        timetraveled += h * slowness[nx * k + ix0];
                }
                timetraveled += .5 * h * slowness[nx * iy + ix];
            } else if (fabs(dy) < eps)
            {
                if (dx < 0)
                {
                    timetraveled += (x0 - (ix0 - .5) * h) * slowness[nx * iy0 + ix0];
                    for (int k = ix0 - 1; k > ix; k--)
                        timetraveled += h * slowness[nx * iy0 + k];
                } else
                {
                    timetraveled += ((ix0 + .5) * h - x0) * slowness[nx * iy0 + ix0];
                    for (int k = ix0 + 1; k < ix; k++)
                        timetraveled += h * slowness[nx * iy0 + k];
                }
                timetraveled += .5 * h * slowness[nx * iy + ix];
            } else
            {
                double xnow = x0;
                double ynow = y0;
                int ixnext, iynext;

                while (ixnow != ix || iynow != iy)
                {
                    double deltax = (ixnow + ((dx > 0)? .5 : (-.5))) * h - xnow;
                    double deltay = (iynow + ((dy > 0)? .5 : (-.5))) * h - ynow;
                    ftemp = deltay * dx/dy;
                    if (fabs(ftemp) < fabs(deltax) - eps)
                    {
                        deltax = ftemp;
                        ixnext = ixnow;
                        iynext = ((dy > 0)? 1 : (-1)) + iynow;
                    } else if (fabs(ftemp) > fabs(deltax) + eps)
                    {
                        deltay = deltax * dy/dx;
                        ixnext = ((dx > 0)? 1 : (-1)) + ixnow;
                        iynext = iynow;
                    } else
                    {
                        ixnext = ((dx > 0)? 1 : (-1)) + ixnow;
                        iynext = ((dy > 0)? 1 : (-1)) + iynow;
                    }

                    timetraveled += sqrt(deltax*deltax + deltay*deltay)
                        * slowness[nx * iynow + ixnow];

                    xnow += deltax;
                    ynow += deltay;
                    ixnow = ixnext;
                    iynow = iynext;
                }

                ftemp = sqrt((ix*h - xnow) * (ix*h - xnow) + (iy*h - ynow) * (iy*h - ynow));
                timetraveled +=  ftemp * slowness[nx * iy + ix];
            }

            ttime[nx * iy + ix] = timetraveled;
        }
    }


    // Outside of the inner circles.

    int radius_max = imax(imax(ix0, nx - 1 - ix0), imax(iy0, ny - 1 - iy0));

    if (*innerradius >= radius_max) return;

    double *ttime_inner = (double *) malloc(sizeof(double) * imax(nx, ny));
    double *ttime_outer = (double *) malloc(sizeof(double) * imax(nx, ny));
    double *slowness_inner = (double *) malloc(sizeof(double) * imax(nx, ny));
    double *slowness_outer = (double *) malloc(sizeof(double) * imax(nx, ny));
    int *sindex = (int *) malloc(sizeof(int) * imax(nx, ny));

    int ix, iy;
    int k, idx1, idx2, idx0;

    for (int radius = (*innerradius) + 1; radius <= radius_max; radius++)
    {
        // Expand to left.
        ix = ix0 - radius;
        if (ix >= 0)
        {
            idx1 = imax(0, iy0 - radius + 1);
                // Index of lower Y.
            idx2 = imin(ny-1, iy0 + radius - 1);
                // Index of upper Y.
            idx0 = nx * idx1 + ix + 1;
                // Index (in array 'slowness') of the lowest grid point
                // in inner circle.
            k = 0;
            for (iy = idx1; iy <= idx2; iy++)
            {
                slowness_inner[k] = slowness[idx0];
                ttime_inner[k] = ttime[idx0];
                slowness_outer[k] = slowness[idx0 - 1];
                k++;
                idx0 += nx;
            }
            travel_time_2d_expand(slowness_inner, slowness_outer,
                    ttime_inner, ttime_outer, sindex, k, h, &n1, &n2);

            *nsingular += n1;
            *nneg += n2;

            k = 0;
            idx0 = nx * idx1 + ix;
                // Index (in array 'slowness') of the lowest grid point
                // in outer circle.
            for (iy = idx1; iy <= idx2; iy++)
            {
                ttime[idx0] = ttime_outer[k];
                idx0 += nx;
                k++;
            }
        }

        // Expand to right.
        ix = ix0 + radius;
        if (ix < nx)
        {
            idx1 = imax(0, iy0 - radius + 1);
                // Index of lower Y.
            idx2 = imin(ny-1, iy0 + radius - 1);
                // Index of upper Y.
            idx0 = nx * idx1 + ix - 1;
                // Index (in array 'slowness') of the lowest grid point
                // in inner circle.
            k = 0;
            for (iy = idx1; iy <= idx2; iy++)
            {
                slowness_inner[k] = slowness[idx0];
                ttime_inner[k] = ttime[idx0];
                slowness_outer[k] = slowness[idx0 + 1];
                k++;
                idx0 += nx;
            }

            travel_time_2d_expand(slowness_inner, slowness_outer,
                    ttime_inner, ttime_outer, sindex, k, h, &n1, &n2);

            *nsingular += n1;
            *nneg += n2;

            k = 0;
            idx0 = nx * idx1 + ix;
                // Index (in array 'slowness') of the lowest grid point
                // in outer circle.
            for (iy = idx1; iy <= idx2; iy++)
            {
                ttime[idx0] = ttime_outer[k];
                idx0 += nx;
                k++;
            }
        }

        // Expand to bottom.
        iy = iy0 - radius;
        if (iy >= 0)
        {
            idx1 = imax(0, ix0 - radius + 1);
                // Index of left X.
            idx2 = imin(nx-1, ix0 + radius - 1);
                // Index of right X.
            idx0 = nx * (iy + 1) + idx1;
                // Index (in array 'slowness') of the left grid point
                // in inner circle.
            k = 0;
            for (ix = idx1; ix <= idx2; ix++)
            {
                slowness_inner[k] = slowness[idx0];
                ttime_inner[k] = ttime[idx0];
                slowness_outer[k] = slowness[idx0 - nx];
                idx0++;
                k++;
            }

            travel_time_2d_expand(slowness_inner, slowness_outer,
                    ttime_inner, ttime_outer, sindex, k, h, &n1, &n2);

            *nsingular += n1;
            *nneg += n2;

            k = 0;
            idx0 = nx * iy + idx1;
                // Index (in array 'slowness') of the left grid point
                // in outer circle.
            for (ix = idx1; ix <= idx2; ix++)
            {
                ttime[idx0] = ttime_outer[k];
                idx0++;
                k++;
            }
        }

        // Expand to top.
        iy = iy0 + radius;
        if (iy < ny)
        {
            idx1 = imax(0, ix0 - radius + 1);
                // Index of left X.
            idx2 = imin(nx-1, ix0 + radius - 1);
                // Index of right X.
            idx0 = nx * (iy - 1) + idx1;
                // Index (in array 'slowness') of the left grid point
                // in inner circle.
            k = 0;
            for (ix = idx1; ix <= idx2; ix++)
            {
                slowness_inner[k] = slowness[idx0];
                ttime_inner[k] = ttime[idx0];
                slowness_outer[k] = slowness[idx0 + nx];
                idx0++;
                k++;
            }

            travel_time_2d_expand(slowness_inner, slowness_outer,
                    ttime_inner, ttime_outer, sindex, k, h, &n1, &n2);

            *nsingular += n1;
            *nneg += n2;

            k = 0;
            idx0 = nx * iy + idx1;
                // Index (in array 'slowness') of the left grid point
                // in outer circle.
            for (ix = idx1; ix <= idx2; ix++)
            {
                ttime[idx0] = ttime_outer[k];
                idx0++;
                k++;
            }
        }


        // Four corner grid points.

        double t0, t1, t2, s;

        // Lower-left corner.
        if (ix0 - radius >= 0 && iy0 - radius >= 0)
        {
            idx0 = nx * (iy0 - radius) + (ix0 - radius);
            t0 = ttime[idx0 + 1 + nx];
            t1 = ttime[idx0 + nx];
            t2 = ttime[idx0 + 1];
            s = (slowness[idx0 + 1 + nx] + slowness[idx0 + nx] +
                slowness[idx0 + 1] + slowness[idx0]) / 4.;
            ftemp = 2. * (h*s) * (h*s) - (t2-t1) * (t2-t1);
            if (ftemp > 0)
                ttime[idx0] = t0 + sqrt(ftemp);
            else
            {
                ttime[idx0] = t0;
                (*nneg) ++;
            }
        }

        // Lower-right corner.
        if (ix0 + radius < nx && iy0 - radius >= 0)
        {
            idx0 = nx * (iy0 - radius) + (ix0 + radius);
            t0 = ttime[idx0 - 1 + nx];
            t1 = ttime[idx0 + nx];
            t2 = ttime[idx0 - 1];
            s = (slowness[idx0 - 1 + nx] + slowness[idx0 + nx] +
                slowness[idx0 - 1] + slowness[idx0]) / 4.;
            ftemp = 2. * (h*s) * (h*s) - (t2-t1) * (t2-t1);
            if (ftemp > 0)
                ttime[idx0] = t0 + sqrt(ftemp);
            else
            {
                ttime[idx0] = t0;
                (*nneg) ++;
            }
        }

        // Upper-left corner.
        if (ix0 - radius >= 0 && iy0 + radius < ny)
        {
            idx0 = nx * (iy0 + radius) + (ix0 - radius);
            t0 = ttime[idx0 + 1 - nx];
            t1 = ttime[idx0 - nx];
            t2 = ttime[idx0 + 1];
            s = (slowness[idx0 + 1 - nx] + slowness[idx0 - nx] +
                slowness[idx0 + 1] + slowness[idx0]) / 4.;
            ftemp = 2. * (h*s) * (h*s) - (t2-t1) * (t2-t1);
            if (ftemp > 0)
                ttime[idx0] = t0 + sqrt(ftemp);
            else
            {
                ttime[idx0] = t0;
                (*nneg) ++;
            }
        }

        // Upper-right corner.
        if (ix0 + radius < nx && iy0 + radius < ny)
        {
            idx0 = nx * (iy0 + radius) + (ix0 + radius);
            t0 = ttime[idx0 - 1 - nx];
            t1 = ttime[idx0 - nx];
            t2 = ttime[idx0 - 1];
            s = (slowness[idx0 - 1 - nx] + slowness[idx0 - nx] +
                slowness[idx0 - 1] + slowness[idx0]) / 4.;
            ftemp = 2. * (h*s) * (h*s) - (t2-t1) * (t2-t1);
            if (ftemp > 0)
                ttime[idx0] = t0 + sqrt(ftemp);
            else
            {
                ttime[idx0] = t0;
                (*nneg) ++;
            }
        }
    }


    free(ttime_inner);
    free(ttime_outer);
    free(slowness_inner);
    free(slowness_outer);
}


