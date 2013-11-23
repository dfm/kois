#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "lightcurve.h"

//
// Compute the area of a "star" occulted by a "planet".
//
// :param double r0:     The radius of the "star".
// :param double p:      The radius of the "planet".
// :param double b:      The impact parameter or center-to-center distance.
//
// :returns double area: The occulted area.
//
double lightcurve_occ_area (double r0, double p, double b)
{
    double r2, p2, b2, k1, k2, k3;

    if (b >= r0 + p) return 0.0;
    else if (b <= r0 - p) return M_PI * p * p;
    else if (b <= p - r0) return M_PI * r0 * r0;

    r2 = r0 * r0;
    p2 = p * p;
    b2 = b * b;

    k1 = acos(0.5 * (b2 + p2 - r2) / b / p);
    k2 = acos(0.5 * (b2 + r2 - p2) / b / r0);
    k3 = sqrt((p+r0-b) * (b+p-r0) * (b-p+r0) * (b+r0+p));

    return p2 * k1 + r2 * k2 - 0.5 * k3;
}

//
// Compute the limb darkened light curve for a transiting and eclipsing
// planet.
//
// :param double p:         The radius of the planet in stellar radii.
// :param int nbins:        The number of bins in the limb darkening profile.
// :param double r[nbins]:  The bin edges in the limb darkening profile
//                          measured in stellar radii.
// :param double ir[nbins]: The (relative) values of the limb darkening
//                          profile in the bins.
// :param int n:            The number of impact parameters to compute.
// :param double b[n]:      The values of impact parameter in stellar radii.
// :param int sgn[n]:       The relative position of the planet. If
//                          ``sgn[n] >= 0`` then the planet is transiting and
//                          if ``sgn[n] < 0`` the planet is being eclipsed.
// :param double lam[n]:    The output relative light curve.
//
void lightcurve_one (double p, int nbins, double *r, double *ir, int n,
                     double *b, int *sgn, double *lam)
{
    int i, j;
    double *areas = malloc(nbins * sizeof(double));

    // First, compute the normalization constant by integrating over the face
    // of the star.
    double norm = ir[0] * r[0] * r[0];
    for (i = 1; i < nbins; ++i)
        norm += ir[i] * (r[i] * r[i] - r[i - 1] * r[i - 1]);
    norm *= M_PI;

    // Compute the fraction of un-occulted flux for each time sample.
    for (i = 0; i < n; ++i) {
        if (sgn[i] >= 0) {
            // The planet is in front of the star.

            // Compute the array of occulted areas.
            for (j = 0; j < nbins; ++j)
                areas[j] = lightcurve_occ_area(r[j], p, b[i]);

            // Do the first order numerical integral over radial bins.
            lam[i] = areas[0] * ir[0];
            for (j = 1; j < nbins; ++j)
                lam[i] += ir[j] * (areas[j] - areas[j - 1]);
            lam[i] = 1.0 - lam[i] / norm;
        } else lam[i] = 1.0;
    }

    // Clean up.
    free(areas);
}

void lightcurve (double n, double *t, int K, double texp, int np,
                 double *periods, double *epochs, double *durations,
                 double *rors, double *impacts, int nbins, double *r,
                 double *ir, double *f)
{
    int i, j, k, ntot = n * K, *sgn = malloc(ntot * sizeof(int));
    double t1, v, b, b2, omb2, t0, P, hp, duration,
           *ttmp = malloc(ntot * sizeof(double)),
           *ftmp = malloc(ntot * sizeof(double)),
           *btmp = malloc(ntot * sizeof(double));

    for (i = 0; i < n; ++i) {
        f[i] = 1.0;
        if (K > 1)
            for (j = 0; j < K; ++j)
                ttmp[i*K+j] = t[i] - 0.5 * texp + j * texp / (K - 1);
        else
            ttmp[i] = t[i];
    }

    for (k = 0; k < np; ++k) {
        // Extract the planet parameters.
        b = impacts[k];
        b2 = b*b;
        omb2 = 1 - b2;
        t0 = epochs[k];
        P = periods[k];
        hp = 0.5*P;
        duration = durations[k];

        // Compute the impact parameter as a function of time.
        for (i = 0; i < ntot; ++i) {
            t1 = (fmod(ttmp[i] - t0 + hp, P) - hp) / duration;
            if (fabs(t1) < 1.0) {
                btmp[i] = sqrt(b2 + 4*omb2*t1*t1);
                sgn[i] = 1;
            } else sgn[i] = -1;
        }

        // Compute the light curve.
        lightcurve_one (rors[k], nbins, r, ir, ntot, btmp, sgn, ftmp);

        // Integrate over exposure time.
        for (i = 0; i < n; ++i) {
            v = 0.0;
            for (j = 0; j < K; ++j)
                v += ftmp[i*K+j];
            f[i] *= v / K;
        }
    }

    free(ftmp);
    free(ttmp);
    free(btmp);
    free(sgn);
}

// int main ()
// {
//     int i, n = 500, np = 1, nbins = 4, K = 3;
//     double f[500], t[500], texp = 0.02,
//            periods[] = {10.0},
//            epochs[] = {5.0},
//            durations[] = {0.5},
//            rors[] = {0.01},
//            impacts[] = {0.5},
//            r[] = {0.25, 0.5, 0.75, 1.0},
//            ir[] = {1.0, 0.75, 0.5, 0.1};
//
//     for (i = 0; i < n; ++i) t[i] = 0.021 * i;
//
//     lightcurve (n, t, K, texp, np, periods, epochs, durations, rors, impacts,
//                 nbins, r, ir, f);
//
//     for (i = 0; i < n; ++i)
//         printf("%e %e\n", t[i], f[i]);
//
//     return 0;
// }
