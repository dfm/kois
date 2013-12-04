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

void lightcurve (int n, double *t, int K, double texp, int np,
                 double *periods, double *epochs, double *durations,
                 double *rors, double *impacts, int nbins, double *r,
                 double *ir, double *f)
{
    int i, j, k, ntot = n * K, *sgn = malloc(ntot * sizeof(int));
    double t1, delta_t, v, b, b2, opr2, t0, P, hp, hpmt0, duration,
           *ttmp,
           *ftmp = malloc(ntot * sizeof(double)),
           *btmp = malloc(ntot * sizeof(double));

    if (K > 1) ttmp = malloc(ntot * sizeof(double));
    else ttmp = t;

    for (i = 0; i < n; ++i) {
        f[i] = 1.0;
        if (K > 1)
            for (j = 0; j < K; ++j)
                ttmp[i*K+j] = t[i] + ((j+0.5)/K - 0.5)*texp;
    }

    for (k = 0; k < np; ++k) {
        // Extract the planet parameters.
        opr2 = 1 + r[k];
        opr2 *= opr2;
        b = impacts[k];
        b2 = b*b;
        t0 = epochs[k];
        P = periods[k];
        hp = 0.5*P;
        hpmt0 = hp - t0;
        duration = durations[k];

        // Skip planets with no transits.
        if (opr2 > b2) {
            // // Compute the velocity based on the b=0 duration.
            // v = 2.0 / duration;
            //
            // // Compute the b=impact transit time.
            // delta_t = sqrt(opr2 - b2) / v;

            // FIXME.
            delta_t = 0.5 * duration;
            v = 2 * sqrt(opr2 - b2) / duration;

            // Compute the impact parameter as a function of time.
            for (i = 0; i < ntot; ++i) {
                // Find the time since transit.
                t1 = fmod(ttmp[i] + hpmt0, P) - hp;

                // If the planet is in transit, compute the impact parameter.
                if (fabs(t1) < delta_t) {
                    t1 *= v;
                    btmp[i] = sqrt(b2 + t1*t1);
                    sgn[i] = 1;
                } else sgn[i] = -1;
            }

            // Compute the light curve.
            lightcurve_one (rors[k], nbins, r, ir, ntot, btmp, sgn, ftmp);

            // Integrate over exposure time.
            if (K > 1) {
                for (i = 0; i < n; ++i) {
                    v = 0.0;
                    for (j = 0; j < K; ++j)
                        v += ftmp[i*K+j];
                    f[i] *= v / K;
                }
            } else {
                for (i = 0; i < n; ++i)
                    f[i] *= ftmp[i];
            }
        }
    }

    if (K > 1) free(ttmp);
    free(ftmp);
    free(btmp);
    free(sgn);
}
