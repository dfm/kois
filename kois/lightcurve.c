#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "lightcurve.h"
#include "quad.h"

void lightcurve (int n, double *t, int K, double texp, int np,
                 double *periods, double *epochs, double *durations,
                 double *rors, double *impacts, double mu1, double mu2,
                 double *f)
{
    int i, j, k, ntot = n * K;
    double t1, delta_t, v, b, b2, opr2, t0, P, hp, hpmt0, duration,
           *ttmp, *ftmp;

    if (K > 1) {
        ttmp = malloc(ntot * sizeof(double));
        ftmp = malloc(ntot * sizeof(double));
    } else {
        ttmp = t;
        ftmp = f;
    }

    for (i = 0; i < ntot; ++i) ftmp[i] = 1.0;
    if (K > 1) {
        for (i = 0; i < n; ++i)
            for (j = 0; j < K; ++j)
                ttmp[i*K+j] = t[i] + ((j+0.5)/K - 0.5)*texp;
    }

    for (k = 0; k < np; ++k) {
        // Extract the planet parameters.
        opr2 = 1 + rors[k];
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
                    ftmp[i] *= ldlc(rors[k], sqrt(b2 + t1*t1), mu1, mu2);
                }
            }
        }
    }

    // Integrate over exposure time.
    if (K > 1) {
        for (i = 0; i < n; ++i) {
            f[i] = 0.0;
            for (j = 0; j < K; ++j)
                f[i] += ftmp[i*K+j] / K;
        }
    }

    if (K > 1) {
        free(ttmp);
        free(ftmp);
    }
}
