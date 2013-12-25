#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "lightcurve.h"
#include "quad.h"

#define LC_FACTOR 0.166666666666666666666666666

double lc_evaluate_one (double t, double p, double t0, double tau, double ror,
                        double b, double mu1, double mu2)
{
    double opr = 1+ror, hp = 0.5*p, x;
    if (b >= opr) return 1.0;
    x = 2*(fmod(t+hp-t0, p)-hp)*sqrt(opr*opr - b*b)/tau;
    x = ldlc(ror, sqrt(b*b + x*x), mu1, mu2);
    return x;
}

double lc_evaluate (double t, int np, double *periods, double *epochs,
                    double *durations, double *rors, double *impacts,
                    double mu1, double mu2)
{
    int k;
    double f = 1.0;
    for (k = 0; k < np; ++k)
        f *= lc_evaluate_one(t, periods[k], epochs[k], durations[k], rors[k],
                             impacts[k], mu1, mu2);
    return f;
}

double lc_integrate (double t, double f0, double texp, int np,
                     double *periods, double *epochs, double *durations,
                     double *rors, double *impacts, double mu1, double mu2,
                     double tol, int max_depth, int depth)
{
    double tp = t+LC_FACTOR*texp,
           tm = t-LC_FACTOR*texp,
           fp = lc_evaluate (tp, np, periods, epochs, durations,
                             rors, impacts, mu1, mu2),
           fm = lc_evaluate (tm, np, periods, epochs, durations,
                             rors, impacts, mu1, mu2);
    double d = fabs(fp- 2*f0 + fm);
    if (d > tol && depth < max_depth) {
        fp = lc_integrate(tp, fp, texp/3.0, np, periods, epochs, durations,
                          rors, impacts, mu1, mu2, tol, depth+1, max_depth);
        f0 = lc_integrate(t, f0, texp/3.0, np, periods, epochs, durations,
                          rors, impacts, mu1, mu2, tol, depth+1, max_depth);
        fm = lc_integrate(tm, fm, texp/3.0, np, periods, epochs, durations,
                          rors, impacts, mu1, mu2, tol, depth+1, max_depth);
    }
    return (f0+fp+fm) / 3.0;
}

void lightcurve (int n, double *t, int K, double texp, int np,
                 double *periods, double *epochs, double *durations,
                 double *rors, double *impacts, double mu1, double mu2,
                 double tol, int max_depth, double *f)
{
    int i;
    double f0;

    for (i = 0; i < n; ++i) {
        f0 = lc_evaluate(t[i], np, periods, epochs, durations, rors, impacts,
                         mu1, mu2);
        f[i] = lc_integrate(t[i], f0, texp, np, periods, epochs, durations,
                            rors, impacts, mu1, mu2, tol, 0, max_depth);
    }
}
