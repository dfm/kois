#ifndef _LIGHTCURVE_H_
#define _LIGHTCURVE_H_

void lightcurve (int n, double *t, int K, double texp, int np,
                 double *periods, double *epochs, double *durations,
                 double *rors, double *impacts, double mu1, double mu2,
                 double tol, int max_depth, double *f);

#endif
// _LIGHTCURVE_H_
