#ifndef _LIGHTCURVE_H_
#define _LIGHTCURVE_H_

void lightcurve (double n, double *t, int K, double texp, int np,
                 double *periods, double *epochs, double *durations,
                 double *rors, double *impacts, int nbins, double *r,
                 double *ir, double *f);

#endif
// _LIGHTCURVE_H_
