#include "global.h"

#ifndef DILL_LINEARREG_H
#define DILL_LINEARREG_H

// y = a + bx
void linearReg(long *X, long y_start, double &a, double &b, int n);
void linearReg_stats(long *X, long y_start, const double a, const double b, int n, double &loss, double &max_diff, double &avg_exp_n_ops, double &avg_linear_n_ops);

#endif //DILL_LINEARREG_H
