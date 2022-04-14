#include "global.h"

#ifndef DILL_LINEARREG_H
#define DILL_LINEARREG_H

//void linearReg_w_expanding(long *X, double &a, double &b, int n, int expanded_n, bool use_simple_strategy=false);
//void linearReg_at_least_four(long *X, double &a, double &b, int n);
//void linearReg_w_simple_strategy(long *X, double &a, double &b, int n);

// y = a + bx
void linearReg(long *X, long y_start, double &a, double &b, int n);
void linearReg_stats(long *X, long y_start, const double a, const double b, int n, double &loss, double &max_diff, double &avg_exp_n_ops, double &avg_linear_n_ops);

//void linearReg_stats(KV *key_values, long y_start, const double a, const double b, int n, double &loss, double &max_diff, double &avg_exp_n_ops, double &avg_linear_n_ops);

#endif //DILL_LINEARREG_H
