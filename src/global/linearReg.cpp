#include "linearReg.h"
#include <cmath>
#include <algorithm>
// y = a + bx

void linearReg(long *X, long y_start, double &a, double &b, int n) {
    if (n == 2) {
        b = 1.0 / (X[1] - X[0]);
        a = y_start - (b * X[0]);
        return;
    } else if (n == 3) {
        long x0 = X[0];
        long x1 = X[1];
        long x2 = X[2];
        long offset = MIN_LONG(x1 - x0, x2 - x1);
        b = 1.0 / offset;
        a = y_start + 1 - (b * X[1]);
        return;
    }
    double nu_b = 0;
    double de_b = 0;
    b = 0;
    double mean_x = 0;
    double mean_y = 0;
    for (int i = 0; i < n; ++i) {
        long y = y_start + i;
        de_b += X[i] * 1.0 * y;
        nu_b += X[i] * 1.0 * X[i];
        mean_x += X[i];
        mean_y += y;
    }

    mean_x /= n;
    mean_y /= n;
    de_b -= mean_x * mean_y * n;
    nu_b -= mean_x * mean_x * n;

    if (nu_b != 0) {
        b = de_b / nu_b;
    }  else {
        b = 0;
    }
    a = mean_y - b * mean_x;
}

void linearReg_stats(long *X, long y_start, const double a, const double b, int n, double &loss, double &max_diff, double &avg_exp_n_ops, double &avg_linear_n_ops) {
    loss = 0;
    max_diff = 0;
    avg_exp_n_ops = 0;
    avg_linear_n_ops = 0;

    for (int i = 0 ; i < n; ++i) {
        long y = y_start + i;

        double pred = a + b * X[i];
        pred = MIN_DOUBLE(MAX_DOUBLE(pred, 0), y_start + n - 1);
        double diff = abs(y - pred);
        if (diff > max_diff) {
            max_diff = diff;
        }
        loss += diff * diff;
    }
    double rmse =  sqrt(loss / n);
    avg_exp_n_ops = 1 + 2 * (rmse > 1 ? log(rmse) / log(2.0) : 0);

    avg_linear_n_ops = rmse / 4;
}


/*
void linearReg_stats(KV *key_values, long y_start, const double a, const double b, int n, double &loss, double &max_diff, double &avg_exp_n_ops, double &avg_linear_n_ops) {
    loss = 0;
    max_diff = 0;
    avg_exp_n_ops = 0;
    avg_linear_n_ops = 0;
    for (int i = 0 ; i < n; ++i) {
        long y = y_start + i;
        double x_i = key_values[i].x;

        double pred = a + b * x_i;
        pred = std::min<double>(std::max<double>(pred, 0), y_start + n - 1);
        double diff = abs(y - pred);
        if (diff > max_diff) {
            max_diff = diff;
        }
        loss += diff * diff;

        double _d = (int)(diff+1);
//        double _d = 1 + diff/4;
//        avg_exp_n_ops += 2 * log(_d) / log(2.0);
        double _n_op = 1 + 2 * (diff > 1 ? log(diff)/log(2.0) : 0);
        avg_exp_n_ops += _n_op;

        avg_linear_n_ops += diff/4;
    }

    avg_exp_n_ops = sqrt(loss / n);
    avg_exp_n_ops = 2 * (avg_exp_n_ops > 4 ? log(avg_exp_n_ops/4) / log(2.0) : 0);
//    avg_exp_n_ops /= n;

    avg_linear_n_ops = sqrt(loss / n) / 4;
//    avg_linear_n_ops /= n;
}
*/


/*
void cal_lr_loss_thread(KV *key_values, int N, const double a, const double b, double *losses, int thread_id, int start, int num) {
    double loss = 0;
    int max_pred = N - 1;
    for (int i = start; i < start + num; ++i) {
        double pred = a + b * key_values[i].x;
        pred = std::min<double>(std::max<double>(pred, 0), max_pred);
        double diff = abs(i - pred);
        loss += diff * diff;
    }
    losses[thread_id] = loss;
}

double parallel_cal_lr_loss(KV *key_values, int N, const double a, const double b) {
    double losses[nThreads];

    // make sure N >= nThreads
    int num = N / nThreads;
    int remainder = N % nThreads;

    vector<thread> ts;
    int start = 0;
    for (int i = 0; i < remainder; ++i) {
        ts.push_back(thread(&cal_lr_loss_thread, key_values, N, a, b, &(losses[0]), i, start, num + 1));
        start += num + 1;
    }

    for (int i = remainder; i < nThreads; ++i) {
        ts.push_back(thread(&cal_lr_loss_thread, key_values, N, a, b, &(losses[0]), i, start, num));
        start += num;
    }

    for_each(ts.begin(), ts.end(), mem_fn(&thread::join));

    double loss = 0;
    for (int i = 0; i < nThreads; ++i) {
        loss += losses[i];
    }
    return loss;
}
*/