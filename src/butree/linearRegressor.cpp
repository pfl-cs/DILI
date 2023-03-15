#include "linearRegressor.h"
#include <iostream>
#include <cstdio>
using namespace std;


const int sampling_factor = 2;

void linearRegressor::copy_from(linearRegressor *rhs) {
    a = rhs->a;
    b = rhs->b;
    _de_b = rhs->_de_b;
    _nu_b = rhs->_nu_b;
    sum_x = rhs->sum_x;
    delta_x = rhs->delta_x;
}

void linearRegressor::init(long *_data, int left_start_idx, int start_idx, int fanout) {

    a = 0;
    b = 0;
    _nu_b = 0;
    _de_b = 0;
    sum_x = 0;
//    sum_y = 0;
    long *data = _data + start_idx;
    long left_min = _data[left_start_idx];
    long x_min = data[0];

    delta_x = x_min - left_min;
    if (fanout <= 0) { return;}


    for (int i = 0; i < fanout; ++i) {
        long x = data[i] - x_min;
        _de_b += x * 1.0 * i;
        _nu_b += 1.0 * x * x;
        sum_x += x;
    }
    cal_ab(fanout, x_min);
}

void linearRegressor::init_w_sampling(long *_data, int left_start_idx, int start_idx, int fanout) {

    a = 0;
    b = 0;
    _nu_b = 0;
    _de_b = 0;
    sum_x = 0;
//    sum_y = 0;
    long *data = _data + start_idx;
    long left_min = _data[left_start_idx];
    long x_min = data[0];

    delta_x = x_min - left_min;
    if (fanout <= 0) { return;}


    int i = 0;
    for (int j = 0; j < fanout; j += sampling_factor) {
        long x = data[j] - x_min;

        _de_b += x * 1.0 * i;
        _nu_b += 1.0 * x * x;
        sum_x += x;
        ++i;
    }
    cal_ab_w_sampling(fanout, x_min);
}


void linearRegressor::cal_ab(int fanout, long x_min) {
    double mean_x = sum_x / fanout;
    double mean_y = (fanout - 1) / 2.0;
    double de_b = _de_b - mean_x * mean_y * fanout;
    double nu_b = _nu_b - mean_x * mean_x * fanout;

    if (nu_b != 0) {
        b = de_b / nu_b;
    }  else {
        b = 0;
    }
    a = mean_y - b * (mean_x + x_min);
}

void linearRegressor::cal_ab_w_sampling(int fanout, long x_min) {
    int n = fanout / sampling_factor;

    double mean_x = sum_x / n;
    double mean_y = (n - 1) / 2.0;
    double de_b = _de_b - mean_x * mean_y * fanout;
    double nu_b = _nu_b - mean_x * mean_x * fanout;

    if (nu_b != 0) {
        b = de_b / nu_b;
    }  else {
        b = 0;
    }
    a = mean_y - b * (mean_x + x_min);

    a *= 2;
    b *= 2;
}


void linearRegressor::merge_and_self_update(linearRegressor *rhs, int left_fan, int right_fan, long x_min) {
    if (left_fan <= 0) {
        if (right_fan > 0) {
            long dx = delta_x;
            copy_from(rhs);
            set_delta_x(dx);
        }
        return;
    }

    if (right_fan > 0) {
        long dx = rhs->delta_x;
        double _relative_sum_x = rhs->sum_x + dx * right_fan;
        double _relative_de_b = rhs->_de_b + dx * 0.5 * right_fan * (right_fan - 1);
        double _relative_nu_b = rhs->_nu_b + 2.0 * dx * rhs->sum_x + right_fan * dx * dx;

        _de_b += _relative_de_b + _relative_sum_x * left_fan;
        _nu_b += _relative_nu_b;
        sum_x += _relative_sum_x;

        cal_ab(left_fan + right_fan, x_min);
    }
}

/*
void linearRegressor::merge(linearRegressor *lhs, linearRegressor *rhs, int left_fan, int right_fan, long x_min) {
    if (left_fan <= 0) {
        if (right_fan > 0) {
            copy_from(rhs);
            set_delta_x(lhs->delta_x);
        }
        return;
    }
    if (right_fan <= 0) {
        copy_from(lhs);
    }

    if (right_fan > 0) {
        this->delta_x = lhs->delta_x;
        long dx = rhs->delta_x;
        double _relative_sum_x = rhs->sum_x + dx * right_fan;
//        double _relative_de_b = rhs->_de_b + dx * 0.5 * right_fan * (right_fan - 1);
        double _relative_de_b = rhs->_de_b + dx * 0.5 * right_fan * (right_fan - 1);
//        double _relative_nu_b = rhs->_nu_b + 2.0 * dx * rhs->sum_x + right_fan * dx * dx;

        _de_b = lhs->_de_b + _relative_de_b + _relative_sum_x * left_fan;
        _nu_b = lhs->_nu_b + rhs->_nu_b;
        sum_x = lhs->sum_x + _relative_sum_x;
        cal_ab(left_fan + right_fan, x_min);
    }
}
*/

void linearRegressor::merge(linearRegressor *lhs, linearRegressor *rhs, int left_fan, int right_fan,
                            long x_min) {
    if (left_fan <= 0) {
        if (right_fan > 0) {
            copy_from(rhs);
            set_delta_x(lhs->delta_x);
        }
        return;
    }
    if (right_fan <= 0) {
        copy_from(lhs);
    }

    if (right_fan > 0) {
        this->delta_x = lhs->delta_x;
        long dx = rhs->delta_x;
        double _relative_sum_x = rhs->sum_x + dx * right_fan;
        double _relative_de_b = rhs->_de_b + dx * 0.5 * right_fan * (right_fan - 1);
        double _relative_nu_b = rhs->_nu_b + 2.0 * dx * rhs->sum_x + right_fan * dx * dx;

        _de_b = lhs->_de_b + _relative_de_b + _relative_sum_x * left_fan;
        _nu_b = lhs->_nu_b + _relative_nu_b;
        sum_x = lhs->sum_x + _relative_sum_x;
        cal_ab(left_fan + right_fan, x_min);
    }
}



void linearRegressor::merge_w_sampling(linearRegressor *lhs, linearRegressor *rhs, int left_fan, int right_fan,
                                       long x_min) {
    if (left_fan <= 0) {
        if (right_fan > 0) {
            copy_from(rhs);
            set_delta_x(lhs->delta_x);
        }
        return;
    }

    if (right_fan <= 0) {
        copy_from(lhs);
    }

    int right_n = right_fan >> 1;
    if (right_n > 0) {
        this->delta_x = lhs->delta_x;
        long dx = rhs->delta_x;
        double _relative_sum_x = rhs->sum_x + dx * right_n;
        double _relative_de_b = rhs->_de_b + dx * 0.5 * right_n * (right_n - 1);
        double _relative_nu_b = rhs->_nu_b + 2.0 * dx * rhs->sum_x + right_n * dx * dx;

        _de_b = lhs->_de_b + _relative_de_b + _relative_sum_x * left_fan;
        _nu_b = lhs->_nu_b + _relative_nu_b;
        sum_x = lhs->sum_x + _relative_sum_x;
        cal_ab_w_sampling(left_fan + right_fan, x_min);
    }
}



double linearRegressor::cal_loss(long *data, int fanout) {
    double loss = 0;
    for (int i = 0; i < fanout; ++i) {
        double pred = a + b * data[i];
//        pred = std::max<double>(0, std::min<double>(pred, fanout - 1));
        double diff = pred - i;
        loss += diff * diff;
    }

    return loss;
}

double linearRegressor::cal_loss_w_sampling(long *data, int fanout) {
    double loss = 0;
    for (int i = 0; i < fanout; i += sampling_factor) {
        double pred = a + b * data[i];
//        pred = std::max<double>(0, std::min<double>(pred, fanout - 1));
        double diff = pred - i;
        loss += diff * diff;
    }

    loss *= 2;

    return loss;
}


void linearRegressor::print_cal_ab(int fanout) {
    double mean_x = sum_x / fanout;
    double mean_y = (fanout - 1) / 2.0;
    double de_b = _de_b - mean_x * mean_y * fanout;
    double nu_b = _nu_b - mean_x * mean_x * fanout;


    if (nu_b != 0) {
        b = de_b / nu_b;
    }  else {
        b = 0;
    }

    printf("de_b = %.2lf, n_b = %.2lf, mean_x = %.2lf, mean_y = %.2lf, b = %.2lf\n", de_b, nu_b, mean_x, mean_y, b);
    printf("mean_x * mean_x * fanout = %.2lf\n", mean_x * mean_x * fanout);
}

void linearRegressor::print(long *data, int fanout) {
    printf("a = %.2lf, b = %.2lf, sum_x = %.2lf, _de_b = %.2lf, _nu_b = %.2lf\n", a, b, sum_x, _de_b, _nu_b);
//    cout << "a = " << a << ", b = " << b << ", sum_x = " << sum_x << ", _de_b = " << _de_b << ", _nu_b = " << _nu_b << endl;
    if (data && fanout > 0) {
        print_cal_ab(fanout);
        for (int i = 0; i < fanout; ++i){
            long x = data[i];
            double pred = a + b * x;
            cout << "x = " << x << ", pred = " << pred << endl;
        }
    }
}