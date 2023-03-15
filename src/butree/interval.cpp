#include "interval.h"
#include <cassert>
#include <iostream>
using namespace std;

const long *interval::data = NULL;
const double *interval::probs = NULL;

double interval::init_lr(bool force_init) {
    if (force_init) {
        if (lr) {
            delete lr;
            lr = NULL;
        }
    }
    if (!lr) {
        lr = new linearRegressor();
        if (lSib) {
            lr->init(data, lSib->start_idx, start_idx, fanout);
        } else {
            lr->init(data, start_idx, start_idx, fanout);
        }
        return lr->cal_loss(data + start_idx, fanout);
    } else {
        return linear_loss;
    }
}

double interval::init_lr_w_sampling(bool force_init) {
    if (force_init) {
        if (lr) {
            delete lr;
            lr = NULL;
        }
    }
    if (!lr) {
        lr = new linearRegressor();
        if (lSib) {
            lr->init_w_sampling(data, lSib->start_idx, start_idx, fanout);
        } else {
            lr->init_w_sampling(data, start_idx, start_idx, fanout);
        }
        return lr->cal_loss_w_sampling(data + start_idx, fanout);
    } else {
        return linear_loss;
    }
}

bool interval::check_lr_delta_x() {
    if (lSib) {
        if (lr && lSib->lr) {
//            if (fanout > 0) {
                bool status = data[lSib->start_idx] + lr->get_delta_x() == data[start_idx];
                if (!status) {
                    cout << "l_start_idx = " << lSib->start_idx << ", start_idx = " << start_idx << endl;
                    cout << "l_data[0] = " << data[lSib->start_idx] << ", data[0] = " << data[start_idx] << ", delta_x = " << lr->get_delta_x() << endl;
                }
                return status;
//            }
        }
    }
    return true;
}

