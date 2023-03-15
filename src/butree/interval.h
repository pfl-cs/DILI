#include "linearRegressor.h"
#include "../global/global.h"
#include <cstdio>
#include <cstdlib>
#include <cstring>
#ifndef DILI_INTERVAL_H
#define DILI_INTERVAL_H

struct interval{
    static const keyType *data;
    static const double *probs;

    bool valid;
    int fanout;
    long lbd;
    long ubd; // lbd < x <= ubd
//    double rj;
//    double tj;

    double merge_metric;

    long start_idx;
    long end_idx;

    interval *rSib;
    interval *lSib;

    linearRegressor *lr;
    double linear_loss;

    double prob_sum;


    //    interval(): lSib(NULL), rSib(NULL), lr(NULL), merge_lr(NULL), linear_loss(0), prob_sum(0), err_q(0),  valid(true) {}
    interval(): lSib(NULL), rSib(NULL), lr(NULL), linear_loss(0), prob_sum(0), valid(true) {}

    void setInvalid() { valid = false; }
    void set_rSib(interval *_rSib) {
        this->rSib = _rSib;
        if (_rSib) {
            _rSib->lSib = this;
        }
    }
    double init_lr(bool force_init=false);
    double init_lr_w_sampling(bool force_init=false);

    bool check_lr_delta_x();

    virtual void init_merge_info() = 0;
    virtual void init_merge_info_w_sampling() = 0;

    virtual void cal_merge_info(int h) = 0;
    virtual void cal_merge_info_w_sampling(int h) = 0;

    virtual bool merge_with_rSib(int h, bool if_merge_lr=false) = 0;


    virtual void free_data() {
        if (lr) { delete lr; lr = NULL; }
    }

    virtual ~interval() {
        if (lr) { delete lr; lr = NULL; }
    }

};

#endif // DILI_INTERVAL_H
