#include "interval.h"

#ifndef DILI_LRINTERVAL_H
#define DILI_LRINTERVAL_H

struct buInterval: public interval {
    linearRegressor *merge_lr;

    buInterval(): interval(), merge_lr(NULL){}

    virtual void init_merge_info();
    virtual void init_merge_info_w_sampling();

    virtual void cal_merge_info(int h);
    virtual void cal_merge_info_w_sampling(int h);

    virtual bool merge_with_rSib(int h, bool if_merge_lr=false);



    virtual void free_data() {
        if (lr) { delete lr; lr = NULL; }
        if (merge_lr) { delete merge_lr; merge_lr = NULL; }
    }

    virtual ~buInterval() {
        if (merge_lr) { delete merge_lr; merge_lr = NULL; }
    }
};

#endif // DILI_LRINTERVAL_H
