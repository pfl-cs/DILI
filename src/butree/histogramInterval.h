#include "interval.h"

#ifndef DILI_HISTOGRAMINTERVAL_H
#define DILI_HISTOGRAMINTERVAL_H

struct histogramInterval: public interval{
    double err_q;
//    double new_err_q;
//    double delta_err_q;
//    double mu_q;

    // data(NULL), probs(NULL) delta_err_q(1e20),
    histogramInterval(): interval(), err_q(0) {}

    virtual void init_merge_info();
    virtual void cal_merge_info(int h);
    virtual bool merge_with_rSib(int h, bool if_merge_lr=false);
private:
    double cal_part_err_q(double _mu_q);

public:
    virtual void free_data() {
        if (lr) { delete lr; lr = NULL; }
    }
};


#endif //DILI_HISTOGRAMINTERVAL_H
