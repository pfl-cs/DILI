#include "buInterval.h"
#include "../global/global.h"
#include <cassert>
#include <iostream>
using namespace std;


// use linear loss
void buInterval::init_merge_info() {
    if (fanout == 0) {
        linear_loss = 0;
    } else {
        assert(lbd < ubd);
        linear_loss = init_lr();
    }
}

void buInterval::init_merge_info_w_sampling() {
    if (fanout == 0) {
        linear_loss = 0;
    } else {
        assert(lbd < ubd);
        linear_loss = init_lr_w_sampling();
    }
}


void buInterval::cal_merge_info(int h) {
    if(lbd >= ubd) {
        cout << "lbd = " << lbd << ", ubd = " << ubd << endl;
    }
    assert(lbd < ubd);
    assert(rSib != NULL);

    if (!merge_lr) {
        merge_lr = new linearRegressor();
    }
    merge_lr->merge(lr, rSib->lr, fanout, rSib->fanout, data[start_idx]);

    if (h == 0) {
        if (fanout + rSib->fanout > fanThreashold) {
            merge_metric = 1e50;
            return;
        }
    }

    double new_linear_loss = merge_lr->cal_loss(data + start_idx, fanout + rSib->fanout);
    double delta_linear_loss = new_linear_loss - linear_loss - rSib->linear_loss;
    merge_metric = delta_linear_loss;
}

void buInterval::cal_merge_info_w_sampling(int h) {
    if(lbd >= ubd) {
        cout << "lbd = " << lbd << ", ubd = " << ubd << endl;
    }
    assert(lbd < ubd);
    assert(rSib != NULL);

    if (!merge_lr) {
        merge_lr = new linearRegressor();
    }
    merge_lr->merge_w_sampling(lr, rSib->lr, fanout, rSib->fanout, data[start_idx]);

    if (h == 0) {
        if (fanout + rSib->fanout > fanThreashold) {
            merge_metric = 1e50;
            return;
        }
    }

    double new_linear_loss = merge_lr->cal_loss_w_sampling(data + start_idx, fanout + rSib->fanout);
    double delta_linear_loss = new_linear_loss - linear_loss - rSib->linear_loss;
    merge_metric = delta_linear_loss;
}


bool buInterval::merge_with_rSib(int h, bool if_merge_lr) {
    assert (rSib != NULL);
    ubd = rSib->ubd;
    if (h == 0) {
        assert(fanout + rSib->fanout <= fanThreashold);
        if (fanout + rSib->fanout > fanThreashold) {
            merge_metric = 2e50;
            return false;
        }
    }

//    assert(rSib->fanout > 0);
    if (rSib->fanout > 0) {
        prob_sum += rSib->prob_sum;
        fanout += rSib->fanout;

        end_idx = rSib->end_idx;
        assert(end_idx == start_idx + fanout);
        linear_loss += (rSib->linear_loss + merge_metric);
//        cout << "lr->delta_x = " << lr->get_delta_x() <<  ", merge_lr->delta_x = "  <<  merge_lr->get_delta_x() << endl;
        lr->copy_from(merge_lr);

//        cout << "lr->delta_x = " << lr->get_delta_x() <<  endl;
//    lr = merge_lr;
//    merge_lr = NULL;
    }

    set_rSib(rSib->rSib);
    if (rSib) {
        assert(end_idx == rSib->start_idx);
        long dx = data[rSib->start_idx] - data[start_idx];
        buInterval *l_rSib = static_cast<buInterval*>(rSib);
        l_rSib->lr->set_delta_x(dx);
        if (l_rSib->merge_lr) {
            l_rSib->merge_lr->set_delta_x(dx);
        }
    }

    return true;
}
