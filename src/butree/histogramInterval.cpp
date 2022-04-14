#include "histogramInterval.h"
#include "../global/global.h"
#include <cassert>
#include <iostream>
using namespace std;

void histogramInterval::init_merge_info() {
    if (fanout == 0) {
        err_q = 0;
    } else {
        assert(lbd < ubd);
        long I_card = (ubd - lbd);
        double mu_q = 1.0 / I_card * prob_sum;
        err_q = cal_part_err_q(mu_q) + (I_card - fanout) * mu_q * mu_q;
    }
}


void histogramInterval::cal_merge_info(int h) {
    if(lbd >= ubd) {
        cout << "lbd = " << lbd << ", ubd = " << ubd << endl;
    }
    assert(lbd < ubd);
    assert(rSib != NULL);
//    assert(rSib->lbd == ubd);

    histogramInterval *h_rSib = static_cast<histogramInterval*>(rSib);

    if (h == 0) {
        if (fanout + rSib->fanout > fanThreashold) {
            merge_metric = 1e50;
            return;
        }
    }

    long I_card = (h_rSib->ubd - lbd);
    double mu_q = 1.0 / I_card * (prob_sum + rSib->prob_sum);
    double new_err_q = cal_part_err_q(mu_q) + h_rSib->cal_part_err_q(mu_q) + (I_card - fanout - h_rSib->fanout) * mu_q * mu_q;
    double delta_err_q = new_err_q - err_q - h_rSib->err_q;
//    merge_metric = delta_err_q;
    merge_metric = new_err_q;


//    double new_err_q = mu_q * (fanout + rSib->fanout);
//    delta_err_q = err_q + rSib->err_q - new_err_q;
}


bool histogramInterval::merge_with_rSib(int h, bool if_merge_lr) {
    assert (rSib != NULL);
    ubd = rSib->ubd;
//    err_q += (rSib->err_q + merge_metric);
    err_q = merge_metric;
//    err_q += (rSib->err_q - delta_err_q);


    if (rSib->fanout <= 0) {
        set_rSib(rSib->rSib);
        return true;
    }

    if (h == 0) {
        assert(fanout + rSib->fanout <= fanThreashold);
        if (fanout + rSib->fanout > fanThreashold) {
            merge_metric *= 2;
            return false;
        }
    }

    prob_sum += rSib->prob_sum;
    int new_fan = fanout + rSib->fanout;
    end_idx = rSib->end_idx;

    if (if_merge_lr) {
        lr->merge_and_self_update(rSib->lr, fanout, rSib->fanout, data[start_idx]);
        linear_loss = lr->cal_loss(data + start_idx, new_fan);
    }
    set_rSib(rSib->rSib);
    if (if_merge_lr && rSib) {
        rSib->lr->set_delta_x(data[rSib->start_idx] - data[start_idx]);
    }

    fanout = new_fan;
    return true;
}

double histogramInterval::cal_part_err_q(double _mu_q) {
    if (fanout <= 0) { return 0; }
    double s = 0;
    if (probs) {
        for (int i = 0; i < fanout; ++i) {
            double diff = probs[i + start_idx] - _mu_q;
            s += diff * diff;
        }
    } else {
        double diff = prob_sum / fanout - _mu_q;
        s = diff * diff * fanout;
    }
    return s;
}
