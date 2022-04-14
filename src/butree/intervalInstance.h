#include "interval.h"
#include "histogramInterval.h"
#include "buInterval.h"
#include <string>
using namespace std;

#ifndef DILI_INTERVALINSTANCE_H
#define DILI_INTERVALINSTANCE_H

struct intervalInstance {
    static interval* newInstance(int type) {
        if (type == 1) { // lr or bu
            interval *i_ptr = new buInterval();
            return i_ptr;
        }
        interval *i_ptr = new histogramInterval();
        return i_ptr;
    }

    static interval* newInstance(const string &type) {
        if (type == "lr" || type == "bu") { // lr
            interval *i_ptr = new buInterval();
            return i_ptr;
        }
        interval *i_ptr = new histogramInterval();
        return i_ptr;
    }
};

#endif // DILI_INTERVALINSTANCE_H
