#include "interval.h"
#include "buInterval.h"
#include <string>
using namespace std;

#ifndef DILI_INTERVALINSTANCE_H
#define DILI_INTERVALINSTANCE_H

struct intervalInstance {
    static interval* newInstance(int type) {
        interval *i_ptr = new buInterval();
        return i_ptr;
    }

    static interval* newInstance(const string &type) {
        interval *i_ptr = new buInterval();
        return i_ptr;
    }
};

#endif // DILI_INTERVALINSTANCE_H
