#include <vector>
#include <iostream>
#include "global_typedef.h"

#ifndef DILI_GLOBAL_H
#define DILI_GLOBAL_H

#define LR_PRED(a, b, key, fanout) std::min<int>(std::max<int>(static_cast<int>(a + b * key), 0), fanout - 1)
//inline int LR_PRED(double a, double b, const long &key, int fanout) { return std::min<int>(std::max<int>(a + b * key, 0), fanout - 1); }
#define MIN_INT(a, b) std::min<int>(a, b)
#define MAX_INT(a, b) std::max<int>(a, b)

#define MIN_LONG(a, b) std::min<long>(a, b)
#define MAX_LONG(a, b) std::max<long>(a, b)

#define MIN_DOUBLE(a, b) std::min<double>(a, b)
#define MAX_DOUBLE(a, b) std::max<double>(a, b)


extern long totalDataSize;
extern long halfN;
extern long query_step;
extern long query_start_idx;

extern const long n_query_keys;
extern const double R1;
extern const double R2;
extern const double R3;

//extern const int minFan;
extern const int fanThreashold;

#define LEAF_MAX_CAPACIY 8192
#define minFan 2

#define MIN_KEY(a, b) std::min<keyType>(a, b)

//extern const int Delta;
//extern const int nThreads;
//extern const int one_in_n;
//extern const double RATIO;
//extern const int maxFan;
//extern const int minFanforSplit;


extern double RHO;
extern int buMinFan;
extern double max_expanding_ratio;
extern double retrain_threshold;

extern long num_adjust_stats;


struct diliNode;
struct fan2Leaf;

struct SearchBound {
    size_t start;
    size_t stop;
};

struct pairEntry {
    keyType key; // key >= 0: ptr is the index of the record in the data array, key == -1: ptr is a child; key < - 1: this position is empty
    union {recordPtr ptr; diliNode *child; fan2Leaf *fan2child; };

    pairEntry(): key(-3), ptr(-3) {}
    inline void setNull() { key = -3; ptr = -3;}


//    inline long getPayload() { return payload;}
//    inline dillNode* getChild() { return child;}
//    inline fan2Leaf* getFan2Child() { return fan2child;}

    inline void assign(keyType _k, recordPtr _p) { key = _k; ptr = _p; }
    inline void setChild(diliNode *_child) { key = -1; child = _child;};
    inline void setFan2Child(fan2Leaf *_fan2child) { key = -2; fan2child = _fan2child;};

};


#endif //DILI_GLOBAL_H
