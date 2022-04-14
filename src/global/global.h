#include <vector>
#include <iostream>

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

#define KEY_TYPE long
#define PAYLOAD_TYPE long

extern const int totalDataSize;
extern const long halfN;
extern const long query_step;
extern const long query_start_idx;

extern const long n_query_keys;
extern const double R1;
extern const double R2;
extern const double R3;

//extern const int minFan;
extern const int fanThreashold;

#define LEAF_MAX_CAPACIY 8192
#define minFan 2

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

typedef std::vector<int> intVec;
typedef std::vector<double> doubleVec;
typedef std::vector<long> longVec;
typedef std::vector<double> d_vector;
typedef std::vector< std::vector<double> > d_matrix;
typedef std::vector< std::vector<long> > l_matrix;
typedef std::vector< std::vector<int> > i_matrix;


struct diliNode;
struct fan2Leaf;

struct keyPayload {
    long key; // key >= 0: paylaod is the payload of key, key == -1: payload is a child; key < - 1: this position is empty
    union {long payload; diliNode *child; fan2Leaf *fan2child; };

    keyPayload(): key(-3), payload(-3) {}
    inline void setNull() { key = -3; payload = -3;}


//    inline long getPayload() { return payload;}
//    inline dillNode* getChild() { return child;}
//    inline fan2Leaf* getFan2Child() { return fan2child;}

    inline void assign(long _k, long _p) { key = _k; payload = _p; }
    inline void setChild(diliNode *_child) { key = -1; child = _child;};
    inline void setFan2Child(fan2Leaf *_fan2child) { key = -2; fan2child = _fan2child;};

};


#endif //DILI_GLOBAL_H
