#include "interval.h"
#include "../global/global.h"
#include <vector>
#include <queue>
#include <set>
#include <string>

#ifndef DILI_INTERVAL_UTILS_H
#define DILI_INTERVAL_UTILS_H

class interval_comparator {
public:
    bool operator() (const interval *lhs, const interval *rhs) const {
        if (lhs->merge_metric == rhs->merge_metric) {
            return (lhs->lbd < rhs->lbd);
        }
        return (lhs->merge_metric < rhs->merge_metric);
//        return (lhs->delta_err_q < rhs->delta_err_q);
//        return (lhs->err_q < rhs->err_q);
    }
};

typedef std::pair<double, long> dl_pair;

class pair_max_comparator {
public:
    bool operator() (const dl_pair &lhs, const dl_pair &rhs) const {
        return (lhs.first < rhs.first);
    }
};

//typedef std::priority_queue<interval*, std::vector<interval*>, interval_comparator> interval_heap;
typedef std::set<interval*, interval_comparator> interval_set;
//typedef std::set<interval*> interval_set;

typedef std::priority_queue< dl_pair, std::vector<dl_pair>, pair_max_comparator > pair_max_heap;

// longVec &complete_borders,
void get_complete_partition_borders(const keyArray X, const doubleArray &probs, long N, int h, long min_border_size, longVec &complete_borders, doubleVec &complete_avg_rmses, int interval_type=0);
//void get_complete_partition_borders(long *X, double *probs, long N, int h, long min_border_size, longVec &complete_borders, int interval_type=0);

void build_mirror(const keyArray &X, const doubleArray &probs, long N, l_matrix &mirror, const std::string &mirror_dir, int interval_type);
void build_ideal_mirror(const keyArray &X, const doubleArray &probs, long N, l_matrix &mirror, const std::string &mirror_dir, int interval_type);
void build_mirror_from_given_layout(const keyArray &X, const doubleArray &probs, long N, l_matrix &mirror, const std::string &mirror_dir, const longVec &layout, int interval_type = 0);

void restore_mirror(const std::string &mirror_dir, l_matrix &mirror, bool ideal=false);
bool restore_complete_borders(const std::string &mirror_dir, const int h, longVec &borders);
bool restore_complete_borders_and_losses(const std::string &mirror_dir, const int h, longVec &borders, doubleVec &losses);

inline int find_nearest(long N, long n_nodes, long next_n_nodes) {
    return 1.0 * N * next_n_nodes / n_nodes;
}
#endif // DILI_INTERVAL_UTILS_H
