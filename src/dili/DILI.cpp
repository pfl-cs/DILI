#include "DILI.h"
#include "../global/linearReg.h"
#include "../global/global.h"
#include "../utils/data_utils.h"

#include <iostream>
#include <vector>
#include <cmath>
#include <queue>
#include <algorithm>
#include <functional>
#include <utility>
#include <cassert>
#include <thread>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <ctime>
#include <regex>
#include <cstdio>

#include <unistd.h>
#include <sys/ioctl.h>
#include <linux/perf_event.h>
#include <asm/unistd.h>

using namespace std;


//let i = returned value,  range_tos[i] <= key < range_tos[i+1]
namespace diliFunc {
    pair<diliNode **, double *>
    create_children(const int &height, diliNode **parents, int n_parents, double *parents_range_froms,
                    keyType *split_keys_for_children, recordPtr *ptrs, int n_keys,
                    int &act_total_N_children) {
        act_total_N_children = 0;
        for (int i = 0; i < n_parents; ++i) {
            act_total_N_children += parents[i]->get_fanout();
        }
        diliNode **_children = new diliNode *[act_total_N_children];

        double *children_range_froms = NULL;
        if (height > 1) {
            children_range_froms = new double[act_total_N_children + 1];
        }

        diliNode *child = NULL;
        int last_idx = 0;
        int cursor = 0;

//        cout << "****height = " << height << ", n_parents = " << n_parents << ", n_keys = " << n_keys << endl;

        for (int i = 0; i < n_parents; ++i) {
            diliNode *parent = parents[i];
            int fanout = parent->fanout;
            parent->pe_data = new pairEntry[fanout];
            double range_from = parents_range_froms[i];
            double parent_range_to = parents_range_froms[i + 1];
//        parent->children = new diliNode*[fanout];
//        double range_from = parent->range_from;

            for (int child_id = 0; child_id < fanout; ++child_id) {
                double range_to = (1.0 * (child_id + 1) - parent->a) / parent->b;
                if (range_to > parent_range_to) {
                    assert(child_id >= (fanout - 1));
                    range_to = parent_range_to;
                }
                if (child_id == fanout - 1) {
                    range_to = parent_range_to;
                }
                int idx = data_utils::array_lower_bound(split_keys_for_children, range_to, 0, n_keys);
                int n_keys_this_child = idx - last_idx;

                if (last_idx > idx) {
                    cout << "child_id = " << child_id << ", range_from = " << range_from << ", range_to = " << range_to
                         << ", parent_range_to = " << parent_range_to << endl;
                }
                assert(idx >= last_idx);

                if (height > 1) {
                    diliNode *int_node = new diliNode(true);
//                int_node->set_range(range_from, range_to);
                    int_node->cal_lr_params(split_keys_for_children + last_idx, n_keys_this_child);
//                int_node->children_init();
                    child = int_node;
                    children_range_froms[cursor] = range_from;
                    _children[cursor++] = child;
                } else {
                    child = new diliNode(false);
                    _children[cursor++] = child;
                }
//            parent->children[child_id] = child;
                parent->pe_data[child_id].setChild(child);

                last_idx = idx;
                range_from = range_to;
            }

//        std::copy(parent->children, parent->children + fanout, _children + cursor);
//        cursor += fanout;
        }
        if (height > 1) {
            children_range_froms[cursor] = parents_range_froms[n_parents];
        }

        assert(last_idx = n_keys);
        return make_pair(_children, children_range_froms);
    }
}

void DILI::save(const string &path) {
    FILE *fp = NULL;

    if (NULL == (fp = fopen(path.c_str(), "wb"))) {
        cout << path << " cannot be created." << endl;
        exit(1);
    }

//    saveNode(root, fp);
    root->save(fp);
    fclose(fp);
}

void DILI::load(const string &path) {
    FILE *fp = NULL;

    if (NULL == (fp = fopen(path.c_str(), "rb"))) {
        cout << path << " cannot be opened." << endl;
        exit(1);
    }
    root = new diliNode(true);
    root->load(fp);
    fclose(fp);
}

void DILI::bulk_load(const keyArray &keys, const recordPtrArray &ptrs, long n_keys) { //}, const string &mirror_dir, const string &layout_conf_path, int interval_type) {
    const int interval_type = 1;
    l_matrix mirror;
    build_ideal_mirror(keys, nullptr, n_keys, mirror, mirror_dir, interval_type);
//    build_mirror(keys, nullptr, n_keys, mirror, mirror_dir, interval_type);

//    cout << "----mirror.layout:------" << endl;
//    for (size_t i = 0; i < mirror.size(); ++i) {
//        cout << mirror[i].size() << " " << endl;
//    }
//    cout << endl;

    cout << "Building " << name() << "......" << endl;
    build_from_mirror(mirror, keys, ptrs, n_keys);
}

void DILI::bulk_load(const std::vector< pair<keyType, recordPtr> > &bulk_load_data) {
    size_t N = bulk_load_data.size();
    keyArray keys = std::make_unique<keyType []>(N + 1);
    recordPtrArray ptrs = std::make_unique<recordPtr []>(N + 1);
    for (size_t i = 0; i < N; ++i) {
        keys[i] = bulk_load_data[i].first;
        ptrs[i] = bulk_load_data[i].second;
    }
    keys[N] = keys[N-1] + 1;
    ptrs[N] = -1;
    bulk_load(keys, ptrs, static_cast<long>(N));
}
