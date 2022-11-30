#include "diliNode.h"
#include "../butree/interval_utils.h"
#include "../global/global.h"
#include "../utils/data_utils.h"
#include <vector>
#include <cstdio>
#include <cstdlib>
#include <string>
#include <iostream>
#include <stack>
#ifndef DILI_DILI_H
#define DILI_DILI_H


namespace diliFunc {
    pair<diliNode **, double *>
    create_children(const int &height, diliNode **parents, int n_parents, double *parents_range_froms,
                    keyType *split_keys_for_children, recordPtr *ptrs, int n_keys,
                    int &act_total_N_children);
}

class DILI {
    diliNode *root;
    string mirror_dir;

public:
    DILI(): root(NULL) {
//        info();
//#ifndef ALLOW_FAN2_NODE
//        std::cout << "2.not define ALLOW_FAN2_NODE" << std::endl;
//#else
//        std::cout << "2.defined ALLOW_FAN2_NODE" << std::endl;
//#endif
        init_insert_aux_vars();
    }
    ~DILI() {
        clear();
    }

    void clear() {
        if (root) {
            delete root;
            root = NULL;
        }
        free_insert_aux_vars();
    }


    void init_insert_aux_vars() {
        dili_auxiliary::init_insert_aux_vars();
    }

    void free_insert_aux_vars() {
        dili_auxiliary::free_insert_aux_vars();
    }

    void set_mirror_dir(const std::string &dir) { mirror_dir = dir; }

    std::string name() {
#ifndef ALLOW_FAN2_NODE
        return "DILI";
#else
        return "simpleDILI";
#endif
    }

    long size() {
        std::stack<diliNode*> s;
        s.push(root);

        long size = 0;
        while (!s.empty()) {
            diliNode* node = s.top(); s.pop();
            if (!(node->is_internal())) {
                continue;
            }

            size += sizeof(diliNode) - 3 * sizeof(int);

            for (int i = 0; i < node->fanout; ++i) {
                pairEntry &kp = node->pe_data[i];
                if (kp.key < 0) {
                    size += sizeof(diliNode *);
                    if (kp.key == -1) {
                        s.push(kp.child);
                    }
                }
            }
        }
        return size;
    }

    long total_size() {
        std::stack<diliNode*> s;
        s.push(root);

        long size = 0;
        while (!s.empty()) {
            diliNode* node = s.top(); s.pop();
//            if (!(node->is_internal())) {
//                continue;
//            }

            size += sizeof(diliNode) - 3 * sizeof(int);
            size += sizeof(pairEntry) * node->fanout;

            for (int i = 0; i < node->fanout; ++i) {
                pairEntry &kp = node->pe_data[i];
                if (kp.key < 0) {
                    if (kp.key == -1) {
                        s.push(kp.child);
                    } else if (kp.key == -2) {
                        size += 4 * sizeof(long);
                    } else {
                        size += 2 * sizeof(long);
                    }
                } else {
                    size += 2 * sizeof(long);
                }
            }
        }
        return size;
    }

    void order_check() {
        keyType *tmp_keys = new keyType [totalDataSize];
        root->collect_all_keys(tmp_keys);
        cout << "---------check, N = " << root->num_nonempty << "----------" << endl;
        data_utils::check(tmp_keys, root->num_nonempty);
        cout << "----------------------------------------" << endl;
        delete[] tmp_keys;
    }

    void build_from_mirror(l_matrix &mirror, const keyArray &all_keys, const recordPtrArray &all_ptrs, long N) {
        size_t H = mirror.size();

//        cout << "+++H = " << H << endl;
        intVec n_nodes_each_level;
        intVec n_nodes_each_level_mirror;
        for (longVec &lv : mirror) {
            n_nodes_each_level_mirror.push_back(lv.size());
        }

        keyType **split_keys_list = new keyType *[H];
        split_keys_list[0] = all_keys.get();

        for (int height = H - 1; height > 0; --height) {
            longVec &lv = mirror[height-1];
            int n_split_keys = lv.size();

            long *split_keys = new long[n_split_keys + 1];
            std::copy(lv.begin(), lv.end(), split_keys);
            split_keys[n_split_keys] = all_keys[N-1] + 1;

            split_keys_list[height] = split_keys;
            data_utils::check(split_keys, n_split_keys+1);
        }
        root = new diliNode(true);
//    root->set_range(0, all_keys[N-1] + 1);
        int n_keys = n_nodes_each_level_mirror[H - 2];
        root->fanout = n_keys + 1;
        long ubd = split_keys_list[H-1][n_keys-1];
        long lbd = split_keys_list[H-1][0];
        root->b = 1.0 * n_keys / (ubd - lbd);
        root->a = -(root->b * lbd);
        n_nodes_each_level.push_back(1);
//    root->children_init();
        root->pe_data = new pairEntry[root->fanout];

//        cout << "lbd = " << lbd << ", ubd = " << ubd << ", n_keys = " << n_keys << ", max_key = " << all_keys[N-1] << endl;

        keyType lastone = split_keys_list[H - 2][n_nodes_each_level_mirror[H-3]-1];
//        cout << "a = " << root->a << ", b = " << root->b
//             << ", lr(x[0]) = " << int(root->a + root->b * split_keys_list[H - 1][0])
//             << ", lr(x[1]) = " << int(root->a + root->b * split_keys_list[H - 1][1])
//             << ", lr(ubd) = " << int(root->a + root->b * ubd)
//             << ", lr(lastone) = " << int(root->a + root->b * lastone) << endl;


//    root->cal_lr_params(split_keys_list[H-2], n_keys, child_fans_list[H-2]);
//    assert(root->range_to > all_keys[N-1]);

        diliNode **parents = new diliNode*[1];
        parents[0] = root;
        int n_parents = 1;
        double *parents_range_froms = new double[2];
        parents_range_froms[0] = 0;
        parents_range_froms[1] = all_keys[N-1] + 1;

        // height: the height of parents
        diliNode **children = NULL;
        int act_total_N_children = 0;

        for (int height = H - 1; height > 0; --height) {
//            cout << "--------height = " << height << "-----------" << endl;
            n_keys = N;

            if (height > 1) {
                n_keys = n_nodes_each_level_mirror[height-2];
            }
//            cout << "1. here is OK." << endl;

            recordPtr *ptrs = NULL;
            if (height == 1) {
                ptrs = all_ptrs.get();
            }

//            cout << "2. here is OK." << endl;

            pair<diliNode **, double *> _pair = diliFunc::create_children(height, parents, n_parents, parents_range_froms, split_keys_list[height-1],
                                                                          ptrs, n_keys, act_total_N_children);
            children = _pair.first;
            double *children_range_froms = _pair.second;
//            cout << "3. here is OK." << endl;
            delete[] parents;

//            cout << "4. here is OK." << endl;
            parents = children;

            delete[] parents_range_froms;
            parents_range_froms = children_range_froms;

            n_parents = act_total_N_children;
            n_nodes_each_level.push_back(act_total_N_children);
        }

//        cout << "-----here is OK." << endl;

        for (int height = 1; height < H; ++height) {
            delete[] split_keys_list[height];
        }
        delete[] split_keys_list;

        for (long i = 0; i < N; ++i) {
            diliNode *leaf = find_leaf(all_keys[i]);
//        leaf->tmp.push_back(all_keys[i]);
            leaf->inc_num_nonempty();
        }
//    long *tmp_keys = new long[totalDataSize];
//    long _j = 0;
//    for (auto i = 0; i < act_total_N_children; ++i) {
//        diliNode *leaf = children[i];
//        if (leaf->tmp.size() > 0) {
//            std::copy((leaf->tmp).begin(), (leaf->tmp).end(), tmp_keys + _j);
//            _j += leaf->tmp.size();
//        }
//    }
//    cout << "---------check, N = " << _j << "----------" << endl;
//    data_utils::check(tmp_keys, _j);
//    cout << "----------------------------------------" << endl;
//    delete[] tmp_keys;

//        cout << "****here is OK." << endl;

        long start_idx = 0;
        bool print = false;
        for (int i = 0; i < act_total_N_children; ++i) {
            diliNode *leaf = children[i];
            int _num_nonempty = leaf->num_nonempty;
            leaf->bulk_loading(all_keys.get() + start_idx, all_ptrs.get() + start_idx, print);
            start_idx += _num_nonempty;
        }
        if (start_idx != N) {
            cout << "error, start_idx = " << start_idx << ", N = " << N << endl;
        }
        assert(start_idx == N);

//        validness_check(all_keys, all_ptrs, N);

        root->trim();
//        cout << "+++1. here is OK." << endl;
        root->cal_num_nonempty();
//        cout << "+++2. here is OK." << endl;
#ifndef ALLOW_FAN2_NODE
        root->simplify();
//        cout << "+++3. here is OK." << endl;
#endif
        root->cal_avg_n_travs();
//        cout << "+++4. here is OK." << endl;
        root->init_after_bulk_load();
//        cout << "+++5. here is OK." << endl;
//        validness_check(all_keys, all_ptrs, N);

//        cout << "layout: [" << n_nodes_each_level[0];
//        for (size_t i = 1; i < n_nodes_each_level.size(); ++i) {
//            cout << ", " << n_nodes_each_level[i];
//        }
//        cout << "]" << endl;
    }

    void validness_check(keyType *keys, recordPtr *ptrs, int n_keys) {
        for (int i = 0; i < n_keys; ++i) {
            recordPtr pred = search(keys[i]);
            if (pred != ptrs[i]) {
                cout << "i = " << i << ", key = " << keys[i] << ", pred = " << pred << ", ptr = " << ptrs[i] << endl;
            }
            assert(pred == ptrs[i]);
        }
    }

    inline bool insert(const keyType &key, const recordPtr &ptr) { return root->insert(key, ptr); };
    inline bool insert(const pair<keyType, recordPtr> &p) { return root->insert(p.first, p.second); };
    inline bool erase(const keyType &key) { return 0 <= (root->erase(key)); }
    inline recordPtr delete_key(const keyType &key) {
        recordPtr ptr = static_cast<recordPtr>(-1);
        root->erase_and_get_ptr(key, ptr);
        return ptr;
    }

    void save(const string &path);
    void load(const string &path);
    diliNode* loadNode(FILE *fp);


    // only called on bulk loading stage
    inline diliNode* find_leaf(const keyType &key) {
        diliNode *node = root->find_child(key);
        while (node->is_internal()) {
            node = node->find_child(key);
        }
        return static_cast<diliNode*>(node);
    }


    inline long search(const keyType &key) {
//        std::cout << "******key = " << key << std::endl;
        diliNode *node = root;
        while (true) {
            int pred = LR_PRED(node->a, node->b, key, node->fanout);
            pairEntry &kp = node->pe_data[pred];
            if (kp.key == key) {
                return kp.ptr;
            } else if (kp.key == -1) {
                node = kp.child;
            } else if (kp.key == -2) {
                fan2Leaf *child = kp.fan2child;
                if (child->k1 == key) {
                    return child->p1;
                }
                if (child->k2 == key) {
                    return child->p2;
                }
                return -1;
            }
            else {
                return -1;
            }
        }
    }

    inline int range_query(const keyType &k1, const keyType &k2, recordPtr *ptrs) { return root->range_query(k1, k2, ptrs); }

    inline long search_w_print(const keyType &key) {
        return root->leaf_find_w_print(key);
    }

    inline int depth(const keyType &key) {
        diliNode *node = root;
        int depth = 0;
        while (true) {
            ++depth;
            int pred = LR_PRED(node->a, node->b, key, node->fanout);
            pairEntry &kp = node->pe_data[pred];
//            if (key == 572167479012l || key == 572167481113l) {
//                cout << "key = " << key << ", is_internal = " << node->is_internal() << ", pred = " << pred << ", fanout = "
//                     << node->fanout << ", num_nonempty = " << node->num_nonempty << ", kp.key = " << kp.key << endl;
//            }
            if (kp.key == key) {
                return depth;
            } else if (kp.key == -1) {
                node = kp.child;
            }
#ifndef ALLOW_FAN2_NODE
            else if (kp.key == -2) {
                return depth + 1;
            }
#endif
            else {
                return depth;
            }
        }
    }

    void bulk_load(const keyArray &keys, const recordPtrArray &ptrs, long n_keys);//, const string &mirror_dir, const string &layout_conf_path, int interval_type=1);
    void bulk_load(const std::vector< pair<keyType, recordPtr> > &bulk_load_data);
    void debug_help(bool print=false){} // will be deprecated
    void stats();
    void check_num_nonempty() { root->check_num_nonempty(); }
    double avg_depth(keyType *keys, long n_keys) {
        double avg_depth = 0;
        for (long i = 0; i < n_keys; ++i) {
            avg_depth += depth(keys[i]);
        }
        return avg_depth / n_keys;
    }
};

/*
inline int array_lower_bound(long *data, double x, int l, int r) {
    while (l < r) {
        int mid = (l + r) >> 1;
//        int mid = l + ((r - l) >> 1);

        if (data[mid] >= x) {
            r = mid;
        } else {
            l = mid + 1;
        }
    }
    return l;
}
 */


#endif //DILI_DILI_H
