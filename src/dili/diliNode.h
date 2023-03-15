#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <queue>
#include <algorithm>
#include <cmath>
#include <cassert>
#include <stack>
#include "../global/global.h"
#include "../global/fan2Leaf.h"
//#include "../global/linearReg.h"

#ifndef DILI_DILINODE_H
#define DILI_DILINODE_H


using namespace std;
struct diliNode;
struct fan2Leaf;

namespace dili_auxiliary {
    extern keyType *retrain_keys;
    extern recordPtr *retrain_ptrs;
    void init_insert_aux_vars();
    void free_insert_aux_vars();
}


inline void linearReg_w_simple_strategy(const keyType *X, double &a, double &b, int n) {
    int left_n = n / 2;
    int right_n = n - left_n - 1;
    keyType x_middle = X[left_n];
    double left_slope = 1.0 * left_n / (x_middle - X[0]);
    double right_slope = 1.0 * right_n / (X[n-1] - x_middle);
    b = MAX_DOUBLE(left_slope, right_slope);
    a = -b * x_middle + left_n;
}

inline void linearReg_at_least_four(const keyType *X, double &a, double &b, int n) {
    double nu_b = 0;
    double de_b = 0;
    double mean_x = 0;
    for (int i = 0; i < n; ++i) {
        de_b += X[i] * 1.0 * i;
        nu_b += X[i] * 1.0 * X[i];
        mean_x += X[i];
    }

    mean_x /= n;
    double mean_y = (n - 1) / 2.0;
    de_b -= mean_x * mean_y * n;
    nu_b -= mean_x * mean_x * n;

    if (nu_b != 0) {
        b = de_b / nu_b;
    }  else {
        b = 0;
    }
    a = mean_y - b * mean_x;
}

inline void linearReg_with_max_b(const keyType *X, double &a, double &b, int n) {
    double nu_b = 0;
    double de_b = 0;
    double mean_x = 0;
    for (int i = 0; i < n; ++i) {
        de_b += X[i] * 1.0 * i;
        nu_b += X[i] * 1.0 * X[i];
        mean_x += X[i];
    }

    mean_x /= n;
    double mean_y = (n - 1) / 2.0;
    de_b -= mean_x * mean_y * n;
    nu_b -= mean_x * mean_x * n;

    if (nu_b != 0) {
        b = de_b / nu_b;
    }  else {
        b = 0;
    }
    a = mean_y - b * mean_x;
}


inline void linearReg_w_expanding(const keyType *X, double &a, double &b, int n, int expanded_n, bool use_simple_strategy) {
    if (!use_simple_strategy) {
//        linearReg_at_least_four(X, a, b, n);
        linearReg_with_max_b(X, a, b, n);
    } else {
        linearReg_w_simple_strategy(X, a, b, n);
    }
    if (expanded_n > n) {
        double expanding_ratio = 1.0 * expanded_n / (n + 1);
        b *= expanding_ratio;
        a = a * expanding_ratio + expanding_ratio;
    }
}

struct diliNode{
    int fanout;
    int meta_info;
    double a;
    double b;
    int num_nonempty;

    pairEntry *pe_data;


    double avg_n_travs_since_last_dist;
    long total_n_travs;
    long last_total_n_travs;
    int last_nn;
    int n_adjust;



    inline bool is_internal() const { return (meta_info & 1); }
    inline void set_leaf_flag() { meta_info &= ~1U;}
    inline void set_int_flag() { meta_info |= 1; }

    inline void set_fanout(int fan) { fanout = fan;}
    inline int get_fanout() { return fanout; }

    inline void set_n_adjust(int n) { meta_info = n << 1; }

    inline int get_n_adjust() { return n_adjust; }
    inline void inc_n_adjust() { ++n_adjust; }


    inline void set_num_nonempty(int n) { num_nonempty = n; }

    inline void init() {
        fanout = std::max<int>(num_nonempty, minFan);
        fanout <<= 1;
//        fanout *= 1.5;

//        fanout += (fanout * n_adjust) / 10;
//        fanout = std::max<int>(num_nonempty, minFan) * (1 + 0.1 * get_n_adjust());
        pe_data = new pairEntry[fanout];
    }



    diliNode(bool _is_internal): a(0), b(0), meta_info(_is_internal), fanout(0), pe_data(NULL), n_adjust(30), num_nonempty(0),
                                 total_n_travs(0), last_total_n_travs(0), last_nn(0), avg_n_travs_since_last_dist(1e10)  {}

    inline void init(const int &_num_nonempty) {
        num_nonempty = _num_nonempty;
        fanout = std::max<int>(_num_nonempty, minFan);
        fanout <<= 1;
        pe_data = new pairEntry[fanout];
    }

    inline void inc_num_nonempty() { ++num_nonempty; }
    inline int get_num_nonempty() { return num_nonempty; }

    int cal_num_nonempty() {
        if (num_nonempty <= 0) {
            assert(num_nonempty == 0);
            for (int i = 0; i < fanout; ++i) {
                pairEntry &pe = pe_data[i];
                if (pe.key >= 0) {
                    ++num_nonempty;
                } else if (pe.key == -1) {
                    num_nonempty += pe.child->cal_num_nonempty();
                } else if (pe.key == -2) {
                    num_nonempty += 2;
                }
            }
        }
        return num_nonempty;
    }

    void init_after_bulk_load() {
        last_total_n_travs = total_n_travs;
        last_nn = num_nonempty;
        for (int i = 0; i < fanout; ++i) {
            pairEntry &pe = pe_data[i];
            if (pe.key == -1) {
                pe.child->init_after_bulk_load();
            }
        }
    }

    void check_num_nonempty() {
        assert(num_nonempty > 1);
        if (!is_internal()) {
            assert(num_nonempty < LEAF_MAX_CAPACIY);
        }
        for (int i = 0; i < fanout; ++i) {
            pairEntry &pe = pe_data[i];
            if (pe.key == -1) {
                pe.child->check_num_nonempty();
            }
        }
    }


    inline recordPtr leaf_find(const keyType &key) const {
        int pred = LR_PRED(a, b, key, fanout);
//        cout << "pred = " << pred << endl;
        pairEntry &pe = pe_data[pred];
//        cout << "pe.key = " << pe.key << endl;
        if (pe.key == key) {
            return pe.ptr;
        }
        if (pe.key == -1) {
            return pe.child->leaf_find(key);
        }
        if (pe.key == -2) {
            fan2Leaf *child = pe.fan2child;
            if (child->k1 == key) {
                return child->p1;
            }
            if (child->k2 == key) {
                return child->p2;
            }
            return -1;
        }
        return -1;
    }

    inline int range_query_from(const keyType &k1, recordPtr *results) const {
        int j = 0;
        int pred = LR_PRED(a, b, k1, fanout);
        pairEntry &first_pe = pe_data[pred];
        if (first_pe.key == -1) {
            j = first_pe.child->range_query_from(k1, results);
        } else if (first_pe.key == -2) {
            fan2Leaf *fan2child = first_pe.fan2child;
            keyType _k1 = fan2child->k1;
            keyType _k2 = fan2child->k2;
            if (_k1 >= k1) {
                results[j++] = fan2child->p1;
                results[j++] = fan2child->p2;
            } else if (_k2 >= k1) {
                results[j++] = fan2child->p2;
            }
        } else if (first_pe.key >= k1) {
            results[j++] = first_pe.ptr;
        }

        for(int i = pred + 1; i < fanout; ++i) {
            pairEntry &pe = pe_data[i];
            if (pe.key >= 0) {
                results[j++] = pe.ptr;
            } else if (pe.key == -1) {
                pe.child->collect_all_ptrs(results+j);
                j += pe.child->num_nonempty;
            } else if (pe.key == -2) {
                fan2Leaf *fan2child = pe.fan2child;
                results[j++] = fan2child->p1;
                results[j++] = fan2child->p2;
            }
        }

        return j;
    }

    inline int range_query_to(const keyType &k2, recordPtr *results) const {
        int j = 0;
        int pred = LR_PRED(a, b, k2, fanout);
        for(int i = 0; i < pred; ++i) {
            pairEntry &pe = pe_data[i];
            if (pe.key >= 0) {
                results[j++] = pe.ptr;
            } else if (pe.key == -1) {
                pe.child->collect_all_ptrs(results+j);
                j += pe.child->num_nonempty;
            } else if (pe.key == -2) {
                fan2Leaf *fan2child = pe.fan2child;
                results[j++] = fan2child->p1;
                results[j++] = fan2child->p2;
            }
        }

        pairEntry &pe = pe_data[pred];
        if (pe.key == -1) {
            j += pe.child->range_query_to(k2, results+j);
        } else if (pe.key == -2) {
            fan2Leaf *fan2child = pe.fan2child;
            keyType _k1 = fan2child->k1;
            keyType _k2 = fan2child->k2;
            if (_k2 < k2) {
                results[j++] = fan2child->p1;
                results[j++] = fan2child->p2;
            } else if (_k1 < k2) {
                results[j++] = fan2child->p1;
            }
        } else if (pe.key < k2) {
            results[j++] = pe.ptr;
        }
        return j;
    }


    inline void collect_all_ptrs(recordPtr *results) const{
        int j = 0;
        for(int i = 0; i < fanout; ++i) {
            pairEntry &pe = pe_data[i];
            if (pe.key >= 0) {
                results[j++] = pe.ptr;
            } else if (pe.key == -1) {
                pe.child->collect_all_ptrs(results+j);
                j += pe.child->num_nonempty;
            } else if (pe.key == -2) {
                fan2Leaf *fan2child = pe.fan2child;
                results[j++] = fan2child->p1;
                results[j++] = fan2child->p2;
            }
        }
    }


    inline int range_query(const keyType &k1, const keyType &k2, recordPtr *results) const {
        int pred1 = LR_PRED(a, b, k1, fanout);
        int pred2 = LR_PRED(a, b, k2, fanout);

        if (pred1 == pred2) {
            pairEntry &pe = pe_data[pred1];
            if (pe.key == -1) {
                return pe.child->range_query(k1, k2, results);
            } else if (pe.key == -2) {
                int n = 0;
                fan2Leaf *leaf = pe.fan2child;
                keyType _k1 = leaf->k1;
                keyType _k2 = leaf->k2;
                if (_k1 >= k1 && _k1 < k2) {
                    results[n++] = leaf->p1;
                }
                if (_k2 >= k1 && _k2 < k2) {
                    results[n++] = leaf->p2;
                }
                return n;
            } else if (pe.key >= k1 && pe.key < k2) {
                results[0] = pe.ptr;
                return 1;
            }
        } else { // pred1 < pred2

            pairEntry &first_pe = pe_data[pred1];
            int n = 0;
            if (first_pe.key == -1) {
                n = first_pe.child->range_query_from(k1, results);
            } else if (first_pe.key == -2) {
                fan2Leaf *leaf = first_pe.fan2child;
                keyType _k1 = leaf->k1;
                keyType _k2 = leaf->k2;
                if (_k1 >= k1) {
                    results[0] = leaf->p1;
                    results[1] = leaf->p2;
                    n = 2;
                } else if (_k2 >= k1) {
                    results[1] = leaf->p2;
                    n = 1;
                }
            } else if (first_pe.key >= k1) {
                results[0] = first_pe.ptr;
                n = 1;
            }

            for (int i = pred1 + 1; i < pred2; ++i) {
                pairEntry &pe = pe_data[i];
                if (pe.key == -1) {
                    pe.child->collect_all_ptrs(results+n);
                    n += pe.child->num_nonempty;
                } else if (pe.key == -2) {
                    fan2Leaf *leaf = pe.fan2child;
                    results[n++] = leaf->p1;
                    results[n++] = leaf->p2;
                } else if (pe.key >= k1 && pe.key < k2) {
                    results[n++] = pe.ptr;
                }
            }

            pairEntry &final_pe = pe_data[pred2];
            if (final_pe.key == -1) {
                n += final_pe.child->range_query_to(k2, results+n);
            } else if (final_pe.key == -2) {
                fan2Leaf *leaf = final_pe.fan2child;
                keyType _k1 = leaf->k1;
                keyType _k2 = leaf->k2;
                if (_k2 < k2) {
                    results[n++] = leaf->p1;
                    results[n++] = leaf->p2;
                } else if (_k1 < k2) {
                    results[n++] = leaf->p1;
                }

            } else if (final_pe.key < k2) {
                results[n++] = final_pe.ptr;
            }
            return n;
        }
    }


    inline diliNode* find_child(const keyType &key) {
        int i = LR_PRED(a, b, key, fanout);
        return pe_data[i].child;
    }


    inline bool if_retrain() { return (!is_internal()) && (total_n_travs * last_nn >= ((last_total_n_travs * num_nonempty) << 1) ); }

    inline void put_three_keys(const keyType &k0, const recordPtr &p0, const keyType &k1, const recordPtr &p1, const keyType &k2, const recordPtr &p2) {
        double offset = fanout / 3.0;
        keyType s = MIN_KEY(k1 - k0, k2 - k1);
        b = offset / s;
        a = 0.5 + offset - b * k1;
        int pos0 = LR_PRED(a, b, k0, fanout);
        int pos1 = LR_PRED(a, b, k1, fanout);
        int pos2 = LR_PRED(a, b, k2, fanout);
        pe_data[pos0].assign(k0, p0);
        pe_data[pos1].assign(k1, p1);
        pe_data[pos2].assign(k2, p2);
        assert(pos0 < pos1 && pos1 < pos2);
        total_n_travs = 3;
        avg_n_travs_since_last_dist = 1;
    }

    inline void put_three_keys(const keyType *_keys, const recordPtr *_ptrs) {
        keyType k0 = _keys[0];
        keyType k1 = _keys[1];
        keyType k2 = _keys[2];
        double offset = fanout / 3.0;
        keyType s = MIN_KEY(k1 - k0, k2 - k1);
        b = offset / s;
        a = 0.5 + offset - b * k1;
        int pos0 = LR_PRED(a, b, k0, fanout);
        int pos1 = LR_PRED(a, b, k1, fanout);
        int pos2 = LR_PRED(a, b, k2, fanout);

        pe_data[pos0].assign(k0, _ptrs[0]);
        pe_data[pos1].assign(k1, _ptrs[1]);
        pe_data[pos2].assign(k2, _ptrs[2]);
        assert(pos0 < pos1 && pos1 < pos2);
        total_n_travs = 3;
        avg_n_travs_since_last_dist = 1;
    }

    void num_nonempty_stats(int &n0, int &n1, int &n2, int &n, long &total_fan, long &n_empty_slos) {
        total_fan += fanout;
        if (num_nonempty == 0) {
            n0 = n = 1;
            n1 = n2 = 0;
            return;
        } else if (num_nonempty == 1) {
            n0 = n2 = 0;
            n1 = n = 1;
        } else if (num_nonempty == 2) {
            n0 = n1 = 0;
            n2 = n = 1;
        } else {
            n = 1;
            n0 = n1 = n2 = 0;
        }
        for (int i = 0; i < fanout; ++i) {
            int cn0 = 0;
            int cn1 = 0;
            int cn2 = 0;
            int cn = 0;
            long c_total_fan = 0;
            long c_n_empty_slots = 0;
            pairEntry &pe = pe_data[i];
            if (pe.key == -1) {
                pe.child->num_nonempty_stats(cn0, cn1, cn2, cn, c_total_fan, c_n_empty_slots);
            }
            if (pe.key < -2) {
                ++n_empty_slos;
            }
            n0 += cn0;
            n1 += cn1;
            n2 += cn2;
            n += cn;
            total_fan += c_total_fan;
            n_empty_slos += c_n_empty_slots;
        }
    }


    ~diliNode(){
        if (pe_data) {
            if (fanout > 0) {
                for (int i = 0; i < fanout; ++i) {
                    if (pe_data[i].key == -1) {
                        diliNode *child = pe_data[i].child;
                        delete child;
                    }
                }
            }
            delete [] pe_data;
            pe_data = NULL;
        }
    }


//-------------------------------------------------------------------------


    void save(FILE *fp) {
        fwrite(&meta_info, sizeof(int), 1, fp);
        fwrite(&a, sizeof(double), 1, fp);
        fwrite(&b, sizeof(double), 1, fp);
        fwrite(&fanout, sizeof(int), 1, fp);
        fwrite(&num_nonempty, sizeof(int), 1, fp);
        fwrite(&avg_n_travs_since_last_dist, sizeof(double), 1, fp);
        fwrite(&total_n_travs, sizeof(long), 1, fp);


        fwrite(&last_total_n_travs, sizeof(long), 1, fp);
        fwrite(&last_nn, sizeof(int), 1, fp);
        fwrite(&n_adjust, sizeof(int), 1, fp);

        for (int i = 0; i < fanout; ++i) {
            pairEntry &pe = pe_data[i];
            keyType key = pe.key;
            fwrite(&(key), sizeof(keyType),1, fp);
            if (key >= 0) {
                fwrite(&(pe.ptr), sizeof(recordPtr),1, fp);
            } else if (key == -1){
                pe.child->save(fp);
            } else if (key == -2) {
                pe.fan2child->save(fp);
            }
        }

        if (fanout <= 0) {
            cout << "!!!!error, fanout = " << fanout << ", meta_info = " << meta_info << endl;
        }
        assert(fanout > 0);

    }

    void load(FILE *fp) {
        fread(&meta_info, sizeof(int), 1, fp);
        fread(&a, sizeof(double), 1, fp);
        fread(&b, sizeof(double), 1, fp);
        fread(&fanout, sizeof(int), 1, fp);
        assert(fanout > 0);

        fread(&num_nonempty, sizeof(int), 1, fp);
        fread(&avg_n_travs_since_last_dist, sizeof(double), 1, fp);
        fread(&total_n_travs, sizeof(long), 1, fp);


        fread(&last_total_n_travs, sizeof(long), 1, fp);
        fread(&last_nn, sizeof(int), 1, fp);
        fread(&n_adjust, sizeof(int), 1, fp);

        pe_data = new pairEntry[fanout];
        keyType key = 0;
        recordPtr ptr = 0;
        for (int i = 0; i < fanout; ++i) {
            fread(&key, sizeof(keyType), 1, fp);
            if (key >= 0) {
                fread(&ptr, sizeof(recordPtr), 1, fp);
                pe_data[i].assign(key, ptr);
            } else if (key == -1){
                diliNode *child = new diliNode(false);
                child->load(fp);
                pe_data[i].setChild(child);
            } else if (key == -2) {
                fan2Leaf *fan2child = new fan2Leaf;
                fan2child->load(fp);
                pe_data[i].setFan2Child(fan2child);
            } else {
                pe_data[i].setNull();
            }
        }
    }


    void cal_lr_params(keyType *keys, int n_keys) { //}, vector<int>& child_fans) {
        assert(n_keys >= 0);
        if (n_keys >= 2) {
            // note y_start = 1 here
//        linearReg(keys, 1, this->a, this->b, n_keys);
//        double first_pred = a + b * keys[0];
//        if (first_pred > 0) {
//            a -= (floor(first_pred) - 1);
//        }
//        double last_pred = a + b * keys[n_keys - 1];
//        fanout = std::min<int>(std::max<int>(static_cast<int>(n_keys * 1.1), n_keys + 1), ceil(last_pred));

            b = 1.0 * n_keys / (keys[n_keys - 1] - keys[0]);
            a = 1.0 - b * keys[0];
            fanout = n_keys + 2;
        } else {
            a = b = 0;
            fanout = 1;
        }
    }

    void trim() {
        if (fanout <= 0) {
            cout << "****error, fanout = " << fanout << ", is_internal = " << is_internal() << endl;
        }
        assert(fanout > 0);
        if (!is_internal() && num_nonempty == 0) {
            return;
        }
        for (int i = 0; i < fanout; ++i) {
            pairEntry &pe = pe_data[i];
            if (!is_internal() && num_nonempty == 0) {
                cout << "i = " << i << ", fan = " << fanout << ", pe.key = " << pe.key << endl;
            }
            if (pe.key == -1) {
                diliNode *child = pe.child;
                child->trim();
            }
        }
        for (int i = 0; i < fanout; ++i) {
            pairEntry &pe = pe_data[i];
            if (pe.key == -1) {
                diliNode *child = pe.child;
                assert(long(child) != -3l);
                assert(child->fanout >= 1);
                if (!(child->is_internal()) && child->num_nonempty == 0) {
                    delete child;
                    pe.setNull();
                } else if ((child->fanout == 1) || (!(child->is_internal()) && child->num_nonempty == 1)) {
                    pe_data[i] = child->pe_data[0];
                    child->fanout = 0;
                    delete child;
                }
            }
        }

    }

    void simplify() {
        for (int i = 0; i < fanout; ++i) {
            pairEntry &pe = pe_data[i];
            if (pe.key == -1) {
                diliNode *child = pe.child;
                if (child->num_nonempty == 2) {
                    pairEntry &cpe = child->pe_data[0];
                    keyType k1 = cpe.key;
                    recordPtr p1 = cpe.ptr;
                    cpe = child->pe_data[1];
                    keyType k2 = cpe.key;
                    recordPtr p2 = cpe.ptr;
                    delete child;

                    fan2Leaf *fan2child = new fan2Leaf(k1, p1, k2, p2);
                    pe.setFan2Child(fan2child);
                } else {
                    child->simplify();
                }
            }
        }
    }

    void bulk_loading(const keyType *keys, const recordPtr *ptrs, bool print) {
        if (num_nonempty == 0) {
            fanout = 2;
            total_n_travs = 0;
            return;
        } else {
            init();

            assert(fanout > 0);
            if (num_nonempty == 1) {
//                pe_data = new pairEntry[fanout];
                a = b = 0;
                pe_data[0].assign(keys[0], ptrs[0]);
                total_n_travs = 1;
                return;
            } else if (num_nonempty == 2) {
                b = 1.0 / (keys[1] - keys[0]);
                a = 0.5 - b * keys[0];
                pe_data[0].assign(keys[0], ptrs[0]);
                pe_data[1].assign(keys[1], ptrs[1]);
                total_n_travs = 2;
                return;
            } else if (num_nonempty == 3) {
                put_three_keys(keys, ptrs);
                return;
            }
        }
        distribute_data(keys, ptrs, print);
    }

    void distribute_data(const keyType *keys, const recordPtr *ptrs, bool print=false) {
        assert(num_nonempty > 3);

        total_n_travs = 0;
        linearReg_w_expanding(keys, a, b, num_nonempty, fanout, false);
//        linearReg_w_expanding(keys, a, b, num_nonempty, fanout, true);
        int last_k_id = 0;
        keyType last_key = keys[0];
        int pos = -1;
        int last_pos = LR_PRED(a, b, last_key, fanout);

        keyType final_key = keys[num_nonempty - 1];
//    int final_pos = LR_PRED(a, b, final_key, fanout);
        if (b < 0 || last_pos == LR_PRED(a, b, final_key, fanout)) {
            linearReg_w_expanding(keys, a, b, num_nonempty, fanout, true);
            last_pos = LR_PRED(a, b, last_key, fanout);
            int final_pos = LR_PRED(a, b, final_key, fanout);
            assert(last_pos != final_pos);
        }

        assert(b >= 0);

        for (int k_id = 1; k_id < num_nonempty; ++k_id) {
            keyType key = keys[k_id];
            assert (key != last_key);
            pos = LR_PRED(a, b, key, fanout);

            assert(pos >= last_pos);

            if (pos != last_pos) {
                if (k_id == last_k_id + 1) {
                    pe_data[last_pos].assign(last_key, ptrs[last_k_id]);
                    ++total_n_travs;
                } else { // need to create a new node
                    int n_keys_this_child = k_id - last_k_id;
                    if (n_keys_this_child == 3) {
                        diliNode *child = new diliNode(false);
                        child->init(3);
                        child->put_three_keys(keys + last_k_id, ptrs + last_k_id);
                        pe_data[last_pos].setChild(child);
                        total_n_travs += 6;
                    }
                    else if (n_keys_this_child == 2) {
                        fan2Leaf *fan2child = new fan2Leaf(keys[last_k_id], ptrs[last_k_id], keys[last_k_id + 1], ptrs[last_k_id + 1]);
                        pe_data[last_pos].setFan2Child(fan2child);
                        total_n_travs += 4;
                    }
                    else {
                        diliNode *child = new diliNode(false);
                        child->init(n_keys_this_child);
                        child->distribute_data(keys + last_k_id, ptrs + last_k_id);
                        pe_data[last_pos].setChild(child);
                        total_n_travs += n_keys_this_child + child->total_n_travs;
                    }
                }
                last_key = key;
                last_pos = pos;
                last_k_id = k_id;
            }
        }

        assert(last_k_id != 0);
        assert(pos >= last_pos);
        if (last_k_id == num_nonempty - 1) {
            ++total_n_travs;
            pe_data[pos].assign(keys[num_nonempty - 1], ptrs[num_nonempty - 1]);
        } else {
            int n_keys_this_child = num_nonempty - last_k_id;
            if (n_keys_this_child == 3) {

                diliNode *child = new diliNode(false);
                child->init(3);
                child->put_three_keys(keys + last_k_id, ptrs + last_k_id);
                pe_data[last_pos].setChild(child);
                total_n_travs += 6;
            }
            else if (n_keys_this_child == 2) {
                fan2Leaf *fan2child = new fan2Leaf(keys[last_k_id], ptrs[last_k_id], keys[last_k_id + 1], ptrs[last_k_id + 1]);
                pe_data[last_pos].setFan2Child(fan2child);
                total_n_travs += 4;
            }

            else {
                diliNode *child = new diliNode(false);
                child->init(n_keys_this_child);
                child->distribute_data(keys + last_k_id, ptrs + last_k_id);
                pe_data[last_pos].setChild(child);
                total_n_travs += n_keys_this_child + child->total_n_travs;
            }
        }

        last_total_n_travs = total_n_travs;
        last_nn = num_nonempty;
    }

    void cal_avg_n_travs() {
        if (!is_internal()) {
            if (num_nonempty <= 2) {
                avg_n_travs_since_last_dist = 1;
                return;
            }
            avg_n_travs_since_last_dist = 1.0 * total_n_travs / num_nonempty;
        }
        for (int i = 0; i < fanout; ++i) {
            pairEntry &pe = pe_data[i];
            if (pe.key == -1) {
                pe.child->cal_avg_n_travs();
            }
        }
    }


    inline bool insert(const keyType &_key, const recordPtr &_ptr) {
        int pred = LR_PRED(a, b, _key, fanout);
        pairEntry &pe = pe_data[pred];
//    if (print) {
//        cout << "_key = " << _key << ", pe.key = " << pe.key << ", fanout = " << fanout << ", pred = " << pred << ", num_nonempty = " << num_nonempty << endl;
//    }
        if (pe.key < -2) {
            pe.assign(_key, _ptr);
            ++num_nonempty;
            ++total_n_travs;
            if (num_nonempty >= LEAF_MAX_CAPACIY) {
                set_int_flag();
            }
            return true;
        } else if (pe.key == -1) {
            diliNode *child = pe.child;
            long child_last_total_n_travs = child->total_n_travs;
            bool if_inserted = child->insert(_key, _ptr);
#ifndef NOT_ADJUST
            if (if_inserted) {
                ++num_nonempty;
                ++total_n_travs;
                total_n_travs += (child->total_n_travs - child_last_total_n_travs);
                if (if_retrain()) {
                    collect_and_clear(dili_auxiliary::retrain_keys, dili_auxiliary::retrain_ptrs);
                    inc_n_adjust();
                    init();
                    distribute_data(dili_auxiliary::retrain_keys, dili_auxiliary::retrain_ptrs);
                    ++num_adjust_stats;
                }

                if (get_n_adjust() >= 4)  {
                    set_int_flag();
                }
            }
#endif
            return if_inserted;
        } else if (pe.key == -2) {
            total_n_travs += 2;
            fan2Leaf *fan2child = pe.fan2child;
            keyType k1 = fan2child->k1;
            recordPtr p1 = fan2child->p1;
            keyType k2 = fan2child->k2;
            recordPtr p2 = fan2child->p2;
            if (_key == k1 || _key == k2) {
                return false;
            }
            ++num_nonempty;

            diliNode *child = new diliNode(false);
            child->init(3);

            if (_key > k2) {
                child->put_three_keys(k1, p1, k2, p2, _key, _ptr);
            } else if (_key < k1) {
                child->put_three_keys(_key, _ptr, k1, p1, k2, p2);
            } else {
                child->put_three_keys(k1, p1, _key, _ptr, k2, p2);
            }

            pe.setChild(child);
            return true;
        } else if (pe.key == _key) {
            return false;
        } else {
            keyType k1, k2;
            recordPtr p1, p2;
            if (pe.key < _key) {
                k1 = pe.key;
                p1 = pe.ptr;
                k2 = _key;
                p2 = _ptr;
            } else {
                k1 = _key;
                p1 = _ptr;
                k2 = pe.key;
                p2 = pe.ptr;
            }
            assert(num_nonempty > 1);
            total_n_travs += 3;
            ++num_nonempty;

            fan2Leaf *fan2child = new fan2Leaf(k1, p1, k2, p2);
            pe.setFan2Child(fan2child);

            return true;
        }
    }


    inline int erase(const keyType &_key) {
        int pred = LR_PRED(a, b, _key, fanout);
        pairEntry &pe = pe_data[pred];
        if (pe.key == _key) {
            pe.setNull();
            --num_nonempty;
            --total_n_travs;
            return num_nonempty;
        }
        else if (pe.key == -2) {
            fan2Leaf *fan2child = pe.fan2child;
            if (fan2child->k1 == _key) {
                pe.assign(fan2child->k2, fan2child->p2);
                --num_nonempty;
                total_n_travs -= 3;
//                delete fan2child;
                return num_nonempty;
            } else if (fan2child->k2 == _key) {
                pe.assign(fan2child->k1, fan2child->p1);
                --num_nonempty;
                total_n_travs -= 3;
//                delete fan2child;
                return num_nonempty;
            } else {
                return -1;
            }
        } else if (pe.key == -1) {
            diliNode *child = pe.child;
            long child_n_travs = child->total_n_travs;
            int flag = child->erase(_key);

            if (flag > 0) {
                total_n_travs -= (child_n_travs - child->total_n_travs + 1);
                --num_nonempty;
                return num_nonempty;
            } else if (flag == 0){
                total_n_travs -= (child_n_travs + 1);
                --num_nonempty;
                pe.setNull();
                return num_nonempty;
            } else {
                return -1;
            }
        } else {
            return -1;
        }
    }
    inline int erase_and_get_ptr(const keyType &_key, recordPtr &ptr) {
        int pred = LR_PRED(a, b, _key, fanout);
        pairEntry &pe = pe_data[pred];
        if (pe.key == _key) {
            ptr = pe.ptr;
            pe.setNull();
            --num_nonempty;
            --total_n_travs;
            return num_nonempty;
        }
        else if (pe.key == -2) {
            fan2Leaf *fan2child = pe.fan2child;
            if (fan2child->k1 == _key) {
                ptr = fan2child->p1;
                pe.assign(fan2child->k2, fan2child->p2);
                --num_nonempty;
                total_n_travs -= 3;
//                delete fan2child;
                return num_nonempty;
            } else if (fan2child->k2 == _key) {
                ptr = fan2child->p2;
                pe.assign(fan2child->k1, fan2child->p1);
                --num_nonempty;
                total_n_travs -= 3;
//                delete fan2child;
                return num_nonempty;
            } else {
                return -1;
            }
        } else if (pe.key == -1) {
            diliNode *child = pe.child;
            long child_n_travs = child->total_n_travs;
            int flag = child->erase_and_get_ptr(_key, ptr);

            if (flag > 0) {
                total_n_travs -= (child_n_travs - child->total_n_travs + 1);
                --num_nonempty;
                return num_nonempty;
            } else if (flag == 0){
                total_n_travs -= (child_n_travs + 1);
                --num_nonempty;
                pe.setNull();
                return num_nonempty;
            } else {
                return -1;
            }
        } else {
            return -1;
        }
    }


    void collect_and_clear(keyType *keys, recordPtr *ptrs) {
        int j = 0;
        for(int i = 0; i < fanout; ++i) {
            pairEntry &pe = pe_data[i];
            if (pe.key >= 0) {
                keys[j] = pe.key;
                ptrs[j++] = pe.ptr;
            } else if (pe.key == -1) {
                diliNode *child = pe.child;
                child->collect_and_clear(keys+j, ptrs+j);
                j += child->num_nonempty;
                delete child;
            }
            else if (pe.key == -2) {
                fan2Leaf *fan2child = pe.fan2child;
                keys[j] = fan2child->k1;
                ptrs[j++] = fan2child->p1;
                keys[j] = fan2child->k2;
                ptrs[j++] = fan2child->p2;
            }
        }
        delete[] pe_data;
        pe_data = NULL;
        assert(j == num_nonempty);
    }

    void collect_all_keys(keyType *keys) {
        assert(b >= 0);
        int j = 0;
        for(int i = 0; i < fanout; ++i) {
            pairEntry &pe = pe_data[i];
            if (pe.key >= 0) {
                keys[j++] = pe.key;
            } else if (pe.key == -1) {
                diliNode *child = pe.child;
                child->collect_all_keys(keys+j);
                j += child->num_nonempty;
            } else if (pe.key == -2) {
                fan2Leaf *fan2child = pe.fan2child;
                keys[j++] = fan2child->k1;
                keys[j++] = fan2child->k2;
            }
        }
        if (j != num_nonempty) {
            cout << "j = " << j << ", num_nonempty = " << num_nonempty << ", is_internal = " << is_internal() << endl;
            for(int i = 0; i < num_nonempty; ++i) {
                pairEntry &pe = pe_data[i];
                if (pe.key >= 0) {
                    cout << "i = " << i << ", pe.key = " << pe.key << endl;
                } else if (pe.key == -1) {
                    diliNode *child = pe.child;
                    child->collect_all_keys(keys+j);
                    cout << "i = " << i << ", child.num_nonempty = " << child->num_nonempty << endl;
                    j += child->num_nonempty;
                } else if (pe.key == -2) {
                    fan2Leaf *fan2child = pe.fan2child;
                    cout << "i = " << i << ", fan2child.num_nonempty = 2" << endl;
                    keys[j++] = fan2child->k1;
                    keys[j++] = fan2child->k2;
                }
            }
        }
        assert(j == num_nonempty);
    }

};


#endif //DILI_DILINODE_H
