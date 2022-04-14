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
    extern std::vector<fan2Leaf*> empty_fan2leaves;
    extern std::vector<diliNode*> empty_fan2nodes;
    extern std::vector<diliNode*> empty_nodes;
    extern long *retrain_keys;
    extern long *retrain_payloads;
    void init_insert_aux_vars();
    void free_insert_aux_vars();
}


inline void linearReg_w_simple_strategy(long *X, double &a, double &b, int n) {
    int left_n = n / 2;
    int right_n = n - left_n - 1;
    long x_middle = X[left_n];
    double left_slope = 1.0 * left_n / (x_middle - X[0]);
    double right_slope = 1.0 * right_n / (X[n-1] - x_middle);
    b = MAX_DOUBLE(left_slope, right_slope);
    a = -b * x_middle + left_n;
}

inline void linearReg_at_least_four(long *X, double &a, double &b, int n) {
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

inline void linearReg_w_expanding(long *X, double &a, double &b, int n, int expanded_n, bool use_simple_strategy) {
    if (!use_simple_strategy) {
        linearReg_at_least_four(X, a, b, n);
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
    double avg_n_travs_since_last_dist;
    long total_n_travs;
    keyPayload *kp_data;

    long last_total_n_travs;
    int last_nn;
    int n_adjust;



    inline bool is_internal() { return meta_info; }
    inline void set_leaf_flag() { meta_info &= ~1U;}
    inline void set_int_flag() { meta_info |= 1; }

    inline void set_fanout(int fan) { fanout = fan;}
    inline int get_fanout() { return fanout; }

    inline void set_n_adjust(int n) { meta_info = n << 1; }

    inline int get_n_adjust() { return n_adjust; }
    inline void inc_n_adjust() { ++n_adjust; }


    inline double get_expanding_ratio() { return MIN_DOUBLE(max_expanding_ratio, 1 + 0.1 * get_n_adjust()); }

    inline void set_num_nonempty(int n) { num_nonempty = n; }

    inline void init() {
        fanout = std::max<int>(num_nonempty, minFan);
        fanout += (fanout * n_adjust) / 10;
//        fanout = std::max<int>(num_nonempty, minFan) * (1 + 0.1 * get_n_adjust());
        kp_data = new keyPayload[fanout];
    }

    diliNode(bool _is_internal): a(0), b(0), meta_info(_is_internal), fanout(0), kp_data(NULL), n_adjust(30), num_nonempty(0),
                                 total_n_travs(0), last_total_n_travs(0), last_nn(0), avg_n_travs_since_last_dist(1e10)  {}

    inline void init(const int &_num_nonempty) {
        num_nonempty = _num_nonempty;
        fanout = std::max<int>(_num_nonempty, minFan);
        fanout += (fanout * n_adjust) / 10;
//        fanout = std::max<int>(_num_nonempty, minFan) * (1 + 0.1 * get_n_adjust());
        kp_data = new keyPayload[fanout];
    }
    inline void inc_num_nonempty() { ++num_nonempty; }
    inline int get_num_nonempty() { return num_nonempty; }

    int cal_num_nonempty() {
        if (num_nonempty <= 0) {
            assert(num_nonempty == 0);
            for (int i = 0; i < fanout; ++i) {
                keyPayload &kp = kp_data[i];
                if (kp.key >= 0) {
                    ++num_nonempty;
                } else if (kp.key == -1) {
                    num_nonempty += kp.child->cal_num_nonempty();
                } else if (kp.key == -2) {
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
            keyPayload &kp = kp_data[i];
            if (kp.key == -1) {
                kp.child->init_after_bulk_load();
            }
        }
    }

    void check_num_nonempty() {
        assert(num_nonempty > 1);
        if (!is_internal()) {
            assert(num_nonempty < LEAF_MAX_CAPACIY);
        }
        for (int i = 0; i < fanout; ++i) {
            keyPayload &kp = kp_data[i];
            if (kp.key == -1) {
                kp.child->check_num_nonempty();
            }
        }
    }


    inline long leaf_find(const long &key) const {
        int pred = LR_PRED(a, b, key, fanout);
//        cout << "pred = " << pred << endl;
        keyPayload &kp = kp_data[pred];
//        cout << "kp.key = " << kp.key << endl;
        if (kp.key == key) {
            return kp.payload;
        }
        if (kp.key == -1) {
            return kp.child->leaf_find(key);
        }
        if (kp.key == -2) {
            fan2Leaf *child = kp.fan2child;
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
    inline long leaf_find_w_print(const long &key) const {
        int pred = LR_PRED(a, b, key, fanout);
        keyPayload &kp = kp_data[pred];
        cout << "key = " << key << ", num_nontempty = " << num_nonempty << ", fanout = " << fanout << ", a = " << a << ", b = " << b << ", pred = " << pred << endl;
        if (kp.key == key) {
            return kp.payload;
        }
        if (kp.key == -1) {
            return kp.child->leaf_find(key);
        }
        if (kp.key == -2) {
            fan2Leaf *child = kp.fan2child;
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

    inline int range_query_from(const long &k1, long *results) const {
        int j = 0;
        int pred = LR_PRED(a, b, k1, fanout);
        keyPayload &first_kp = kp_data[pred];
        if (first_kp.key == -1) {
            j = first_kp.child->range_query_from(k1, results);
        } else if (first_kp.key == -2) {
            fan2Leaf *fan2child = first_kp.fan2child;
            long _k1 = fan2child->k1;
            long _k2 = fan2child->k2;
            if (_k1 >= k1) {
                results[j++] = fan2child->p1;
                results[j++] = fan2child->p2;
            } else if (_k2 >= k1) {
                results[j++] = fan2child->p2;
            }
        } else if (first_kp.key >= k1) {
            results[j++] = first_kp.payload;
        }

        for(int i = pred + 1; i < fanout; ++i) {
            keyPayload &kp = kp_data[i];
            if (kp.key >= 0) {
                results[j++] = kp.payload;
            } else if (kp.key == -1) {
                kp.child->collect_all_payloads(results+j);
                j += kp.child->num_nonempty;
            } else if (kp.key == -2) {
                fan2Leaf *fan2child = kp.fan2child;
                results[j++] = fan2child->p1;
                results[j++] = fan2child->p2;
            }
        }

        return j;
    }

    inline int range_query_to(const long &k2, long *results) const {
        int j = 0;
        int pred = LR_PRED(a, b, k2, fanout);
        for(int i = 0; i < pred; ++i) {
            keyPayload &kp = kp_data[i];
            if (kp.key >= 0) {
                results[j++] = kp.payload;
            } else if (kp.key == -1) {
                kp.child->collect_all_payloads(results+j);
                j += kp.child->num_nonempty;
            } else if (kp.key == -2) {
                fan2Leaf *fan2child = kp.fan2child;
                results[j++] = fan2child->p1;
                results[j++] = fan2child->p2;
            }
        }

        keyPayload &kp = kp_data[pred];
        if (kp.key == -1) {
            j += kp.child->range_query_to(k2, results+j);
        } else if (kp.key == -2) {
            fan2Leaf *fan2child = kp.fan2child;
            long _k1 = fan2child->k1;
            long _k2 = fan2child->k2;
            if (_k2 < k2) {
                results[j++] = fan2child->p1;
                results[j++] = fan2child->p2;
            } else if (_k1 < k2) {
                results[j++] = fan2child->p1;
            }
        } else if (kp.key < k2) {
            results[j++] = kp.payload;
        }
        return j;
    }

    inline void collect_all_payloads(long *results) const{
        int j = 0;
        for(int i = 0; i < fanout; ++i) {
            keyPayload &kp = kp_data[i];
            if (kp.key >= 0) {
                results[j++] = kp.payload;
            } else if (kp.key == -1) {
                kp.child->collect_all_payloads(results+j);
                j += kp.child->num_nonempty;
            } else if (kp.key == -2) {
                fan2Leaf *fan2child = kp.fan2child;
                results[j++] = fan2child->p1;
                results[j++] = fan2child->p2;
            }
        }
    }

    inline int range_query(const long &k1, const long &k2, long *results) const {
        int pred1 = LR_PRED(a, b, k1, fanout);
        int pred2 = LR_PRED(a, b, k2, fanout);

        if (pred1 == pred2) {
            keyPayload &kp = kp_data[pred1];
            if (kp.key == -1) {
                return kp.child->range_query(k1, k2, results);
            } else if (kp.key == -2) {
                int n = 0;
                fan2Leaf *leaf = kp.fan2child;
                long _k1 = leaf->k1;
                long _k2 = leaf->k2;
                if (_k1 >= k1 && _k1 < k2) {
                    results[n++] = leaf->p1;
                }
                if (_k2 >= k1 && _k2 < k2) {
                    results[n++] = leaf->p2;
                }
                return n;
            } else if (kp.key >= k1 && kp.key < k2) {
                results[0] = kp.payload;
                return 1;
            }
        } else { // pred1 < pred2

            keyPayload &first_kp = kp_data[pred1];
            int n = 0;
            if (first_kp.key == -1) {
                n = first_kp.child->range_query_from(k1, results);
            } else if (first_kp.key == -2) {
                fan2Leaf *leaf = first_kp.fan2child;
                long _k1 = leaf->k1;
                long _k2 = leaf->k2;
                if (_k1 >= k1) {
                    results[0] = leaf->p1;
                    results[1] = leaf->p2;
                    n = 2;
                } else if (_k2 >= k1) {
                    results[1] = leaf->p2;
                    n = 1;
                }
            } else if (first_kp.key >= k1) {
                results[0] = first_kp.payload;
                n = 1;
            }

            for (int i = pred1 + 1; i < pred2; ++i) {
                keyPayload &kp = kp_data[i];
                if (kp.key == -1) {
                    kp.child->collect_all_payloads(results+n);
                    n += kp.child->num_nonempty;
                } else if (kp.key == -2) {
                    fan2Leaf *leaf = kp.fan2child;
                    results[n++] = leaf->p1;
                    results[n++] = leaf->p2;
                } else if (kp.key >= k1 && kp.key < k2) {
                    results[n++] = kp.payload;
                }
            }

            keyPayload &final_kp = kp_data[pred2];
            if (final_kp.key == -1) {
                n += final_kp.child->range_query_to(k2, results+n);
            } else if (final_kp.key == -2) {
                fan2Leaf *leaf = final_kp.fan2child;
                long _k1 = leaf->k1;
                long _k2 = leaf->k2;
                if (_k2 < k2) {
                    results[n++] = leaf->p1;
                    results[n++] = leaf->p2;
                } else if (_k1 < k2) {
                    results[n++] = leaf->p1;
                }

            } else if (final_kp.key < k2) {
                results[n++] = final_kp.payload;
            }
            return n;
        }
    }


    inline diliNode* find_child(const long &key) {
        int i = LR_PRED(a, b, key, fanout);
        return kp_data[i].child;
    }

    inline diliNode* find_child_w_print(const long &key) {
        int i = LR_PRED(a, b, key, fanout);
        cout << "key = " << key << ", pred = " << i << endl;
        return kp_data[i].child;
    }

    inline int get_child_id(const long &key) { return LR_PRED(a, b, key, fanout); }
    inline bool if_retrain() { return (!is_internal()) && (total_n_travs * last_nn >= ((last_total_n_travs * num_nonempty) << 1) ); }
//    inline bool if_retrain() { return 1.0 * total_n_travs / num_nonempty >= 2 * avg_n_travs_since_last_dist; }

    inline void put_two_keys(long k0, long p0, long k1, long p1) {
        double offset = fanout / 3.0;
        b = offset / (k1 - k0);
        a = offset - b * k0;

//        b = 1.0 / (k1 - k0);
//        a = 0.5 - b * k0;

        int pos0 = LR_PRED(a, b, k0, fanout);
        int pos1 = LR_PRED(a, b, k1, fanout);
        kp_data[pos0].assign(k0, p0);
        kp_data[pos1].assign(k1, p1);
        assert(pos0 < pos1);
        total_n_travs = 2;
//        avg_n_travs_since_last_dist = 1;
    }

    inline void put_three_keys(const long &k0, const long &p0, const long &k1, const long &p1, const long &k2, const long &p2) {
        double offset = fanout / 3.0;
        long s = MIN_LONG(k1 - k0, k2 - k1);
        b = offset / s;
        a = 0.5 + offset - b * k1;
        int pos0 = LR_PRED(a, b, k0, fanout);
        int pos1 = LR_PRED(a, b, k1, fanout);
        int pos2 = LR_PRED(a, b, k2, fanout);
        kp_data[pos0].assign(k0, p0);
        kp_data[pos1].assign(k1, p1);
        kp_data[pos2].assign(k2, p2);
        assert(pos0 < pos1 && pos1 < pos2);
        total_n_travs = 3;
        avg_n_travs_since_last_dist = 1;
    }

    inline void put_three_keys(long *_keys, long *_payloads) {
        long k0 = _keys[0];
        long k1 = _keys[1];
        long k2 = _keys[2];
        double offset = fanout / 3.0;
        long s = MIN_LONG(k1 - k0, k2 - k1);
        b = offset / s;
        a = 0.5 + offset - b * k1;
        int pos0 = LR_PRED(a, b, k0, fanout);
        int pos1 = LR_PRED(a, b, k1, fanout);
        int pos2 = LR_PRED(a, b, k2, fanout);

        kp_data[pos0].assign(k0, _payloads[0]);
        kp_data[pos1].assign(k1, _payloads[1]);
        kp_data[pos2].assign(k2, _payloads[2]);
        assert(pos0 < pos1 && pos1 < pos2);
        total_n_travs = 3;
        avg_n_travs_since_last_dist = 1;
    }



    inline bool order_check() const{

        long last_key = 0;
        int i = 0;
        for (; i < fanout; ++i) {
            long key = kp_data[i].key;
            if (key >= 0) {
                last_key = key;
                break;
            }
        }
        ++i;
        for (; i < fanout; ++i) {
            long key = kp_data[i].key;
            if (key >= 0) {
                if (last_key >= key) {
                    cout << "error at " << i << ", last_key = " << last_key << ", keys[i] = " << key << endl;
                    return false;
                }
            }
        }
        return true;
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
            keyPayload &kp = kp_data[i];
            if (kp.key == -1) {
                kp.child->num_nonempty_stats(cn0, cn1, cn2, cn, c_total_fan, c_n_empty_slots);
            }
            if (kp.key < -2) {
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
        if (kp_data) {
            if (fanout > 0) {
                for (int i = 0; i < fanout; ++i) {
                    if (kp_data[i].key == -1) {
                        diliNode *child = kp_data[i].child;
                        delete child;
                    }
                }
            }
            delete [] kp_data;
            kp_data = NULL;
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
            keyPayload &kp = kp_data[i];
            long key = kp.key;
            fwrite(&(key), sizeof(long),1, fp);
            if (key >= 0) {
                fwrite(&(kp.payload), sizeof(long),1, fp);
            } else if (key == -1){
                kp.child->save(fp);
            } else if (key == -2) {
                kp.fan2child->save(fp);
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

        kp_data = new keyPayload[fanout];
        long key = 0;
        long payload = 0;
        for (int i = 0; i < fanout; ++i) {
            fread(&key, sizeof(long), 1, fp);
            if (key >= 0) {
                fread(&payload, sizeof(long), 1, fp);
                kp_data[i].assign(key, payload);
            } else if (key == -1){
                diliNode *child = new diliNode(false);
                child->load(fp);
                kp_data[i].setChild(child);
            } else if (key == -2) {
                fan2Leaf *fan2child = new fan2Leaf;
                fan2child->load(fp);
                kp_data[i].setFan2Child(fan2child);
            } else {
                kp_data[i].setNull();
            }
        }
    }


    void cal_lr_params(long *keys, int n_keys) { //}, vector<int>& child_fans) {
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
            keyPayload &kp = kp_data[i];
            if (!is_internal() && num_nonempty == 0) {
                cout << "i = " << i << ", fan = " << fanout << ", kp.key = " << kp.key << endl;
            }
            if (kp.key == -1) {
                diliNode *child = kp.child;
                child->trim();
            }
        }
        for (int i = 0; i < fanout; ++i) {
            keyPayload &kp = kp_data[i];
            if (kp.key == -1) {
                diliNode *child = kp.child;
                assert(long(child) != -3l);
                assert(child->fanout >= 1);
                if (!(child->is_internal()) && child->num_nonempty == 0) {
                    delete child;
                    kp.setNull();
                } else if ((child->fanout == 1) || (!(child->is_internal()) && child->num_nonempty == 1)) {
                    kp_data[i] = child->kp_data[0];
                    child->fanout = 0;
                    delete child;
                }
            }
        }

    }

    void simplify() {
        for (int i = 0; i < fanout; ++i) {
            keyPayload &kp = kp_data[i];
            if (kp.key == -1) {
                diliNode *child = kp.child;
                if (child->num_nonempty == 2) {
                    keyPayload &ckp = child->kp_data[0];
                    long k1 = ckp.key;
                    long p1 = ckp.payload;
                    ckp = child->kp_data[1];
                    long k2 = ckp.key;
                    long p2 = ckp.payload;
                    delete child;

                    fan2Leaf *fan2child = new fan2Leaf(k1, p1, k2, p2);
                    kp.setFan2Child(fan2child);
                } else {
                    child->simplify();
                }
            }
        }
    }

    void bulk_loading(long *keys, long *payloads, bool print) {
        if (num_nonempty == 0) {
            fanout = 2;
            total_n_travs = 0;
            return;
        } else {
            init();
            assert(fanout > 0);
            if (num_nonempty == 1) {
                kp_data = new keyPayload[fanout];
                a = b = 0;
                kp_data[0].assign(keys[0], payloads[0]);
                total_n_travs = 1;
                return;
            } else if (num_nonempty == 2) {
#ifndef ALLOW_FAN2_NODE
                b = 1.0 / (keys[1] - keys[0]);
                a = 0.5 - b * keys[0];
                kp_data[0].assign(keys[0], payloads[0]);
                kp_data[1].assign(keys[1], payloads[1]);
                total_n_travs = 2;
#else
                put_two_keys(keys[0], payloads[0], keys[1], payloads[1]);
#endif
                return;
            } else if (num_nonempty == 3) {
                put_three_keys(keys, payloads);
                return;
            }
        }

        distribute_data(keys, payloads, print);
    }


    void distribute_data(long *keys, long *payloads, bool print=false) {
        if (print || fanout <= 0) {
            cout << "fanout = " << fanout << ", num_nonempty = " << num_nonempty << ", minFan = " << minFan
                 << ", expanding_ratio = " << get_expanding_ratio() << ", max_expanding_ratio = " << max_expanding_ratio
                 << ", n_adjust = " << get_n_adjust() << endl;
        }
        assert(num_nonempty > 3);

        total_n_travs = 0;
        linearReg_w_expanding(keys, a, b, num_nonempty, fanout, false);
        int last_k_id = 0;
        long last_key = keys[0];
        int pos = -1;
        int last_pos = LR_PRED(a, b, last_key, fanout);

        long final_key = keys[num_nonempty - 1];
//    int final_pos = LR_PRED(a, b, final_key, fanout);
        if (b < 0 || last_pos == LR_PRED(a, b, final_key, fanout)) {
            linearReg_w_expanding(keys, a, b, num_nonempty, fanout, true);
            last_pos = LR_PRED(a, b, last_key, fanout);
            int final_pos = LR_PRED(a, b, final_key, fanout);

            if (last_pos == final_pos) {
                for (int i = 0; i < num_nonempty; ++i) {
                    cout << "i = " << i << ", key = " << keys[i] << ", doublepred = " << a + b * keys[i] << ", pred = " <<  LR_PRED(a, b, keys[i], fanout) << endl;
                }
                cout << "first_key = " << last_key << ", final_key = " << final_key << ", first_pos = " << last_pos << ", final_pos = " << final_pos << ", num_nonempty = " << num_nonempty << endl;
            }
            assert(last_pos != final_pos);
        }
        if (b < 0) {
            cout << "error!!!!!!b = " << b << ", fanout = " << fanout << ", num_nonempty = " << num_nonempty << endl;
            for (int _i = 1; _i < num_nonempty; ++_i) {
                assert(keys[_i] > keys[_i-1]);
            }
        }
        assert(b >= 0);

        for (int k_id = 1; k_id < num_nonempty; ++k_id) {
            long key = keys[k_id];
            assert (key != last_key);
            pos = LR_PRED(a, b, key, fanout);
//            if (print) {
//                cout << "k_id = " << k_id << ", key = " << key << ", pos = " << pos << ", last_pos = " << last_pos << ", last_k_id = " << last_k_id << endl;
//            }
            if (pos < last_pos) {
                cout << "**********check**********" << endl;
                for (int _i = 1; _i < num_nonempty; ++_i) {
                    assert(keys[_i] > keys[_i-1]);
                }
                cout << "**********check finished**********" << endl;
            }
            assert(pos >= last_pos);

            if (pos != last_pos) {
                if (k_id == last_k_id + 1) {
                    kp_data[last_pos].assign(last_key, payloads[last_k_id]);
                    ++total_n_travs;
                } else { // need to create a new node
                    int n_keys_this_child = k_id - last_k_id;
                    if (n_keys_this_child == 3) {
//                        diliNode *child;
//                        if (!(dili_auxiliary::empty_nodes.empty())) {
//                            child = dili_auxiliary::empty_nodes.back();
//                            dili_auxiliary::empty_nodes.pop_back();
//                            assert(child->num_nonempty == 3);
//                        } else {
//                            child = new diliNode(false);
//                            child->init(3);
//                        }
                        diliNode *child = new diliNode(false);
                        child->init(3);
                        child->put_three_keys(keys + last_k_id, payloads + last_k_id);
                        kp_data[last_pos].setChild(child);
                        total_n_travs += 6;
                    }

#ifndef ALLOW_FAN2_NODE
                    else if (n_keys_this_child == 2) {
                        fan2Leaf *fan2child;
                        if (!(dili_auxiliary::empty_fan2leaves.empty())) {
                            fan2child = dili_auxiliary::empty_fan2leaves.back();
                            dili_auxiliary::empty_fan2leaves.pop_back();
                        } else {
                            fan2child = new fan2Leaf();
                        }
                        fan2child->k1 = keys[last_k_id];
                        fan2child->p1 = payloads[last_k_id];
                        fan2child->k2 = keys[last_k_id + 1];
                        fan2child->p2 = payloads[last_k_id + 1];
                        kp_data[last_pos].setFan2Child(fan2child);
                        total_n_travs += 4;
                    }
#else
                        else if (n_keys_this_child == 2) {
//                    diliNode *fan2node;
//                    if (!(dili_auxiliary::empty_fan2nodes.empty())) {
//                        fan2node = dili_auxiliary::empty_fan2nodes.back();
//                        dili_auxiliary::empty_fan2nodes.pop_back();
//                        assert(fan2node->num_nonempty == 2);
//                    } else {
//                        fan2node = new diliNode(false);
//                        fan2node->init(2);
//                    }
                    diliNode *fan2node = new diliNode(false);
                    fan2node->init(2);
                    fan2node->put_two_keys(keys[last_k_id], payloads[last_k_id], keys[last_k_id + 1],
                                           payloads[last_k_id + 1]);
                    kp_data[last_pos].setChild(fan2node);
                    total_n_travs += 4;
                }
#endif
                    else {
                        diliNode *child = new diliNode(false);
                        child->init(n_keys_this_child);
                        child->distribute_data(keys + last_k_id, payloads + last_k_id);
                        kp_data[last_pos].setChild(child);
                        total_n_travs += n_keys_this_child + child->total_n_travs;
                    }
                }
                last_key = key;
                last_pos = pos;
                last_k_id = k_id;
            }
        }

//        if (print) {
//            cout << "final_key = " << keys[num_nonempty-1] << ", pos = " << pos << ", last_pos = " << last_pos << ", last_k_id = " << last_k_id << ", n_keys_this_child = " << num_nonempty - last_k_id << endl;
//        }
        assert(last_k_id != 0);
        assert(pos >= last_pos);
        if (last_k_id == num_nonempty - 1) {
            ++total_n_travs;
            kp_data[pos].assign(keys[num_nonempty - 1], payloads[num_nonempty - 1]);
        } else {
//        diliNode *child = new diliNode(false);
//        child->set_leaf_flag();
//        child->set_num_nonempty(n_keys_this_child);
//        child->bulk_loading(keys + last_k_id, payloads + last_k_id);
//        kp_data[pos].setChild(child);
//        total_n_travs += n_keys_this_child + child->total_n_travs;

            int n_keys_this_child = num_nonempty - last_k_id;
            if (n_keys_this_child == 3) {
//                diliNode *child;
//                if (!(dili_auxiliary::empty_nodes.empty())) {
//                    child = dili_auxiliary::empty_nodes.back();
//                    dili_auxiliary::empty_nodes.pop_back();
//                    assert(child->num_nonempty == 3);
//                } else {
//                    child = new diliNode(false);
//                    child->init(3);
//                }

                diliNode *child = new diliNode(false);
                child->init(3);
                child->put_three_keys(keys + last_k_id, payloads + last_k_id);
                kp_data[last_pos].setChild(child);
                total_n_travs += 6;
            }

#ifndef ALLOW_FAN2_NODE
            else if (n_keys_this_child == 2) {
                fan2Leaf *fan2child;
                if (!(dili_auxiliary::empty_fan2leaves.empty())) {
                    fan2child = dili_auxiliary::empty_fan2leaves.back();
                    dili_auxiliary::empty_fan2leaves.pop_back();
                } else {
                    fan2child = new fan2Leaf();
                }
                fan2child->k1 = keys[last_k_id];
                fan2child->p1 = payloads[last_k_id];
                fan2child->k2 = keys[last_k_id + 1];
                fan2child->p2 = payloads[last_k_id + 1];
                kp_data[last_pos].setFan2Child(fan2child);
                total_n_travs += 4;
            }
#else
            else if (n_keys_this_child == 2) {
//                diliNode *fan2node;
//                if (!(dili_auxiliary::empty_fan2nodes.empty())) {
//                    fan2node = dili_auxiliary::empty_fan2nodes.back();
//                    dili_auxiliary::empty_fan2nodes.pop_back();
//                    assert(fan2node->num_nonempty == 2);
//                } else {
//                    fan2node = new diliNode(false);
//                    fan2node->init(2);
//                }

                diliNode *fan2node = new diliNode(false);
                fan2node->init(2);
                fan2node->put_two_keys(keys[last_k_id], payloads[last_k_id], keys[last_k_id + 1],
                                       payloads[last_k_id + 1]);
                kp_data[last_pos].setChild(fan2node);
                total_n_travs += 4;
            }
#endif
            else {
                diliNode *child = new diliNode(false);
                child->init(n_keys_this_child);
                child->distribute_data(keys + last_k_id, payloads + last_k_id);
                kp_data[last_pos].setChild(child);
                total_n_travs += n_keys_this_child + child->total_n_travs;
            }
        }

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
            keyPayload &kp = kp_data[i];
            if (kp.key == -1) {
                kp.child->cal_avg_n_travs();
            }
        }
    }
    bool collect_and_check(long x0) {
        int j = 0;
        long last_key = x0;
        for(int i = 0; i < fanout; ++i) {
            keyPayload &kp = kp_data[i];
            if (kp.key >= 0) {
                if (kp.key <= last_key) {
                    cout << "error!!!!! case 1, x0 = " << x0 << ", last_key = " << last_key << ", key = " << kp.key << ", num_nonempty = " << num_nonempty << endl;
                    return false;
                }
                last_key = kp.key;
                ++j;
            } else if (kp.key == -1) {
                diliNode *child = kp.child;
                child->collect_and_check(last_key);
                j += child->num_nonempty;
            } else if (kp.key == -2) {
                j += 2;
                fan2Leaf *fan2child = kp.fan2child;
                long k1 = fan2child->k1;
                long k2 = fan2child->k2;
                if (k1 >= k2 || last_key >= k1) {
                    cout << "error!!!!! case 2, x0 = " << x0 << ", last_key = " << last_key << ", k1 = " << k1 << ", k2 = " << k2 << ", num_nonempty = " << num_nonempty <<  endl;

                    for(int t = 0; t < fanout; ++t) {
                        keyPayload &kp = kp_data[t];
                        if (kp.key >= 0) {
                            cout << kp.key << " ";
                        }
                    }
                    cout << endl;
                    return false;
                }
            }
        }
        if (j != num_nonempty) {
            cout << "error!!!!!j = " << j << ", num_nonempty = " << num_nonempty << endl;
        }
        return (j == num_nonempty);
    }


#ifndef ALLOW_FAN2_NODE
    inline bool insert(const long &_key, const long &_payload) {
        int pred = LR_PRED(a, b, _key, fanout);
        keyPayload &kp = kp_data[pred];
//    if (print) {
//        cout << "_key = " << _key << ", kp.key = " << kp.key << ", fanout = " << fanout << ", pred = " << pred << ", num_nonempty = " << num_nonempty << endl;
//    }
        if (kp.key < -2) {
            kp.assign(_key, _payload);
            ++num_nonempty;
            ++total_n_travs;
            if (num_nonempty >= LEAF_MAX_CAPACIY) {
                set_int_flag();
            }
            return true;
        } else if (kp.key == -1) {
            diliNode *child = kp.child;
            long child_last_total_n_travs = child->total_n_travs;
            bool if_inserted = child->insert(_key, _payload);
            if (if_inserted) {
                ++num_nonempty;
                ++total_n_travs;
                total_n_travs += (child->total_n_travs - child_last_total_n_travs);
//                if (num_nonempty >= LEAF_MAX_CAPACIY) {
//                    set_int_flag();
//                }


//                if (if_retrain()) {
////                    collect_and_clear(dili_auxiliary::retrain_keys, dili_auxiliary::retrain_payloads);
//                    inc_n_adjust();
////                    distribute_data(dili_auxiliary::retrain_keys, dili_auxiliary::retrain_payloads);
//                }
            }
            return if_inserted;
        } else if (kp.key == -2) {
//        not define ALLOW_FAN2_NODE
            total_n_travs += 2;
            fan2Leaf *fan2child = kp.fan2child;
            long k1 = fan2child->k1;
            long p1 = fan2child->p1;
            long k2 = fan2child->k2;
            long p2 = fan2child->p2;
            if (_key == k1 || _key == k2) {
                return false;
            }
            ++num_nonempty;
//            if (num_nonempty >= LEAF_MAX_CAPACIY) {
//                set_int_flag();
//            }


//        if (print) {
//            cout << "****empty_nodes.size() = " << dili_auxiliary::empty_nodes.size() << ", k1 = " << k1 << ", k2 = " << k2 << ", key = " << _key << endl;
//        }
            /*
            diliNode *child;
            if (!(dili_auxiliary::empty_nodes.empty())) {
                child = dili_auxiliary::empty_nodes.back();
                dili_auxiliary::empty_nodes.pop_back();
                assert(child->num_nonempty == 3);
            } else {
                child = new diliNode(false);
                child->init(3);
            }*/
            diliNode *child = new diliNode(false);
            child->init(3);

            if (_key > k2) {
                child->put_three_keys(k1, p1, k2, p2, _key, _payload);
            } else if (_key < k1) {
                child->put_three_keys(_key, _payload, k1, p1, k2, p2);
            } else {
                child->put_three_keys(k1, p1, _key, _payload, k2, p2);
            }

            kp.setChild(child);
//            dili_auxiliary::empty_fan2leaves.push_back(fan2child);
            return true;
        } else if (kp.key == _key) {
            return false;
        } else {
            long k1, k2;
            long p1, p2;
            if (kp.key < _key) {
                k1 = kp.key;
                p1 = kp.payload;
                k2 = _key;
                p2 = _payload;
            } else {
                k1 = _key;
                p1 = _payload;
                k2 = kp.key;
                p2 = kp.payload;
            }
            assert(num_nonempty > 1);
            total_n_travs += 3;
            ++num_nonempty;


//            if(num_nonempty >= LEAF_MAX_CAPACIY) {
//                set_int_flag();
//            }


//    if (if_retrain()) {
//        int tmp = 1;
//        collect_and_clear(dili_auxiliary::retrain_keys, dili_auxiliary::retrain_payloads, k1, p1, k2, p2, pred);
//        inc_n_adjust();
//        distribute_data(dili_auxiliary::retrain_keys, dili_auxiliary::retrain_payloads);
//    } else {

//            fan2Leaf *fan2child;
//            if (!(dili_auxiliary::empty_fan2leaves.empty())) {
//                fan2child = dili_auxiliary::empty_fan2leaves.back();
//                dili_auxiliary::empty_fan2leaves.pop_back();
//                fan2child->set(k1, p1, k2, p2);
//            } else {
//                fan2child = new fan2Leaf(k1, p1, k2, p2);
//            }

            fan2Leaf *fan2child = new fan2Leaf(k1, p1, k2, p2);
            kp.setFan2Child(fan2child);

            return true;
        }
    }
#else
    inline bool insert(const long &_key, const long &_payload) {

        int pred = LR_PRED(a, b, _key, fanout);
        keyPayload &kp = kp_data[pred];
        if (kp.key < -2) {
            kp.assign(_key, _payload);
            ++num_nonempty;
            ++total_n_travs;

            if (num_nonempty >= LEAF_MAX_CAPACIY) {
                set_int_flag();
            }

            return true;
        } else if (kp.key == -1) {
            diliNode *child = kp.child;
            long child_last_total_n_travs = child->total_n_travs;
            bool if_inserted = child->insert(_key, _payload);
            if (if_inserted) {
                ++num_nonempty;
                ++total_n_travs;
                total_n_travs += (child->total_n_travs - child_last_total_n_travs);
                if (num_nonempty >= LEAF_MAX_CAPACIY) {
                    set_int_flag();
                }

                if(!is_internal() && (num_nonempty > (last_nn << 1))) {
                    if (if_retrain()) {
//                    collect_and_clear(dili_auxiliary::retrain_keys, dili_auxiliary::retrain_payloads);
                        inc_n_adjust();
//                    distribute_data(dili_auxiliary::retrain_keys, dili_auxiliary::retrain_payloads);
                    }
                }



            }
            return if_inserted;
        }  else {
            long k1, k2;
            long p1, p2;
            if (kp.key < _key) {
                k1 = kp.key;
                p1 = kp.payload;
                k2 = _key;
                p2 = _payload;
            } else {
                k1 = _key;
                p1 = _payload;
                k2 = kp.key;
                p2 = kp.payload;
            }
//            assert(num_nonempty > 1);
            total_n_travs += 3;
            ++num_nonempty;

            if(num_nonempty >= LEAF_MAX_CAPACIY) {
                set_int_flag();
            }


//    if (if_retrain()) {
//        int tmp = 1;
//        collect_and_clear(dili_auxiliary::retrain_keys, dili_auxiliary::retrain_payloads, k1, p1, k2, p2, pred);
//        inc_n_adjust();
//        distribute_data(dili_auxiliary::retrain_keys, dili_auxiliary::retrain_payloads);
//    } else {

//            diliNode *fan2node;
//            if (!(dili_auxiliary::empty_fan2nodes.empty())) {
//                fan2node = dili_auxiliary::empty_fan2nodes.back();
//                dili_auxiliary::empty_fan2nodes.pop_back();
//                assert(fan2node->num_nonempty == 2);
//            } else {
//                fan2node = new diliNode(false);
//                fan2node->init(2);
//            }

            diliNode *fan2node = new diliNode(false);
            fan2node->init(2);
            fan2node->put_two_keys(k1, p1, k2, p2);
            kp.setChild(fan2node);
            return true;
        }
    }
#endif

    inline bool leaf_insert(const long &_key, const long &_payload) {
        int pred = LR_PRED(a, b, _key, fanout);
        keyPayload &kp = kp_data[pred];
//        if (_key == 957900359748l) {
//            cout << "key = " << _key << ", is_internal = " << is_internal() << ", pred = " << pred << ", num_nonempty = " << num_nonempty << ", fanout = "
//                 << fanout << ", kp.key = " << kp.key << ", empty_fan2leaves.size = " << dili_auxiliary::empty_fan2leaves.size() << endl;
//        }
        if (kp.key < -2) {
            kp.assign(_key, _payload);
            ++num_nonempty;
            ++total_n_travs;
            if (num_nonempty >= LEAF_MAX_CAPACIY) {
                set_int_flag();
            }
            return true;
        } else if (kp.key == -1) {
            diliNode *child = kp.child;
            long child_last_total_n_travs = child->total_n_travs;
            bool if_inserted = child->leaf_insert(_key, _payload);
            if (if_inserted) {
                ++num_nonempty;
                ++total_n_travs;
                total_n_travs += (child->total_n_travs - child_last_total_n_travs);
                if (num_nonempty >= LEAF_MAX_CAPACIY) {
                    set_int_flag();
                }
//                if (if_retrain()) {
//                    cout << "+++++1. prepare to retrain. key = " << _key << ", num_nonempty = " << num_nonempty
//                            << ", fanout = " << fanout << ", pred = " << pred << endl;
//                    collect_and_clear(dili_auxiliary::retrain_keys, dili_auxiliary::retrain_payloads);
//                    cout << "+++++2. here is OK." << endl;
//                    inc_n_adjust();
//                    cout << "+++++3. here is OK." << endl;
//                    distribute_data(dili_auxiliary::retrain_keys, dili_auxiliary::retrain_payloads);
//                    cout << "+++++4. here is OK." << endl;
//                }
            }
            return if_inserted;
        } else if (kp.key == -2) {
            total_n_travs += 2;
            fan2Leaf *fan2child = kp.fan2child;
            long k1 = fan2child->k1;
            long p1 = fan2child->p1;
            long k2 = fan2child->k2;
            long p2 = fan2child->p2;
            if (_key == k1 || _key == k2) {
                return false;
            }
            ++num_nonempty;
            if (num_nonempty >= LEAF_MAX_CAPACIY) {
                set_int_flag();
            }
            /*
            diliNode *child;
            if (!(dili_auxiliary::empty_nodes.empty())) {
                child = dili_auxiliary::empty_nodes.back();
                dili_auxiliary::empty_nodes.pop_back();
                assert(child->num_nonempty == 3);
            } else {
                child = new diliNode(false);
                child->init(3);
            }*/
            diliNode *child = new diliNode(false);
            child->init(3);

            if (_key > k2) {
                child->put_three_keys(k1, p1, k2, p2, _key, _payload);
            } else if (_key < k1) {
                child->put_three_keys(_key, _payload, k1, p1, k2, p2);
            } else {
                child->put_three_keys(k1, p1, _key, _payload, k2, p2);
            }

            kp.setChild(child);
            dili_auxiliary::empty_fan2leaves.push_back(fan2child);
            return true;
        } else if (kp.key == _key) {
            return false;
        } else {
            long k1, k2;
            long p1, p2;
            if (kp.key < _key) {
                k1 = kp.key;
                p1 = kp.payload;
                k2 = _key;
                p2 = _payload;
            } else {
                k1 = _key;
                p1 = _payload;
                k2 = kp.key;
                p2 = kp.payload;
            }
            assert(num_nonempty > 1);
            total_n_travs += 3;
            ++num_nonempty;
            if(num_nonempty >= LEAF_MAX_CAPACIY) {
                set_int_flag();
            }

#ifndef ALLOW_FAN2_NODE
            fan2Leaf *fan2child;
            if (!(dili_auxiliary::empty_fan2leaves.empty())) {
                fan2child = dili_auxiliary::empty_fan2leaves.back();
                dili_auxiliary::empty_fan2leaves.pop_back();
                fan2child->set(k1, p1, k2, p2);
            } else {
                fan2child = new fan2Leaf(k1, p1, k2, p2);
            }

//        fan2Leaf *fan2child = new fan2Leaf(k1, p1, k2, p2);
            kp.setFan2Child(fan2child);
#else
            diliNode *fan2node;
            if (!(dili_auxiliary::empty_fan2nodes.empty())) {
                fan2node = dili_auxiliary::empty_fan2nodes.back();
                dili_auxiliary::empty_fan2nodes.pop_back();
                assert(fan2node->num_nonempty == 2);
            } else {
                fan2node = new diliNode(false);
                fan2node->init(2);
            }

            fan2node->put_two_keys(k1, p1, k2, p2);
            kp.setChild(fan2node);
#endif
            return true;
        }
    }




#ifndef ALLOW_FAN2_NODE
    fan2Leaf* convert_to_fan2leaf() { // num_nonemtpy == 2
        long k1, p1;
        long k2, p2;
        bool flag = true;
        for (int i = 0; i < fanout; ++i) {
            keyPayload &kp = kp_data[i];
            if (kp.key == -2) {
                fan2Leaf *leaf = kp.fan2child;
                kp.setNull();
                return leaf;
            } else if (kp.key >= 0) {
                if (flag) {
                    k1 = kp.key;
                    p1 = kp.payload;
                    flag = false;
                } else {
                    k2 = kp.key;
                    p2 = kp.payload;
                    return new fan2Leaf(k1, p1, k2, p2);
                }
            }
        }
        return NULL;
    }

    inline int erase(const long &_key) {
        int pred = LR_PRED(a, b, _key, fanout);
        keyPayload &kp = kp_data[pred];
        if (kp.key == _key) {
            kp.setNull();
            --num_nonempty;
            --total_n_travs;
            return num_nonempty;
        }
        else if (kp.key == -2) {
            fan2Leaf *fan2child = kp.fan2child;
            if (fan2child->k1 == _key) {
                kp.assign(fan2child->k2, fan2child->p2);
                --num_nonempty;
                total_n_travs -= 3;
//                delete fan2child;
                return num_nonempty;
            } else if (fan2child->k2 == _key) {
                kp.assign(fan2child->k1, fan2child->p1);
                --num_nonempty;
                total_n_travs -= 3;
//                delete fan2child;
                return num_nonempty;
            } else {
                return -1;
            }
        } else if (kp.key == -1) {
            diliNode *child = kp.child;
            long child_n_travs = child->total_n_travs;
            int flag = child->erase(_key);
//            if (flag > 2) {
//                total_n_travs -= (child_n_travs - child->total_n_travs + 1);
//                --num_nonempty;
//                return num_nonempty;
//            } else if (flag == 2) { // adjust
//                fan2Leaf *fan2child = child->convert_to_fan2leaf();
//                delete child;
//                total_n_travs -= (child_n_travs - 1); // total_n_travs -= (child_n_travs + 1 - 2);
//                kp.setFan2Child(fan2child);
//                --num_nonempty;
//                return num_nonempty;
//            } else {
//                return -1;
//            }
            if (flag > 0) {
                total_n_travs -= (child_n_travs - child->total_n_travs + 1);
                --num_nonempty;
                return num_nonempty;
            } else if (flag == 0){
                total_n_travs -= (child_n_travs + 1);
                --num_nonempty;
                kp.setNull();
                return num_nonempty;
            } else {
                return -1;
            }
        }
    }
#else
    inline int erase(const long &_key) {
        int pred = LR_PRED(a, b, _key, fanout);
        keyPayload &kp = kp_data[pred];
        if (kp.key == _key) {
            kp.setNull();
            --num_nonempty;
            --total_n_travs;
            return num_nonempty;
        }
        else if (kp.key == -1) {
            diliNode *child = kp.child;
            long child_n_travs = child->total_n_travs;
            int flag = kp.child->erase(_key);
            if (flag > 1) {
                total_n_travs -= (child_n_travs - child->total_n_travs + 1);
                --num_nonempty;
                return num_nonempty;
            } else if (flag == 1) {
                total_n_travs -= (child_n_travs + 1);
                diliNode *fan1_child = kp.child;
                keyPayload *child_kps = fan1_child->kp_data;
                int child_fan = fan1_child->fanout;
                for (int i = 0; i < child_fan; ++i) {
                    if (child_kps[i].key >= 0) {
                        kp = child_kps[i];
                        break;
                    }
                }
                delete fan1_child;
                --num_nonempty;
                return num_nonempty;
            } else {
                return -1;
            }
        }
        return -1;
    }
#endif

    void collect_and_clear(long *keys, long *payloads) {
        int j = 0;
        for(int i = 0; i < fanout; ++i) {
            keyPayload &kp = kp_data[i];
            if (kp.key >= 0) {
                keys[j] = kp.key;
                payloads[j++] = kp.payload;
            } else if (kp.key == -1) {
                diliNode *child = kp.child;
                child->collect_and_clear(keys+j, payloads+j);
                j += child->num_nonempty;
                delete child;
            }
#ifndef ALLOW_FAN2_NODE
            else if (kp.key == -2) {
                fan2Leaf *fan2child = kp.fan2child;
                keys[j] = fan2child->k1;
                payloads[j++] = fan2child->p1;
                keys[j] = fan2child->k2;
                payloads[j++] = fan2child->p2;
                dili_auxiliary::empty_fan2leaves.push_back(fan2child);
            }
#endif
        }
        delete[] kp_data;
        kp_data = NULL;
        assert(j == num_nonempty);
    }

    void collect_and_clear(long *keys, long *payloads, long k1, long p1, long k2, long p2, int pos) {
        int j = 0;
        for(int i = 0; i < pos; ++i) {
            keyPayload &kp = kp_data[i];
            if (kp.key >= 0) {
                keys[j] = kp.key;
                payloads[j++] = kp.payload;
            } else if (kp.key == -1) {
                diliNode *child = kp.child;
                child->collect_and_clear(keys+j, payloads+j);
                j += child->num_nonempty;
                delete child;
            }
#ifndef ALLOW_FAN2_NODE
            else if (kp.key == -2) {
                fan2Leaf *fan2child = kp.fan2child;
                keys[j] = fan2child->k1;
                payloads[j++] = fan2child->p1;
                keys[j] = fan2child->k2;
                payloads[j++] = fan2child->p2;
                delete fan2child;
            }
#endif
        }
        keys[j] = k1;
        payloads[j++] = p1;
        keys[j] = k2;
        payloads[j++] = p2;
        for(int i = pos + 1; i < fanout; ++i) {
            keyPayload &kp = kp_data[i];
            if (kp.key >= 0) {
                keys[j] = kp.key;
                payloads[j++] = kp.payload;
            } else if (kp.key == -1) {
                diliNode *child = kp.child;
                child->collect_and_clear(keys+j, payloads+j);
                j += child->num_nonempty;
                delete child;
            }
#ifndef ALLOW_FAN2_NODE
            else if (kp.key == -2) {
                fan2Leaf *fan2child = kp.fan2child;
                keys[j] = fan2child->k1;
                payloads[j++] = fan2child->p1;
                keys[j] = fan2child->k2;
                payloads[j++] = fan2child->p2;
                delete fan2child;
            }
#endif
        }
        delete[] kp_data;
        kp_data = NULL;
        assert(j == num_nonempty);
    }

    void collect_all_keys(long *keys) {
        assert(b >= 0);
        int j = 0;
        for(int i = 0; i < fanout; ++i) {
            keyPayload &kp = kp_data[i];
            if (kp.key >= 0) {
                keys[j++] = kp.key;
            } else if (kp.key == -1) {
                diliNode *child = kp.child;
                child->collect_all_keys(keys+j);
                j += child->num_nonempty;
            } else if (kp.key == -2) {
                fan2Leaf *fan2child = kp.fan2child;
                keys[j++] = fan2child->k1;
                keys[j++] = fan2child->k2;
            }
        }
        if (j != num_nonempty) {
            cout << "j = " << j << ", num_nonempty = " << num_nonempty << ", is_internal = " << is_internal() << endl;
            for(int i = 0; i < num_nonempty; ++i) {
                keyPayload &kp = kp_data[i];
                if (kp.key >= 0) {
                    cout << "i = " << i << ", kp.key = " << kp.key << endl;
                } else if (kp.key == -1) {
                    diliNode *child = kp.child;
                    child->collect_all_keys(keys+j);
                    cout << "i = " << i << ", child.num_nonempty = " << child->num_nonempty << endl;
                    j += child->num_nonempty;
                } else if (kp.key == -2) {
                    fan2Leaf *fan2child = kp.fan2child;
                    cout << "i = " << i << ", fan2child.num_nonempty = 2" << endl;
                    keys[j++] = fan2child->k1;
                    keys[j++] = fan2child->k2;
                }
            }
        }
        assert(j == num_nonempty);
    }


    void debug_help() {}

};


#endif //DILI_DILINODE_H
