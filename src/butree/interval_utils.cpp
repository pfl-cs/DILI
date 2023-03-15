#include "intervalInstance.h"
#include "interval_utils.h"
#include "../global/global.h"
#include "../utils/data_utils.h"
#include "../utils/file_utils.h"
#include <set>
#include <queue>
#include <algorithm>
#include <functional>
#include <utility>
#include <cassert>
#include <cmath>
#include <iostream>
#include <vector>
using namespace std;



void check_intervals(interval *i_ptr) {
    int count = 1;
    while (i_ptr) {
        interval *lSib = i_ptr->lSib;
        if (lSib) {
            assert(lSib->rSib == i_ptr);
            assert(lSib->start_idx + lSib->fanout == i_ptr->start_idx);
            assert(lSib->check_lr_delta_x());
        }

        i_ptr = i_ptr->rSib;
        ++count;
    }
//    cout << "****count = " << count << endl;
}

double cal_loss(long N, long n_intervals, double total_linear_loss, bool print=false) {
    double avg_linear_loss = sqrt(total_linear_loss / N);
    avg_linear_loss = 1 + ((avg_linear_loss > 1) ? (2 * log(avg_linear_loss) / log(2.0)) : 0);
    assert(avg_linear_loss >= 1);

    double base = 1.0 * N / n_intervals;
    double n_stages = log(N) / log(base);
//    double loss = ((4 + R1) * n_stages - 1) + (1 + R2) * avg_linear_loss;


//    return loss;
//
    double s = 0;
    double x = n_stages;
    int h = 0;
    double item = 0;
    const double r = 0.1;
    while (true) {
        item = pow(r, h);
        ++h;
        if (x <= 1) {
            s += item * x;
            break;
        }
        x -= 1;
        s += item;
    }


    double loss = ((2 + R1) * n_stages - 1) + s * (1 + R2) * avg_linear_loss;
    if (print) {
        cout << "n_stages = " << n_stages << ", linear_loss = " << avg_linear_loss << ", xl_loss = "
             << (1 + R2) * avg_linear_loss << ", other_loss = " << (2 + R1) * n_stages - 1 << ", s = " << s << endl;
    }
    return loss;
}

//longVec &complete_borders,
double partition(const keyArray &X, const doubleArray &probs, long N, int h, longVec &best_borders, long min_border_size, longVec &complete_borders, doubleVec &complete_losses, int interval_type){
    interval::data = X.get();
    interval::probs = probs.get();

    double single_prob = 1.0; // 1.0 / N;

    set<long> border_set;
    keyType last_ubd = X[0] - 1;

    interval *lSib = intervalInstance::newInstance(interval_type);
    lSib->fanout = 2;
    if (probs) {
        lSib->prob_sum = probs[0] + probs[1];
    } else {
        lSib->prob_sum = single_prob * 2;
    }
    lSib->start_idx = 0;
    lSib->end_idx = 2;
    lSib->lbd = last_ubd;
    last_ubd = lSib->ubd = X[1] + 1;
    lSib->init_merge_info();

    interval *leftest_interval = lSib;
    long s_j = 1;

    size_t last_size = 0;
    interval_set cs;
    long data_idx = 0;

    cout << "step-1 starts." << endl;
    auto start = chrono::system_clock::now();
    for (long i = 2; i < N; i += 2) {
        keyType x = X[i];
//        long x0 = x - 1;
//        long x2 = x + 1;
        keyType x2 = X[i + 1] + 1;
        int fan = 2;
        if (i + 3 >= N) {
            fan = N - i;
            x2 = X[i + fan - 1] + 1;
        }

//        if (x > last_ubd) {
//            interval *i_ptr = intervalInstance::newInstance(interval_type);
//            i_ptr->fanout = 0;
//            i_ptr->prob_sum = 0;
//            i_ptr->start_idx = i;
//            i_ptr->end_idx = i;
//            i_ptr->lbd = last_ubd;
//            last_ubd = i_ptr->ubd = x;
//            i_ptr->init_merge_info();
//
//            lSib->rSib = i_ptr;
//            i_ptr->lSib = lSib;
//            lSib->cal_merge_info();
//            cs.insert(lSib);
//            lSib = i_ptr;
//            last_size = cs.size();
//            ++s_j;
//        }


        interval *i_ptr = intervalInstance::newInstance(interval_type);
        i_ptr->fanout = fan;
        if (probs) {
            lSib->prob_sum = probs[i] + probs[i + 1];
            if (fan > 2) {
                lSib->prob_sum += probs[i + 1];
            }
        } else {
            lSib->prob_sum = single_prob * fan;
        }
        i_ptr->start_idx = i;
        i_ptr->end_idx = i + fan;
//        i_ptr->lbd = std::max<long>(last_ubd, x0);
        i_ptr->lbd = x;
        last_ubd = i_ptr->ubd = x2;

        lSib->rSib = i_ptr;
        i_ptr->lSib = lSib;

        i_ptr->init_merge_info();
        lSib->cal_merge_info(h);
        cs.insert(lSib);
        if (cs.size() != (last_size + 1)) {
            cout << "cs.size() = " << cs.size() << endl;
        }
        assert(cs.size() == last_size + 1);
        last_size = cs.size();
        ++s_j;
        lSib = i_ptr;
        if (fan > 2) {
            break;
        }
    }

    auto stop = chrono::system_clock::now();
    auto duration = chrono::duration_cast<chrono::seconds>(stop - start);
    double response_time = static_cast<double>(duration.count());
    cout << "step-1 takes " << response_time << " seconds." << endl;
//    border_set_check(cs);
//    check_intervals(leftest_interval);

//    interval *tmp_ptr = leftest_interval;
//    interval *tmp_l_ptr = leftest_interval;
//    while (tmp_ptr) {
//        if (tmp_ptr->start_idx == 20) {
//            break;
//        }
//        tmp_ptr = tmp_ptr->rSib;
//    }


    long max_s_j = N / 32;
    long least_num_merge = s_j - max_s_j;

//    cout << "step-2 finished. s_j = " << s_j << ", least_num_merge = " << least_num_merge << ", cs.size() = " << cs.size() << endl;

    long n_merges = 0;
    last_size = cs.size();

    cout << "step-2 starts." << endl;
    start = chrono::system_clock::now();
    while (n_merges < least_num_merge) {
        bool merge_flag = true;
        interval *i_ptr = *(cs.begin());
        cs.erase(i_ptr);

//        if (n_merges >= 4129) {
//            cout << "******i_ptr->start_idx = " << i_ptr->start_idx << ", i_ptr->delta_x = " << i_ptr->lr->get_delta_x()
//                    << ", i_ptr->merge_lr->delta_x = " << static_cast<buInterval*>(i_ptr)->merge_lr->get_delta_x() << endl;
//        }

        interval *original_rSib = i_ptr->rSib;
        assert(original_rSib);
        cs.erase(original_rSib);
        merge_flag = i_ptr->merge_with_rSib(h);

        if (merge_flag) {
            original_rSib->rSib = NULL;
            original_rSib->lSib = NULL;
            original_rSib->setInvalid();
            original_rSib->free_data();
            delete original_rSib;

            if (i_ptr->rSib) {
                i_ptr->cal_merge_info(h);
                cs.insert(i_ptr);
            }

            interval *lSib = i_ptr->lSib;

            if (lSib) {
                assert(lSib->rSib == i_ptr);
                cs.erase(lSib);

                lSib->cal_merge_info(h);
                cs.insert(lSib);
            }
            ++n_merges;
//            if (n_merges % 1000000 == 0) {
////            check_intervals(leftest_interval);
//                cout << "step-3, n_merges = " << n_merges << endl;
//            }

            assert(last_size == cs.size() + 1);
            last_size = cs.size();
        }  else {
            cs.insert(original_rSib);
            cs.insert(i_ptr);
        }
    }
    stop = chrono::system_clock::now();
    duration = chrono::duration_cast<chrono::seconds>(stop - start);
    response_time = static_cast<double>(duration.count());
    cout << "step-2 takes " << response_time << " seconds." << endl;
//    check_intervals(leftest_interval);


//    cout << "cs.size() = " << cs.size() << endl;
    interval *i_ptr = leftest_interval;

    double total_linear_loss = 0;

    cout << "step-3 starts." << endl;
    start = chrono::system_clock::now();
    long sum_fan = 0;
    while (i_ptr) {
        border_set.insert(X[i_ptr->start_idx]); // border_set.insert(i_ptr->lbd);
//        assert(!(i_ptr->lr));
        i_ptr->linear_loss = i_ptr->init_lr();
        total_linear_loss += i_ptr->linear_loss;
        sum_fan += i_ptr->fanout;
        i_ptr = i_ptr->rSib;

    }
    stop = chrono::system_clock::now();
    duration = chrono::duration_cast<chrono::seconds>(stop - start);
    response_time = static_cast<double>(duration.count());
    cout << "step-3 takes " << response_time << " seconds." << endl;
//    border_set.insert(X[N-1] + 1);

//    cout << "step-4 finished. border_set.size = " << border_set.size() << ", sum_fan = " << sum_fan << endl;

    check_intervals(leftest_interval);

    double best_loss = cal_loss(N, border_set.size() - 1, total_linear_loss, true);

    n_merges = 0;
    longVec deleted_borders;
    doubleVec deleted_losses;

    long best_idx = 0;
//    pair_max_heap pmh;

    long min_size = N / 10000;
    min_size = std::max<long>(min_size, 100);
    if (min_border_size > 0) {
        min_size = min_border_size;
    }

    cout << "step-4 starts." << endl;
    start = chrono::system_clock::now();
    while (border_set.size() > min_size) {
        interval *i_ptr = *(cs.begin());
        cs.erase(i_ptr);

        interval *original_rSib = i_ptr->rSib;
        cs.erase(original_rSib);

        long border_to_delete = original_rSib->lbd;
        double _total_loss = total_linear_loss - i_ptr->linear_loss - original_rSib->linear_loss;
        bool merge_flag = i_ptr->merge_with_rSib(h, true);
        if (merge_flag) {
            border_set.erase(border_to_delete);
            deleted_borders.push_back(border_to_delete);
//            total_linear_loss -= i_ptr->linear_loss + original_rSib->linear_loss;
            total_linear_loss = _total_loss + i_ptr->linear_loss;

            original_rSib->rSib = NULL;
            original_rSib->setInvalid();
            original_rSib->free_data();
            delete original_rSib;

            if (i_ptr->rSib) {
                i_ptr->cal_merge_info(h);
                cs.insert(i_ptr);
            }

            interval *lSib = i_ptr->lSib;
            if (lSib) {
                cs.erase(lSib);

                lSib->cal_merge_info(h);
                cs.insert(lSib);
            }


            ++n_merges;
            bool print = false;
            if (n_merges % 10000 == 0) {
                print = true;
                check_intervals(leftest_interval);
            }
            double curr_loss = cal_loss(N, border_set.size() - 1, total_linear_loss, print);
            deleted_losses.push_back(curr_loss);


//        if (pmh.size() < 5) {
//            pmh.push(make_pair(curr_loss, n_merges));
//        } else {
//            dl_pair max_p = pmh.top();
//            if (curr_loss < max_p.first) {
//                pmh.pop();
//                pmh.push(make_pair(curr_loss, n_merges));
//            }
//        }


            if (best_loss > curr_loss) {
                best_idx = n_merges;
                best_loss = curr_loss;
            }

//            if (n_merges % 10000 == 0) {
//                cout << "step-5, n_merges = " << n_merges << ", n_intervals = " << border_set.size() - 1 << ", loss = "
//                     << curr_loss << ", best_loss = " << best_loss << endl;
//            }
        } else {
            cs.insert(i_ptr);
            cs.insert(original_rSib);
        }
    }
    stop = chrono::system_clock::now();
    duration = chrono::duration_cast<chrono::seconds>(stop - start);
    response_time = static_cast<double>(duration.count());
    cout << "step-4 takes " << response_time << " seconds." << endl;


    // complete_borders should not be sorted.
    complete_borders.clear();
    complete_borders.insert(complete_borders.end(), border_set.begin(), border_set.end());
    complete_borders.insert(complete_borders.end(), deleted_borders.rbegin(), deleted_borders.rend());
    complete_losses.insert(complete_losses.end(), deleted_losses.rbegin(), deleted_losses.rend());

    border_set.insert(deleted_borders.begin() + best_idx, deleted_borders.end());
//    cout << "best_border_set.size() = " << border_set.size() << ", best_idx = " << best_idx << endl;


    best_borders.insert(best_borders.end(), border_set.begin(), border_set.end());
    std::sort(best_borders.begin(), best_borders.end());

//    int x = 5;
//    while (!pmh.empty()) {
//        dl_pair p = pmh.top();
//        cout << "top-" << x << "-idx = " << p.second << endl;
//        pmh.pop();
//        --x;
//    }
    i_ptr = leftest_interval;
    while (i_ptr) {
        interval *tmp = i_ptr->rSib;
        delete i_ptr;
        i_ptr = tmp;
    }

    return best_loss;
}

void get_complete_partition_borders_for_lowest_layer(const keyArray &X, long N, int h, long min_fanout, long max_fanout, longVec &complete_borders, doubleVec &complete_avg_rmses, int interval_type){
    interval::data = X.get();

    double single_prob = 1.0; // 1.0 / N;

    set<long> border_set;
    keyType last_ubd = X[0] - 1;

    interval *lSib = intervalInstance::newInstance(interval_type);
    const int init_fan = 4;

    lSib->fanout = init_fan;
    lSib->prob_sum = single_prob * init_fan;

    lSib->start_idx = 0;
    lSib->end_idx = init_fan;
    lSib->lbd = last_ubd;
    last_ubd = lSib->ubd = X[1] + 1;
    lSib->init_merge_info();

    interval *leftest_interval = lSib;
    long s_j = 1;

    size_t last_size = 0;
    interval_set cs;
    long data_idx = 0;

    cout << "step-1 starts." << endl;
    auto start = chrono::system_clock::now();
    for (long i = init_fan; i < N; i += init_fan) {
        keyType x = X[i];
        int fan = init_fan;
        if (i + init_fan + 1 >= N) {
            fan = N - i;
        }
        keyType x2 = X[i + fan - 1] + 1;

        interval *i_ptr = intervalInstance::newInstance(interval_type);
        i_ptr->fanout = fan;

        lSib->prob_sum = single_prob * fan;

        i_ptr->start_idx = i;
        i_ptr->end_idx = i + fan;
//        i_ptr->lbd = std::max<long>(last_ubd, x0);
        i_ptr->lbd = x;
        last_ubd = i_ptr->ubd = x2;

        lSib->rSib = i_ptr;
        i_ptr->lSib = lSib;

        i_ptr->init_merge_info();
        lSib->cal_merge_info(h);

        cs.insert(lSib);
//        if (cs.size() != (last_size + 1)) {
//            cout << "cs.size() = " << cs.size() << ", last_size = " << last_size << endl;
//        }
//        assert(cs.size() == last_size + 1);
        last_size = cs.size();
        ++s_j;
        lSib = i_ptr;
        if (fan > init_fan) {
            break;
        }
    }
    if (s_j != cs.size() + 1) {
        cout << "s_j = " << s_j << ", cs.size() = " << cs.size() << endl;
    }

    assert (s_j == cs.size() + 1);

    auto stop = chrono::system_clock::now();
    auto duration = chrono::duration_cast<chrono::seconds>(stop - start);
    double response_time = static_cast<double>(duration.count());
    cout << "step-1 takes " << response_time << " seconds." << endl;


    assert(s_j == cs.size() + 1);

    long max_s_j = N / min_fanout;
    long least_num_merge = s_j - max_s_j;

//    cout << "step-2 finished. s_j = " << s_j << ", least_num_merge = " << least_num_merge << ", cs.size() = " << cs.size() << endl;

    long n_merges = 0;
    last_size = cs.size();

    cout << "step-2 starts." << endl;
    start = chrono::system_clock::now();
    while (n_merges < least_num_merge) {
        interval *i_ptr = *(cs.begin());
        cs.erase(i_ptr);

        interval *original_rSib = i_ptr->rSib;
        assert(original_rSib);
        cs.erase(original_rSib);

        bool merge_flag = i_ptr->merge_with_rSib(h);

        if (merge_flag) {
            original_rSib->rSib = NULL;
            original_rSib->lSib = NULL;
            original_rSib->setInvalid();
            original_rSib->free_data();
            delete original_rSib;

            if (i_ptr->rSib) {
                i_ptr->cal_merge_info(h);
                cs.insert(i_ptr);
            }


            interval *lSib = i_ptr->lSib;

            if (lSib) {
                assert(lSib->rSib == i_ptr);
                cs.erase(lSib);

                lSib->cal_merge_info(h);
                cs.insert(lSib);
            }
            ++n_merges;
//            if (n_merges % 100000 == 0) {
//                cout << "step-3, n_merges = " << n_merges << endl;
//            }

            assert(last_size == cs.size() + 1);
            last_size = cs.size();
        } else {
            cs.insert(original_rSib);
            cs.insert(i_ptr);
        }
    }
    stop = chrono::system_clock::now();
    duration = chrono::duration_cast<chrono::seconds>(stop - start);
    response_time = static_cast<double>(duration.count());
    cout << "step-2 takes " << response_time << " seconds." << endl;

//    cout << "cs.size() = " << cs.size() << endl;
    if (max_s_j != cs.size() + 1) {
        cout << "max_s_j = " << max_s_j << ", cs.size = " << cs.size() << endl;
    }
    assert(max_s_j == cs.size() + 1);
    interval *i_ptr = leftest_interval;

    cout << "step-3 starts." << endl;
    start = chrono::system_clock::now();

    long sum_fan = 0;
    double total_linear_loss = 0;
    while (i_ptr) {
        border_set.insert(X[i_ptr->start_idx]);
//        sum_fan += i_ptr->fanout;
        i_ptr->linear_loss = i_ptr->init_lr();
        total_linear_loss += i_ptr->linear_loss;
        sum_fan += i_ptr->fanout;
        i_ptr = i_ptr->rSib;
    }
    stop = chrono::system_clock::now();
    duration = chrono::duration_cast<chrono::seconds>(stop - start);
    response_time = static_cast<double>(duration.count());
    cout << "step-3 takes " << response_time << " seconds." << endl;

    if (max_s_j != border_set.size()) {
        cout << "max_s_j = " << max_s_j << ", border_set.size = " << border_set.size() << endl;
    }
    assert(border_set.size() == max_s_j);
//    border_set.insert(X[N-1] + 1);
//    cout << "step-4 finished. border_set.size = " << border_set.size() << ", sum_fan = " << sum_fan << endl;

    check_intervals(leftest_interval);

    n_merges = 0;
    longVec deleted_borders;
    long min_size = N / max_fanout; // assert(max_fanout <= fanThreshold / 2);
//    min_size = std::max<long>(min_size, 100);
    if (min_size <= 0) {
        min_size = 2;
    }

    cout << "step-4 starts." << endl;
    start = chrono::system_clock::now();

    doubleVec all_rmses;
    while (border_set.size() > min_size) {
        interval *i_ptr = *(cs.begin());
        cs.erase(i_ptr);

        interval *original_rSib = i_ptr->rSib;
        cs.erase(original_rSib);

        long border_to_delete = original_rSib->lbd;

        double _total_loss = total_linear_loss - i_ptr->linear_loss - original_rSib->linear_loss;
        bool merge_flag = i_ptr->merge_with_rSib(h, true);

        if (merge_flag) {
            border_set.erase(border_to_delete);
            deleted_borders.push_back(border_to_delete);
            total_linear_loss = _total_loss + i_ptr->linear_loss;

            double rmse = sqrt(total_linear_loss / N);
            all_rmses.push_back(rmse);

            original_rSib->rSib = NULL;
            original_rSib->setInvalid();
            original_rSib->free_data();
            delete original_rSib;

            if (i_ptr->rSib) {
                i_ptr->cal_merge_info(h);
                cs.insert(i_ptr);
            }

            interval *lSib = i_ptr->lSib;
            if (lSib) {
                cs.erase(lSib);
                lSib->cal_merge_info(h);
                cs.insert(lSib);
            }

            ++n_merges;
//            if (n_merges % 100000 == 0) {
//                cout << "step-5, n_merges = " << n_merges << ", n_intervals = " << border_set.size() - 1 << endl;
//            }
        } else {
            cs.insert(original_rSib);
            cs.insert(i_ptr);
        }
    }
    stop = chrono::system_clock::now();
    duration = chrono::duration_cast<chrono::seconds>(stop - start);
    response_time = static_cast<double>(duration.count());
    cout << "step-4 takes " << response_time << " seconds." << endl;

    assert(min_size == border_set.size());
    assert(min_size + deleted_borders.size() == max_s_j);

    // complete_borders should not be sorted.
    complete_borders.clear();
    complete_borders.insert(complete_borders.end(), border_set.begin(), border_set.end());
    complete_borders.insert(complete_borders.end(), deleted_borders.rbegin(), deleted_borders.rend());
    complete_borders.push_back(min_size);

    complete_avg_rmses.clear();
    for (int i = 0; i < min_size; ++i) {
        complete_avg_rmses.push_back(1e50);
    }
    complete_avg_rmses.insert(complete_avg_rmses.end(), all_rmses.rbegin(), all_rmses.rend());

    i_ptr = leftest_interval;
    while (i_ptr) {
        interval *tmp = i_ptr->rSib;
        delete i_ptr;
        i_ptr = tmp;
    }
}

void get_complete_partition_borders_w_sampling(const keyArray &X, long N, int h, long min_fanout, long max_fanout, longVec &complete_borders, doubleVec &complete_avg_rmses, int interval_type){
    interval::data = X.get();

    double single_prob = 1.0; // 1.0 / N;

    set<long> border_set;
    keyType last_ubd = X[0] - 1;

    interval *lSib = intervalInstance::newInstance(interval_type);
    const int init_fan = 6;

    lSib->fanout = init_fan;
    lSib->prob_sum = single_prob * init_fan;

    lSib->start_idx = 0;
    lSib->end_idx = init_fan;
    lSib->lbd = last_ubd;
    last_ubd = lSib->ubd = X[1] + 1;
    lSib->init_merge_info_w_sampling();

    interval *leftest_interval = lSib;
    long s_j = 1;

    size_t last_size = 0;
    interval_set cs;
    long data_idx = 0;

    cout << "step-1 starts." << endl;
    auto start = chrono::system_clock::now();
    for (long i = init_fan; i < N; i += init_fan) {
        keyType x = X[i];
        int fan = init_fan;
        if (i + init_fan + 1 >= N) {
            fan = N - i;
        }
        keyType x2 = X[i + fan - 1] + 1;

        interval *i_ptr = intervalInstance::newInstance(interval_type);
        i_ptr->fanout = fan;

        lSib->prob_sum = single_prob * fan;

        i_ptr->start_idx = i;
        i_ptr->end_idx = i + fan;
//        i_ptr->lbd = std::max<long>(last_ubd, x0);
        i_ptr->lbd = x;
        last_ubd = i_ptr->ubd = x2;

        lSib->rSib = i_ptr;
        i_ptr->lSib = lSib;

        i_ptr->init_merge_info_w_sampling();
        lSib->cal_merge_info_w_sampling(h);

        cs.insert(lSib);
//        if (cs.size() != (last_size + 1)) {
//            cout << "cs.size() = " << cs.size() << ", last_size = " << last_size << endl;
//        }
//        assert(cs.size() == last_size + 1);
        last_size = cs.size();
        ++s_j;
        lSib = i_ptr;
        if (fan > init_fan) {
            break;
        }
    }
    if (s_j != cs.size() + 1) {
        cout << "s_j = " << s_j << ", cs.size() = " << cs.size() << endl;
    }

    assert (s_j == cs.size() + 1);

    auto stop = chrono::system_clock::now();
    auto duration = chrono::duration_cast<chrono::seconds>(stop - start);
    double response_time = static_cast<double>(duration.count());
    cout << "step-1 takes " << response_time << " seconds." << endl;


    assert(s_j == cs.size() + 1);

    long max_s_j = N / min_fanout;
    long least_num_merge = s_j - max_s_j;

//    cout << "step-2 finished. s_j = " << s_j << ", least_num_merge = " << least_num_merge << ", cs.size() = " << cs.size() << endl;

    long n_merges = 0;
    last_size = cs.size();

    cout << "step-2 starts." << endl;
    start = chrono::system_clock::now();
    while (n_merges < least_num_merge) {
        interval *i_ptr = *(cs.begin());
        cs.erase(i_ptr);

        interval *original_rSib = i_ptr->rSib;
        assert(original_rSib);
        cs.erase(original_rSib);

        bool merge_flag = i_ptr->merge_with_rSib(h);

        if (merge_flag) {
            original_rSib->rSib = NULL;
            original_rSib->lSib = NULL;
            original_rSib->setInvalid();
            original_rSib->free_data();
            delete original_rSib;

            if (i_ptr->rSib) {
                i_ptr->cal_merge_info_w_sampling(h);
                cs.insert(i_ptr);
            }


            interval *lSib = i_ptr->lSib;

            if (lSib) {
                assert(lSib->rSib == i_ptr);
                cs.erase(lSib);

                lSib->cal_merge_info_w_sampling(h);
                cs.insert(lSib);
            }
            ++n_merges;
//            if (n_merges % 100000 == 0) {
//                cout << "step-3, n_merges = " << n_merges << endl;
//            }

            assert(last_size == cs.size() + 1);
            last_size = cs.size();
        } else {
            cs.insert(original_rSib);
            cs.insert(i_ptr);
        }
    }
    stop = chrono::system_clock::now();
    duration = chrono::duration_cast<chrono::seconds>(stop - start);
    response_time = static_cast<double>(duration.count());
    cout << "step-2 takes " << response_time << " seconds." << endl;

//    cout << "cs.size() = " << cs.size() << endl;
    if (max_s_j != cs.size() + 1) {
        cout << "max_s_j = " << max_s_j << ", cs.size = " << cs.size() << endl;
    }
    assert(max_s_j == cs.size() + 1);
    interval *i_ptr = leftest_interval;

    cout << "step-3 starts." << endl;
    start = chrono::system_clock::now();

    long sum_fan = 0;
    double total_linear_loss = 0;
    while (i_ptr) {
        border_set.insert(X[i_ptr->start_idx]);
//        sum_fan += i_ptr->fanout;
        i_ptr->linear_loss = i_ptr->init_lr_w_sampling();
        total_linear_loss += i_ptr->linear_loss;
        sum_fan += i_ptr->fanout;
        i_ptr = i_ptr->rSib;
    }
    stop = chrono::system_clock::now();
    duration = chrono::duration_cast<chrono::seconds>(stop - start);
    response_time = static_cast<double>(duration.count());
    cout << "step-3 takes " << response_time << " seconds." << endl;

    if (max_s_j != border_set.size()) {
        cout << "max_s_j = " << max_s_j << ", border_set.size = " << border_set.size() << endl;
    }
    assert(border_set.size() == max_s_j);
//    border_set.insert(X[N-1] + 1);
//    cout << "step-4 finished. border_set.size = " << border_set.size() << ", sum_fan = " << sum_fan << endl;

    check_intervals(leftest_interval);

    n_merges = 0;
    longVec deleted_borders;
    long min_size = N / max_fanout; // assert(max_fanout <= fanThreshold / 2);
//    min_size = std::max<long>(min_size, 100);
    if (min_size <= 0) {
        min_size = 2;
    }

    cout << "step-4 starts." << endl;
    start = chrono::system_clock::now();

    doubleVec all_rmses;
    while (border_set.size() > min_size) {
        interval *i_ptr = *(cs.begin());
        cs.erase(i_ptr);

        interval *original_rSib = i_ptr->rSib;
        cs.erase(original_rSib);

        long border_to_delete = original_rSib->lbd;

        double _total_loss = total_linear_loss - i_ptr->linear_loss - original_rSib->linear_loss;
        bool merge_flag = i_ptr->merge_with_rSib(h, true);

        if (merge_flag) {
            border_set.erase(border_to_delete);
            deleted_borders.push_back(border_to_delete);
            total_linear_loss = _total_loss + i_ptr->linear_loss;

            double rmse = sqrt(total_linear_loss / N);
            all_rmses.push_back(rmse);

            original_rSib->rSib = NULL;
            original_rSib->setInvalid();
            original_rSib->free_data();
            delete original_rSib;

            if (i_ptr->rSib) {
                i_ptr->cal_merge_info_w_sampling(h);
                cs.insert(i_ptr);
            }

            interval *lSib = i_ptr->lSib;
            if (lSib) {
                cs.erase(lSib);
                lSib->cal_merge_info_w_sampling(h);
                cs.insert(lSib);
            }

            ++n_merges;
//            if (n_merges % 100000 == 0) {
//                cout << "step-5, n_merges = " << n_merges << ", n_intervals = " << border_set.size() - 1 << endl;
//            }
        } else {
            cs.insert(original_rSib);
            cs.insert(i_ptr);
        }
    }
    stop = chrono::system_clock::now();
    duration = chrono::duration_cast<chrono::seconds>(stop - start);
    response_time = static_cast<double>(duration.count());
    cout << "step-4 takes " << response_time << " seconds." << endl;

    assert(min_size == border_set.size());
    assert(min_size + deleted_borders.size() == max_s_j);

    // complete_borders should not be sorted.
    complete_borders.clear();
    complete_borders.insert(complete_borders.end(), border_set.begin(), border_set.end());
    complete_borders.insert(complete_borders.end(), deleted_borders.rbegin(), deleted_borders.rend());
    complete_borders.push_back(min_size);

    complete_avg_rmses.clear();
    for (int i = 0; i < min_size; ++i) {
        complete_avg_rmses.push_back(1e50);
    }
    complete_avg_rmses.insert(complete_avg_rmses.end(), all_rmses.rbegin(), all_rmses.rend());

    i_ptr = leftest_interval;
    while (i_ptr) {
        interval *tmp = i_ptr->rSib;
        delete i_ptr;
        i_ptr = tmp;
    }
}


void get_complete_partition_borders(const keyArray &X, const doubleArray &probs, long N, int h, long min_fanout, long max_fanout, longVec &complete_borders, doubleVec &complete_avg_rmses, int interval_type){
    interval::data = X.get();
    interval::probs = probs.get();

    double single_prob = 1.0; // 1.0 / N;

    set<long> border_set;
    keyType last_ubd = X[0] - 1;

    interval *lSib = intervalInstance::newInstance(interval_type);
    lSib->fanout = 2;
    if (probs) {
        lSib->prob_sum = probs[0] + probs[1];
    } else {
        lSib->prob_sum = single_prob * 2;
    }
    lSib->start_idx = 0;
    lSib->end_idx = 2;
    lSib->lbd = last_ubd;
    last_ubd = lSib->ubd = X[1] + 1;
    lSib->init_merge_info();

    interval *leftest_interval = lSib;
    long s_j = 1;

    size_t last_size = 0;
    interval_set cs;
    long data_idx = 0;

    cout << "step-1 starts." << endl;
    auto start = chrono::system_clock::now();
    for (long i = 2; i < N; i += 2) {
        keyType x = X[i];
//        long x0 = x - 1;
//        long x2 = x + 1;
        keyType x2 = X[i + 1] + 1;
        int fan = 2;
        if (i + 3 >= N) {
            fan = N - i;
            x2 = X[i + fan - 1] + 1;
        }

        interval *i_ptr = intervalInstance::newInstance(interval_type);
        i_ptr->fanout = fan;
        if (probs) {
            lSib->prob_sum = probs[i] + probs[i + 1];
            if (fan > 2) {
                lSib->prob_sum += probs[i + 1];
            }
        } else {
            lSib->prob_sum = single_prob * fan;
        }
        i_ptr->start_idx = i;
        i_ptr->end_idx = i + fan;
//        i_ptr->lbd = std::max<long>(last_ubd, x0);
        i_ptr->lbd = x;
        last_ubd = i_ptr->ubd = x2;

        lSib->rSib = i_ptr;
        i_ptr->lSib = lSib;

        i_ptr->init_merge_info();
        lSib->cal_merge_info(h);
//        auto it = cs.find(lSib);
//        if (it != cs.end()) {
//            interval *tmp = *it;
//            cout << "i = " << i << ", tmp.merge_metric = " << tmp->merge_metric << ", tmp.lbd = " << tmp->lbd << ", lSib.merge_metric = " << lSib->merge_metric << ", lSib.lbd = " << lSib->lbd << endl;
//        }
        cs.insert(lSib);
        if (cs.size() != (last_size + 1)) {
            cout << "cs.size() = " << cs.size() << ", last_size = " << last_size << endl;
        }
        assert(cs.size() == last_size + 1);
        last_size = cs.size();
        ++s_j;
        lSib = i_ptr;
        if (fan > 2) {
            break;
        }
    }
    if (s_j != cs.size() + 1) {
        cout << "s_j = " << s_j << ", cs.size() = " << cs.size() << endl;
    }

    auto stop = chrono::system_clock::now();
    auto duration = chrono::duration_cast<chrono::seconds>(stop - start);
    double response_time = static_cast<double>(duration.count());
    cout << "step-1 takes " << response_time << " seconds." << endl;


    assert(s_j == cs.size() + 1);

    long max_s_j = N / min_fanout;
    long least_num_merge = s_j - max_s_j;

//    cout << "step-2 finished. s_j = " << s_j << ", least_num_merge = " << least_num_merge << ", cs.size() = " << cs.size() << endl;

    long n_merges = 0;
    last_size = cs.size();

    cout << "step-2 starts." << endl;
    start = chrono::system_clock::now();
    while (n_merges < least_num_merge) {
        interval *i_ptr = *(cs.begin());
        cs.erase(i_ptr);

        interval *original_rSib = i_ptr->rSib;
        assert(original_rSib);
        cs.erase(original_rSib);

        bool merge_flag = i_ptr->merge_with_rSib(h);

        if (merge_flag) {
            original_rSib->rSib = NULL;
            original_rSib->lSib = NULL;
            original_rSib->setInvalid();
            original_rSib->free_data();
            delete original_rSib;

            if (i_ptr->rSib) {
                i_ptr->cal_merge_info(h);
                cs.insert(i_ptr);
            }


            interval *lSib = i_ptr->lSib;

            if (lSib) {
                assert(lSib->rSib == i_ptr);
                cs.erase(lSib);

                lSib->cal_merge_info(h);
                cs.insert(lSib);
            }
            ++n_merges;
//            if (n_merges % 100000 == 0) {
//                cout << "step-3, n_merges = " << n_merges << endl;
//            }

            assert(last_size == cs.size() + 1);
            last_size = cs.size();
        } else {
            cs.insert(original_rSib);
            cs.insert(i_ptr);
        }
    }
    stop = chrono::system_clock::now();
    duration = chrono::duration_cast<chrono::seconds>(stop - start);
    response_time = static_cast<double>(duration.count());
    cout << "step-2 takes " << response_time << " seconds." << endl;

//    cout << "cs.size() = " << cs.size() << endl;
    if (max_s_j != cs.size() + 1) {
        cout << "max_s_j = " << max_s_j << ", cs.size = " << cs.size() << endl;
    }
    assert(max_s_j == cs.size() + 1);
    interval *i_ptr = leftest_interval;

    cout << "step-3 starts." << endl;
    start = chrono::system_clock::now();

    long sum_fan = 0;
    double total_linear_loss = 0;
    while (i_ptr) {
        border_set.insert(X[i_ptr->start_idx]);
//        sum_fan += i_ptr->fanout;
        i_ptr->linear_loss = i_ptr->init_lr();
        total_linear_loss += i_ptr->linear_loss;
        sum_fan += i_ptr->fanout;
        i_ptr = i_ptr->rSib;
    }
    stop = chrono::system_clock::now();
    duration = chrono::duration_cast<chrono::seconds>(stop - start);
    response_time = static_cast<double>(duration.count());
    cout << "step-3 takes " << response_time << " seconds." << endl;

    if (max_s_j != border_set.size()) {
        cout << "max_s_j = " << max_s_j << ", border_set.size = " << border_set.size() << endl;
    }
    assert(border_set.size() == max_s_j);
//    border_set.insert(X[N-1] + 1);
//    cout << "step-4 finished. border_set.size = " << border_set.size() << ", sum_fan = " << sum_fan << endl;

    check_intervals(leftest_interval);

    n_merges = 0;
    longVec deleted_borders;
    long min_size = N / max_fanout; // assert(max_fanout <= fanThreshold / 2);
//    min_size = std::max<long>(min_size, 100);
    if (min_size <= 0) {
        min_size = 2;
    }

    cout << "step-4 starts." << endl;
    start = chrono::system_clock::now();

    doubleVec all_rmses;
    while (border_set.size() > min_size) {
        interval *i_ptr = *(cs.begin());
        cs.erase(i_ptr);

        interval *original_rSib = i_ptr->rSib;
        cs.erase(original_rSib);

        long border_to_delete = original_rSib->lbd;

        double _total_loss = total_linear_loss - i_ptr->linear_loss - original_rSib->linear_loss;
        bool merge_flag = i_ptr->merge_with_rSib(h, true);

        if (merge_flag) {
            border_set.erase(border_to_delete);
            deleted_borders.push_back(border_to_delete);
            total_linear_loss = _total_loss + i_ptr->linear_loss;

            double rmse = sqrt(total_linear_loss / N);
            all_rmses.push_back(rmse);

            original_rSib->rSib = NULL;
            original_rSib->setInvalid();
            original_rSib->free_data();
            delete original_rSib;

            if (i_ptr->rSib) {
                i_ptr->cal_merge_info(h);
                cs.insert(i_ptr);
            }

            interval *lSib = i_ptr->lSib;
            if (lSib) {
                cs.erase(lSib);
                lSib->cal_merge_info(h);
                cs.insert(lSib);
            }

            ++n_merges;
//            if (n_merges % 100000 == 0) {
//                cout << "step-5, n_merges = " << n_merges << ", n_intervals = " << border_set.size() - 1 << endl;
//            }
        } else {
            cs.insert(original_rSib);
            cs.insert(i_ptr);
        }
    }
    stop = chrono::system_clock::now();
    duration = chrono::duration_cast<chrono::seconds>(stop - start);
    response_time = static_cast<double>(duration.count());
    cout << "step-4 takes " << response_time << " seconds." << endl;

    assert(min_size == border_set.size());
    assert(min_size + deleted_borders.size() == max_s_j);

    // complete_borders should not be sorted.
    complete_borders.clear();
    complete_borders.insert(complete_borders.end(), border_set.begin(), border_set.end());
    complete_borders.insert(complete_borders.end(), deleted_borders.rbegin(), deleted_borders.rend());
    complete_borders.push_back(min_size);

    complete_avg_rmses.clear();
    for (int i = 0; i < min_size; ++i) {
        complete_avg_rmses.push_back(1e50);
    }
    complete_avg_rmses.insert(complete_avg_rmses.end(), all_rmses.rbegin(), all_rmses.rend());

    i_ptr = leftest_interval;
    while (i_ptr) {
        interval *tmp = i_ptr->rSib;
        delete i_ptr;
        i_ptr = tmp;
    }
}

void build_mirror(const keyArray &X, const doubleArray &probs, long N, l_matrix &mirror, const string &mirror_dir, int interval_type) {
    mirror.clear();
    restore_mirror(mirror_dir, mirror);

    if (mirror.size() > 0 && mirror[mirror.size()-1].size() == 1) {
        return;
    }
    int curr_n_nodes = -1;

    if (mirror.size() <= 0) {
        mirror.push_back(longVec());
        longVec h0_borders;
        doubleVec h0_losses;
        partition(X, probs, N, 0, mirror[0], N / 200, h0_borders, h0_losses, interval_type);
        string h0_path = mirror_dir + "/h0.dat";
        data_utils::save_vec_data(h0_path.c_str(), mirror[0]);
    }
    curr_n_nodes = mirror[mirror.size() - 1].size();
    cout << "leaf nodes have been created" << endl;


    int height = mirror.size();
    longVec root_info({X[0]});
    while (true) {
        // build root_cand

        longVec &highest_bounds = mirror[mirror.size()-1];
//        keyType *topX = new keyType [curr_n_nodes+1];
        keyArray topX = make_unique<keyType []>(curr_n_nodes+1);
        std::copy(highest_bounds.begin(), highest_bounds.end(), topX.get());
        topX[curr_n_nodes] = X[N-1] + 1;

        cout << "stage-" <<  (height+1) << " from bottom. curr_n_nodes = " << curr_n_nodes << endl;
        if (height > 1 && curr_n_nodes <= 200001) {
            mirror.push_back(root_info);
            assert(mirror.size() == height + 1);
            string hi_path = mirror_dir + "/h" + to_string(height) + ".dat";
            data_utils::save_vec_data(hi_path.c_str(), root_info);

            break;
        }

        linearRegressor root_lr;
        root_lr.init(topX.get(), 0, 0, curr_n_nodes);
        double root_loss = root_lr.cal_loss(topX.get(), curr_n_nodes);
        root_loss = sqrt(root_loss / curr_n_nodes);
        root_loss = (root_loss > 1) ? log(root_loss) / log(2.0) : 0;
//        root_loss = 1 + ((root_loss > 1) ? (2 * log(root_loss) / log(2.0)) : 0);
//        root_loss *= (1 + R2);
        cout << "root_loss = " << root_loss << endl;

        mirror.push_back(longVec());
        longVec hi_borders;
        doubleVec hi_losses;
        double int_loss = partition(topX, nullptr, curr_n_nodes, height, mirror[mirror.size()-1], -1, hi_borders, hi_losses, interval_type);

        cout << "int_loss = " << int_loss << endl;
        // check if one more stage needs to be added.
        ++height;
        if (height > 2 && root_loss < int_loss) {
            mirror.pop_back();
            mirror.push_back(root_info);
            assert(mirror.size() == height);

            string hi_path = mirror_dir + "/h" + to_string(height-1) + ".dat";
            data_utils::save_vec_data(hi_path.c_str(), root_info);
            break;
        } else {
            longVec &lv = mirror[mirror.size()-1];
            curr_n_nodes = lv.size();
            assert(mirror.size() == height);
            string hi_path = mirror_dir + "/h" + to_string(height-1) + ".dat";
            data_utils::save_vec_data(hi_path.c_str(), mirror[mirror.size()-1]);
        }
    }
}


void build_mirror_from_given_layout(const keyArray &X, const doubleArray &probs, long N, l_matrix &mirror, const string &mirror_dir, const longVec &layout, int interval_type) {
    mirror.clear();
    restore_mirror(mirror_dir, mirror);

    if (mirror.size() > 0 && mirror[mirror.size()-1].size() == 1) {
        return;
    }


    if (mirror.size() <= 0) {
        longVec h0_borders;
        doubleVec h0_rmses;
        bool restore_status = restore_complete_borders(mirror_dir, 0, h0_borders);
        if (!restore_status) {
//            get_complete_partition_borders(X, probs, N, 0, 32, fanThreashold/2, h0_borders, h0_rmses, interval_type);
            get_complete_partition_borders(X, probs, N, 0, buMinFan, fanThreashold/2, h0_borders, h0_rmses, interval_type);
            string h0_borders_path = mirror_dir + "/h0_borders";
            data_utils::save_vec_data(h0_borders_path.c_str(), h0_borders);

            string h0_rmses_path = mirror_dir + "/h0_rmses";
            data_utils::save_vec_data(h0_rmses_path.c_str(), h0_rmses);
        }

        long n_nodes = layout[0];
        long h0_min_size = h0_borders.back();
        h0_borders.pop_back();
        cout << "h0_borsers.size() = " << h0_borders.size() << endl;

        mirror.push_back(longVec());

        assert(h0_min_size < n_nodes);
        mirror[0].insert(mirror[0].end(), h0_borders.begin(), h0_borders.begin() + n_nodes);
        sort(mirror[0].begin(), mirror[0].end());

//        string h0_path = mirror_dir + "/h0.dat";
//        data_utils::save_vec_data(h0_path.c_str(), mirror[0]);
    }
    int curr_n_nodes = mirror[mirror.size() - 1].size();
    cout << "leaf nodes have been created. n_nodes = " << curr_n_nodes << ", mirror.size = " << mirror.size() << endl;

    int H = layout.size();
    for (int i = 1; i < H - 1; ++i) {
        longVec &bounds = mirror[i-1];
        keyArray topX = make_unique<keyType []>(curr_n_nodes+1);
        std::copy(bounds.begin(), bounds.end(), topX.get());
        topX[curr_n_nodes] = X[N-1] + 1;

        longVec hi_borders;
        doubleVec hi_rmses;
        bool restore_status = restore_complete_borders(mirror_dir, i, hi_borders);
        if (!restore_status) {
            get_complete_partition_borders(topX, NULL, curr_n_nodes, i, 100, fanThreashold, hi_borders, hi_rmses, interval_type);
//            get_complete_partition_borders(topX, NULL, curr_n_nodes, i, buMinFan * 2, fanThreashold, hi_borders, hi_rmses, interval_type);
        }
        long hi_min_size = hi_borders.back();
        hi_borders.pop_back();
        mirror.push_back(longVec());
        cout << "hi_borsers.size() = " << hi_borders.size() << endl;

        long n_nodes = layout[i];

        assert(hi_min_size < n_nodes);
        mirror[i].insert(mirror[i].end(), hi_borders.begin(), hi_borders.begin() + n_nodes);
        sort(mirror[i].begin(), mirror[i].end());

//        string hi_path = mirror_dir + "/h" + to_string(i) + ".dat";
//        data_utils::save_vec_data(hi_path.c_str(), mirror[i]);
    }

    longVec root_info({X[0]});
    mirror.push_back(root_info);
}


double loss_est_complex(int h, long n_nodes, long next_n_nodes, long N_at_h0, doubleVec &h0_rmses, double _rho, double rmse_current_height, bool print=false) {

    double avg_fan = 1.0 * n_nodes / next_n_nodes;
    double _n_stages = log(n_nodes) / log(avg_fan);
    int n_stages = static_cast<int>(_n_stages);
    double residual = _n_stages - n_stages;
    long idx = static_cast<long>(N_at_h0 / avg_fan);
    double next_level_rmse = 1e50;
    if (idx < h0_rmses.size()) {
        next_level_rmse = h0_rmses[idx];
    }

    double s = 1;
    for (int i = 0; i <= n_stages; ++i) {
        s *= _rho;
    }
    double last_item = s * residual;
    s = (s - _rho) / (_rho - 1) + last_item;

    double rh = 1;
    for (int i = 0; i < h; ++i) {
        rh *= _rho;
    }
    assert(n_stages >= 1);


    assert(s >= 0);
    assert(next_level_rmse >= 0);
    const int n_pairs_each_cache_line = 4;

    double loss_curr_h_lin = rmse_current_height / n_pairs_each_cache_line;
    double loss_next_lin = s * next_level_rmse / n_pairs_each_cache_line;

    double loss_curr_h_bin = ((rmse_current_height > 2) ? (2 * log(rmse_current_height) / log(2.0)) : 2);
    double loss_nex_bin = s * ((next_level_rmse > 2) ? (2 * log(next_level_rmse) / log(2.0)) : 2);

//    double loss_curr_h = std::min<double>(loss_curr_h_lin, loss_curr_h_bin);
//    double loss_next = std::min<double>(loss_next_lin, loss_nex_bin);

    double loss_curr_h = loss_curr_h_bin;
    double loss_next = loss_nex_bin;

    double avg_linear_loss = (1 + R2) * rh * (loss_curr_h + loss_next);
//    double avg_linear_loss = rh * (rmse_current_height + s * next_level_rmse);
//    avg_linear_loss = ((avg_linear_loss > 1) ? (2 * log(avg_linear_loss) / log(2.0)) : 0);
    double find_leaf_loss = ((2 + R1) * (_n_stages + 1));
    double loss = find_leaf_loss + avg_linear_loss;

    if (print) {// || (h == 0 && n_nodes == 6060606)) {
        cout << "----------------" << endl
             << "h = " << h << ", n_nodes = " << n_nodes << ", next_n_nodes = " << next_n_nodes
             << ", s = " << s << ", rh = " << rh << ", n_stages = " << _n_stages << endl
             << ", rmse_curr_h = " << rmse_current_height << ", rmse_next = " << next_level_rmse << endl
             << ", loss_curr_h_lin = " << loss_curr_h_lin << ", loss_curr_h_bin = " << loss_curr_h_bin << endl
             << ", loss_next_lin = " << loss_next_lin << ", loss_next_bin = " << loss_nex_bin << endl
             << ", avg_linear_loss = " << avg_linear_loss << ", find_leaf_loss = " << find_leaf_loss << ", total_loss = " << loss << endl
             << "----------------" << endl;
    }

    return loss;
}


double loss_est(int h, long n_nodes, long last_n_nodes, double _rho, double rmse_current_height, bool print=false) {

    double avg_fan = 1.0 * last_n_nodes / n_nodes;
    double _n_stages = log(n_nodes) / log(avg_fan);
//    if (h == 0) {
//        _n_stages = log(n_nodes) / log(avg_fan * 2);
//    }

    if (print) {
        cout << "h = " << h << ", _n_stages = " << _n_stages << endl;
    }
    int n_stages = static_cast<int>(_n_stages);

    double residual = _n_stages - n_stages;
    double next_level_rmse = rmse_current_height;

    double s = 1;
    if (n_stages == 0) {
        for (int i = 0; i <= n_stages; ++i) {
            s *= _rho;
        }
        double last_item = s * residual;
        s = (s - _rho) / (_rho - 1) + last_item;
    } else {
        s = _rho * residual;
    }

    double rh = 1;
    for (int i = 0; i < h; ++i) {
        rh *= _rho;
    }
//    assert(n_stages >= 1);


    assert(s >= 0);
    assert(next_level_rmse >= 0);
    const int n_pairs_each_cache_line = 4;

    double loss_curr_h_lin = rmse_current_height / n_pairs_each_cache_line;
    double loss_next_lin = s * next_level_rmse / n_pairs_each_cache_line;

//    double loss_curr_h_bin = ((rmse_current_height > 2) ? (2 * log(rmse_current_height) / log(2.0)) : 2);
//    double loss_nex_bin = s * ((next_level_rmse > 2) ? (2 * log(next_level_rmse) / log(2.0)) : 2);

    double loss_curr_h_bin = 2 * log(1 + rmse_current_height) / log(2.0);
    double loss_nex_bin = 2 * log(1 + next_level_rmse) / log(2.0);

//    double loss_curr_h = std::min<double>(loss_curr_h_lin, loss_curr_h_bin);
//    double loss_next = std::min<double>(loss_next_lin, loss_nex_bin);

    double loss_curr_h = loss_curr_h_bin;
    double loss_next = loss_nex_bin;

    double avg_linear_loss = (1 + R2) * rh * (loss_curr_h + loss_next);
//    double avg_linear_loss = rh * (rmse_current_height + s * next_level_rmse);
//    avg_linear_loss = ((avg_linear_loss > 1) ? (2 * log(avg_linear_loss) / log(2.0)) : 0);
    double find_leaf_loss = ((2 + R1) * (_n_stages + 1));
    double loss = find_leaf_loss + avg_linear_loss;

    if (print) {// || (h == 0 && n_nodes == 6060606)) {
        cout << "----------------" << endl
             << "h = " << h << ", n_nodes = " << n_nodes << ", last_n_nodes = " << last_n_nodes
             << ", s = " << s << ", rh = " << rh << ", n_stages = " << _n_stages << endl
             << ", rmse_curr_h = " << rmse_current_height << ", rmse_next = " << next_level_rmse << endl
             << ", loss_curr_h_lin = " << loss_curr_h_lin << ", loss_curr_h_bin = " << loss_curr_h_bin << endl
             << ", loss_next_lin = " << loss_next_lin << ", loss_next_bin = " << loss_nex_bin << endl
             << ", avg_linear_loss = " << avg_linear_loss << ", find_leaf_loss = " << find_leaf_loss << ", total_loss = " << loss << endl
             << "----------------" << endl;
    }

    return loss;
}
/*
 *
 */

long estimate_ideal_layout_complex(int h, long N_last, doubleVec &rmses, long N_at_h0, doubleVec &h0_rmses,
                           long min_leaf_fan, long max_leaf_fan, long min_int_fan, long max_int_fan, double _rho, long &best_next_n_nodes, double &best_loss) {

    long start_idx = N_last / max_leaf_fan;
    long end_idx = N_last / min_leaf_fan;
//    assert(end_idx < rmses.size());
    end_idx = std::min<long>(end_idx, rmses.size() - 1);


    long best_n_nodes = start_idx;
    best_loss = 1e50;
    best_next_n_nodes = start_idx / max_int_fan;

    int tmp = 1;
    for (int i = 0; i < h; ++i) {
        tmp *= 10;
    }

    int step = 1000 / tmp;
//    cout << "start_idx = " << start_idx << ", end_idx = " << end_idx << ", step = " << step << endl;


    for (long i = start_idx; i <= end_idx; i += step) { //i: curr_n_nodes
//        double avg_fan = 1.0 * N / i;
////        long leaf_fan = static_cast<long>(avg_fan);
////        long _min_int_fan = std::max<long>(leaf_fan, min_int_fan);
////        long _max_int_fan = std::min<long>(max_int_fan, i / leaf_fan);
//
//        long idx = static_cast<long>(1.0 * N_at_h0 / avg_fan);
//        double rmse = 1e50;
//        if (idx < h0_rmses.size()) {
//            rmse = h0_rmses[idx];
//        }

        long best_next_n_nodes_i = i / max_int_fan;
        double best_loss_i = 1e50;


        long last_next_n_nodes = -1;
        for (long int_fan = max_int_fan; int_fan >= min_int_fan; --int_fan) {

            long next_n_nodes = i / int_fan;
            if (next_n_nodes > 0 && next_n_nodes != last_next_n_nodes) {
//                double loss = loss_est(i, N, h, rmses, next_n_nodes, r, rmse);
                double loss = loss_est_complex(h, i, next_n_nodes, N_at_h0, h0_rmses, _rho, rmses[i]);
                if (loss < best_loss_i) {
                    best_next_n_nodes_i = next_n_nodes;
                    best_loss_i = loss;

//                    if (rmses[i] > 60) {
//                        cout << "+++++++++i = " << i << ", best_loss_i = " << best_loss_i << ", rmse = " << rmses[i] << ", best_next_n_nodes_i = " << best_next_n_nodes_i << endl;
//                    }
                }
            }
            last_next_n_nodes = next_n_nodes;
        }
//        if (i % step == 0) {
//        cout << "**i = " << i << ", loss = " << best_loss_i << ", rmse = " << rmses[i] << ", best_next_n_nodes_i = " << best_next_n_nodes_i << endl;
//        }
        if (best_loss_i < best_loss) {
            best_loss = best_loss_i;
            best_n_nodes = i;
            best_next_n_nodes = best_next_n_nodes_i;
//            if (h == 0) {
            loss_est_complex(h, i, best_next_n_nodes, N_at_h0, h0_rmses, _rho, rmses[i], true);
//            }

//            cout << "i = " << i << ", loss = " << best_loss << ", rmse = " << rmses[i] << ", best_next_n_nodes = " << best_next_n_nodes << endl;
//            cout << "i = " << i << ", loss = " << best_loss_i << ", rmse = " << rmse << ", min_next = " << i / max_int_fan << ", max_next = " << i / min_int_fan << endl;
        }
    }

    cout << "best_n_nodes = " << best_n_nodes << ", best_loss = " << best_loss << endl;

    start_idx = std::max<long>(start_idx, best_n_nodes - step);
    end_idx = std::min<long>(end_idx, best_n_nodes + step);
    for (long i = start_idx; i <= end_idx; ++i) {
//        double avg_fan = 1.0 * N / i;
//
//        long idx = static_cast<long>(1.0 * N_at_h0 / avg_fan);
//        double rmse = 1e50;
//        if (idx < h0_rmses.size()) {
//            rmse = h0_rmses[idx];
//        }

//        long _min_int_fan = std::max<long>(leaf_fan, min_int_fan);
//        long _max_int_fan = std::min<long>(max_int_fan, i / leaf_fan);

        long best_next_n_nodes_i = i / max_int_fan;
        double best_loss_i = 1e50;

        long last_next_n_nodes = -1;
        for (long int_fan = max_int_fan; int_fan >= min_int_fan; --int_fan) {
            long next_n_nodes = i / int_fan;
            if (next_n_nodes > 0 && next_n_nodes != last_next_n_nodes) {
//            double loss = loss_est(i, N, h, rmses, next_n_nodes, r, rmse);
                double loss = loss_est_complex(h, i, next_n_nodes, N_at_h0, h0_rmses, _rho, rmses[i]);
                if (loss < best_loss_i) {
                    best_next_n_nodes_i = next_n_nodes;
                    best_loss_i = loss;
                }
            }
            last_next_n_nodes = next_n_nodes;
        }
//        if (i % 1000 == 0) {
//            cout << "i = " << i << ", loss = " << best_loss_i << ", rmse = " << rmse << ", min_next = " << i / max_int_fan << ", max_next = " << i / min_int_fan << endl;
//        }
        if (best_loss_i < best_loss) {
            best_loss = best_loss_i;
            best_n_nodes = i;
            best_next_n_nodes = best_next_n_nodes_i;
        }
    }



//    for (long i = start_idx; i <= end_idx; ++i) {
//        double min_fan = 1.0 * N / i;
//        long l = i / max_int_fan;
//        long r = i / min_fan;
//        double rmse = rmses[i];
//
//        long best_next_n_nodes_i = l;
//        double best_loss_i = 1e50;
//
//        while (true) {
//            long middle = (l + r) / 2;
//            double middle_loss = loss_est(i, N, rmses, middle, r, rmse);
//
//            long middle_left = middle - 1;
//            double middle_left_loss = loss_est(i, N, rmses, middle_left, r, rmse);
//
//            long middle_right = middle + 1;
//            double middle_right_loss = loss_est(i, N, rmses, middle_right, r, rmse);
//
//            if (middle_loss < middle_left_loss && middle_loss < middle_right_loss) {
//                best_next_n_nodes_i = middle;
//                best_loss_i = middle_loss;
//                break;
//            } else if (middle_loss < middle_right_loss) {
//                l = middle_right;
//
//            } else {
//                r = middle_left;
//            }
//            if (l < r - 2) {
//                best_next_n_nodes_i = middle;
//                best_loss_i = middle_loss;
//                break;
//            }
//        }
//        if (i % 1000 == 0) {
//            cout << "i = " << i << ", loss = " << best_loss_i << endl;
//        }
//        if (best_loss_i < best_loss) {
//            best_loss = best_loss_i;
//            best_n_nodes = i;
//            best_next_n_nodes = best_next_n_nodes_i;
//        }
//    }

    return best_n_nodes;
}


long estimate_ideal_layout(int h, long N_last, doubleVec &rmses, long N_at_h0, doubleVec &h0_rmses,
        long min_leaf_fan, long max_leaf_fan, long min_int_fan, long max_int_fan, double _rho, long &best_next_n_nodes, double &best_loss) {

    long start_idx = N_last / max_leaf_fan;
    long end_idx = N_last / min_leaf_fan;
    end_idx = MIN_LONG(end_idx, rmses.size() - 1);


    long best_n_nodes = start_idx;
    best_loss = 1e50;

    int tmp = 1;
    for (int i = 0; i < h; ++i) {
        tmp *= 10;
    }
    int step = 1000 / tmp;
//    cout << "start_idx = " << start_idx << ", end_idx = " << end_idx << ", step = " << step
//         << ", N_last = " << N_last << ", rmses.size() = " << rmses.size() << ", min_leaf_fan = " << min_leaf_fan << endl;

    for (long i = start_idx; i <= end_idx; i += step) { //i: curr_n_nodes
        bool print = false;
        if (i == 1990000 || i == 10000000) {
            print = true;
        }
        double loss = loss_est(h, i, N_last, _rho, rmses[i], print);
        if (loss < best_loss) {
            best_loss = loss;
            best_n_nodes = i;
//            cout << "best_n_nodes = " << best_n_nodes << ", best_loss = " << best_loss << endl;
        }
        if (i == 1990000 || i == 10000000) {
            cout << "-----i = " << i << ", rmse = " << rmses[i] << ", loss = " << loss << ", h = " << h << endl;
        }
    }

    start_idx = std::max<long>(start_idx, best_n_nodes - step);
    end_idx = std::min<long>(end_idx, best_n_nodes + step);
    for (long i = start_idx; i <= end_idx; ++i) {
        double loss = loss_est(h, i, N_last, _rho, rmses[i]);
        if (loss < best_loss) {
            best_loss = loss;
            best_n_nodes = i;
//            cout << "best_n_nodes = " << best_n_nodes << ", best_loss = " << best_loss << endl;
        }
    }
    return best_n_nodes;
}


void build_ideal_mirror(const keyArray &X, const doubleArray &probs, long N, l_matrix &mirror, const string &mirror_dir, int interval_type) {
    mirror.clear();
    int dir_status = file_utils::path_status(mirror_dir);

    if (dir_status > 1) {
        cout << mirror_dir << " exists and it is not a directory path." << endl;
        return;
    }
    else if (dir_status == 0) {
        file_utils::detect_and_create_dir(mirror_dir);
    } else {
        restore_mirror(mirror_dir, mirror, true);
    }

    if (mirror.size() > 0 && mirror[mirror.size()-1].size() == 1) {
        return;
    }
//    const double r = 0.2;

//    assert(wh == 0.1);
    longVec h0_borders;
    doubleVec h0_rmses;
    if (mirror.size() <= 0) {
        bool restore_status = restore_complete_borders_and_losses(mirror_dir, 0, h0_borders, h0_rmses);
        if (!restore_status) {
//            get_complete_partition_borders(X, probs, N, 0, 32, fanThreashold/2, h0_borders, h0_rmses, interval_type);
            cout << "Building BU-Tree...... This step may take several minutes but will be executed once only." << endl;

            auto start = chrono::system_clock::now();
            get_complete_partition_borders(X, probs, N, 0, buMinFan, fanThreashold/2, h0_borders, h0_rmses, interval_type);
//            get_complete_partition_borders_for_lowest_layer(X, N, 0, buMinFan, fanThreashold/2, h0_borders, h0_rmses, interval_type);
//            get_complete_partition_borders_w_sampling(X, N, 0, buMinFan, fanThreashold/2, h0_borders, h0_rmses, interval_type);
            auto stop = chrono::system_clock::now();
            auto duration = chrono::duration_cast<chrono::seconds>(stop - start);
            double response_time = static_cast<double>(duration.count());
            cout << "get_complete_partition_borders takes " << response_time << " seconds." << endl;
            string h0_borders_path = mirror_dir + "/h0_borders";
            data_utils::save_vec_data(h0_borders_path.c_str(), h0_borders);

            string h0_rmses_path = mirror_dir + "/h0_rmses";
            data_utils::save_vec_data(h0_rmses_path.c_str(), h0_rmses);
        }

        int xcount = 0;
        for (size_t i = 0; i < h0_rmses.size(); ++i) {
            if(h0_rmses[i] > 1e10) {
                ++xcount;
            }
        }
//        cout << "######, xcount = " << xcount << endl;

        long h0_min_size = h0_borders.back();
        h0_borders.pop_back();
//        cout << "h0_borsers.size() = " << h0_borders.size() << ", h0_rmses.size = " << h0_rmses.size() << endl;

        long est_next_n_nodes = 0;
        double h0_loss = 0;
        long est_n_nodes = estimate_ideal_layout(0, N, h0_rmses, N, h0_rmses, buMinFan, 200, 100, 10000, RHO, est_next_n_nodes, h0_loss);

        mirror.push_back(longVec());
        assert(h0_min_size < est_n_nodes);

        cout << "h0, est_n_nodes = " << est_n_nodes << endl;
        mirror[0].insert(mirror[0].end(), h0_borders.begin(), h0_borders.begin() + est_n_nodes);
        sort(mirror[0].begin(), mirror[0].end());
//        string h0_path = mirror_dir + "/h0.dat";
//        data_utils::save_vec_data(h0_path.c_str(), mirror[0]);
    }

//    cout << "BU leaf nodes have been created. #nodes = " << mirror[0].size() << endl;

    longVec root_info({X[0]});
    while (true){
        int curr_h = mirror.size();
        longVec &highest_bounds = mirror[curr_h-1];
        int curr_n_nodes = highest_bounds.size();

        keyArray topX = make_unique<keyType []>(curr_n_nodes+1);
        std::copy(highest_bounds.begin(), highest_bounds.end(), topX.get());
        topX[curr_n_nodes] = X[N-1] + 1;

        if (curr_h > 1 && curr_n_nodes <= 10000) {
            mirror.push_back(root_info);
            break;
        }

        linearRegressor root_lr;
        root_lr.init(topX.get(), 0, 0, curr_n_nodes);
        double root_rmse = root_lr.cal_loss(topX.get(), curr_n_nodes);
        root_rmse = sqrt(root_rmse / curr_n_nodes);
        double root_lin_loss = root_rmse / 4;
        double root_bin_loss = (root_rmse > 2) ? 2 * log(root_rmse) / log(2.0) : 2;
        double root_loss = std::min<double>(root_lin_loss, root_bin_loss);
        root_loss *= (1 + R2) * pow(RHO, curr_h);
//        cout << "curr_h = " << curr_h << ", root_loss = " << root_loss << endl;

        longVec hi_borders;
        doubleVec hi_rmses;
        bool restore_status = restore_complete_borders(mirror_dir, curr_h, hi_borders);
        if (!restore_status) {
            get_complete_partition_borders(topX, nullptr, curr_n_nodes, curr_h, 100, fanThreashold/2, hi_borders, hi_rmses, interval_type);
//            get_complete_partition_borders(topX, NULL, curr_n_nodes, curr_h, buMinFan * 2, fanThreashold/2, hi_borders, hi_rmses, interval_type);
        }

        long hi_min_size = hi_borders.back();
        hi_borders.pop_back();
//        cout << "h" << curr_h << "_borsers.size() = " << hi_borders.size() << ", h" << curr_h << "_rmses.size = " << hi_rmses.size() << endl;

        long est_next_n_nodes = 0;
        double hi_loss = 0;
        long est_n_nodes = estimate_ideal_layout(curr_h, curr_n_nodes, hi_rmses, N, h0_rmses, 100, 10000, 100, 10000, RHO, est_next_n_nodes, hi_loss);
//        long est_n_nodes = estimate_ideal_layout(curr_h, curr_n_nodes, hi_rmses, N, h0_rmses, buMinFan, 10000, 100, 10000, RHO, est_next_n_nodes, hi_loss);


        if (curr_h > 1 && root_loss < hi_loss) {
            mirror.push_back(root_info);
            break;
        } else {
            mirror.push_back(longVec());
//            assert(hi_min_size < est_n_nodes);

//            cout << "h" << curr_h <<", est_n_nodes = " << est_n_nodes << endl;
            mirror[curr_h].insert(mirror[curr_h].end(), hi_borders.begin(), hi_borders.begin() + est_n_nodes);
            sort(mirror[curr_h].begin(), mirror[curr_h].end());
        }
    }

}

void restore_mirror(const string &mirror_dir, l_matrix &mirror, bool ideal) {
    int h = 0;
    while (true) {
        string hi_path = mirror_dir + "/h" + to_string(h) + ".dat";
        if (ideal) {
            hi_path = mirror_dir + "/h" + to_string(h) + "_ideal.dat";
        }


        int hi_status = file_utils::path_status(hi_path);
//        cout << "hi_path = " << hi_path << ", status = " << hi_status << endl;
        if (hi_status > 1) {
            mirror.push_back(longVec());
            data_utils::load_vec_data(hi_path.c_str(), mirror[h]);
            longVec &lv = mirror[h];
//            long *tmp_data = new long[lv.size()];
//            std::copy(lv.begin(), lv.end(), tmp_data);
//            check_order(tmp_data, lv.size());
//            cout << "h = " << h << ", lv.size = " << lv.size() << endl;
            ++h;
        } else {
            break;
        }
    }
}


bool restore_complete_borders(const string &mirror_dir, const int h, longVec &borders) {
    string hi_path = mirror_dir + "/h" + to_string(h) + "_borders";
    int hi_status = file_utils::path_status(hi_path);
    if (hi_status > 1) {
        borders.clear();
        data_utils::load_vec_data(hi_path.c_str(), borders);
        return true;
    }
    return false;
}

bool restore_complete_borders_and_losses(const string &mirror_dir, const int h, longVec &borders, doubleVec &losses) {
    bool status = restore_complete_borders(mirror_dir, h, borders);
    if(!status) {
        return false;
    }
    string hi_path = mirror_dir + "/h" + to_string(h) + "_rmses";
    int hi_status = file_utils::path_status(hi_path);
    if (hi_status > 1) {
        losses.clear();
        data_utils::load_vec_data(hi_path.c_str(), losses);
        return true;
    }
    return false;
}
