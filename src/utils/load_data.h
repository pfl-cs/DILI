
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <cassert>
#include "../global/global.h"
#include "data_utils.h"
using namespace std;

#ifndef DILI_LOAD_DATA_H
#define DILI_LOAD_DATA_H

namespace dataset {
    keyArray allKeys = nullptr;
    recordPtrArray allPtrs = nullptr;
    keyArray keys_in_diff_dist = nullptr;
    recordPtrArray ptrs_in_diff_dist = nullptr;
    keyArray eraseKeys = nullptr;
    recordPtrArray erasePtrs = nullptr;
    int64Array random_indices = nullptr;
    int64Array first_half_sorted_indices = nullptr;
    int64Array first_n_quaters_sorted_indices = nullptr;
    int64Array second_half_random_indices = nullptr;

    keyArray query_keys = nullptr;
    keyArray query_ranges = nullptr;
    recordPtrArray ptrs_of_ranges = nullptr;
    int64Array query_lens = nullptr;
    recordPtrArray range_query_ptrs = nullptr;
    recordPtrArray query_ptrs = nullptr;
    keyArray auxiliary_keys = nullptr;


//    keyType *allKeys = NULL;
//    recordPtr *allPtrs = NULL;
//    keyType *keys_in_diff_dist = NULL;
//    recordPtr *ptrs_in_diff_dist = NULL;
//    keyType *eraseKeys = NULL;
//    recordPtr *erasePtrs = NULL;
//    long *random_indices = NULL;
//    long *first_half_sorted_indices = NULL;
//    long *first_n_quaters_sorted_indices = NULL;
//    long *second_half_random_indices = NULL;
//
//    keyType *query_keys = NULL;
//    keyType *query_ranges = NULL;
//    recordPtr *ptrs_of_ranges = NULL;
//    long *query_lens = NULL;
//    recordPtr *range_query_ptrs = NULL;
//    recordPtr *query_ptrs = NULL;
//    keyType *auxiliary_keys = NULL;
}

/*
long *query_range_lens = NULL;
long *query_ranges = NULL;
*/

using namespace dataset;


namespace load_data {

    void _params_setting() {
        halfN = totalDataSize / 2;
        query_step = halfN / n_query_keys;
        query_start_idx = query_step / 2;
    }

    void allocate() {
//        cout << "calling load_data::allocate()" << endl;

        query_keys = make_unique<keyType []>(n_query_keys);
        query_ptrs = make_unique<recordPtr []>(n_query_keys);
//        random_indices = new long[n_query_keys];
//        first_half_sorted_indices = new long[halfN];
//        allKeys = new keyType [totalDataSize + 2];
//        second_half_random_indices = new long[halfN];
    }

    void free_auxiliary_variables() {
        random_indices = nullptr;
        first_half_sorted_indices = nullptr;
        first_n_quaters_sorted_indices = nullptr;
        second_half_random_indices = nullptr;
        auxiliary_keys = nullptr;
    }

    void free_space_except_query_data() {
        free_auxiliary_variables();
        allKeys = nullptr;
        allPtrs = nullptr;
        eraseKeys = nullptr;
        erasePtrs = nullptr;
        query_ranges = nullptr;
        query_lens = nullptr;
        range_query_ptrs = nullptr;
        ptrs_of_ranges = nullptr;
        keys_in_diff_dist = nullptr;
        ptrs_in_diff_dist = nullptr;
    }

    void free_space() {
        free_space_except_query_data();
        query_keys = nullptr;
        query_ptrs = nullptr;
    }


    void bulk_loading_data_load(const string &data_path, const string &random_indices_path,
                                const string &first_half_sorted_indices_path) {
        data_utils::load_data_to_uptr(random_indices_path.c_str(), data_utils::INT64_FLAG);
        random_indices = std::move(data_utils::array_int64);

        data_utils::load_data_to_uptr(first_half_sorted_indices_path.c_str(), data_utils::INT64_FLAG);
        first_half_sorted_indices = std::move(data_utils::array_int64);

        totalDataSize = data_utils::load_data_to_uptr(data_path.c_str(), data_utils::KEY_FLAG);
//        std::cout << "totalDataSize = " << totalDataSize << std::endl;
        assert(data_utils::array_key);
        allKeys = std::move(data_utils::array_key);
        halfN = totalDataSize / 2;

        keyType upper_bound = allKeys[totalDataSize - 1] + 1;
        allKeys[totalDataSize] = upper_bound;

        allPtrs = make_unique<recordPtr []>(totalDataSize + 2);
        for (long i = 0; i < totalDataSize; ++i) {
            allPtrs[i] = i;
        }

        for (int i = 0; i < n_query_keys; ++i) {
            long j = random_indices[i] * query_step + query_start_idx;
            long idx = first_half_sorted_indices[j];
            query_keys[i] = allKeys[idx];
            query_ptrs[i] = allPtrs[idx];
        }

        free_auxiliary_variables();

    }


    void range_query_data_load(const string &data_path, const string &range_query_start_idxes_path, const string &range_query_lens_path) {
        query_ranges = make_unique<keyType []>(2 * n_query_keys);
        query_keys = nullptr;
        query_ptrs = nullptr;

        data_utils::load_data_to_uptr(range_query_start_idxes_path.c_str(), data_utils::INT64_FLAG);
        random_indices = std::move(data_utils::array_int64);

        data_utils::load_data_to_uptr(range_query_lens_path.c_str(), data_utils::INT64_FLAG);
        query_lens = std::move(data_utils::array_int64);

        totalDataSize = data_utils::load_data_to_uptr(data_path.c_str(), data_utils::KEY_FLAG);
        allKeys = std::move(data_utils::array_key);

        allPtrs = make_unique<recordPtr []>(totalDataSize + 2);

        keyType upper_bound = allKeys[totalDataSize - 1] + 1;
        allKeys[totalDataSize] = upper_bound;

        for (long i = 0; i < totalDataSize; ++i) {
            allPtrs[i] = i;
        }
        for (int i = 0; i < n_query_keys; ++i) {
            long j = random_indices[i];
            long k = j + query_lens[i];
            query_ranges[2 * i] = allKeys[j];
            query_ranges[2 * i + 1] = allKeys[k];
        }

        range_query_ptrs = make_unique<recordPtr []>(totalDataSize);
    }

    void different_card_data_load(const string &data_path, const string &random_indices_path, const string &first_half_sorted_indices_path,
                             const string &first_n_quaters_sorted_indices_path, int n_quarters) {
        assert(n_quarters >= 1 && n_quarters < 4);
        data_utils::load_data_to_uptr(random_indices_path.c_str(), data_utils::INT64_FLAG);
        random_indices = std::move(data_utils::array_int64);

        long _n = data_utils::load_data_to_uptr(first_half_sorted_indices_path.c_str(), data_utils::INT64_FLAG);
//        long quater_n = _n / n_quarters;
//        long halfN = quater_n * 2;
        first_half_sorted_indices = std::move(data_utils::array_int64);

        data_utils::load_data_to_uptr(first_n_quaters_sorted_indices_path.c_str(), data_utils::INT64_FLAG);
        first_n_quaters_sorted_indices = std::move(data_utils::array_int64);

        totalDataSize = data_utils::load_data_to_uptr(data_path.c_str(), data_utils::KEY_FLAG);
        long halfN = totalDataSize / 2;
        auxiliary_keys = std::move(data_utils::array_key);

        keyType upper_bound = auxiliary_keys[totalDataSize - 1] + 1;
        auxiliary_keys[totalDataSize] = upper_bound;

        long maxIdx = n_quarters * halfN / 2;
//        if (_n != maxIdx) {
//            cout << "_n = " << _n << ", maxIdx = " << maxIdx << endl;
//        }
//        assert(_n == maxIdx);

        allKeys = make_unique<keyType []>(totalDataSize + 2);
        allPtrs = make_unique<recordPtr []>(totalDataSize + 2);
        for (long i = 0; i < maxIdx; ++i) {
            long idx = first_n_quaters_sorted_indices[i];
            allKeys[i] = auxiliary_keys[idx];
            allPtrs[i] = idx;
        }
        allKeys[maxIdx] = upper_bound;
        allPtrs[maxIdx] = -1;


        if (n_quarters >= 2) {
            for (int i = 0; i < n_query_keys; ++i) {
                long j = random_indices[i] * query_step + query_start_idx;
                long idx = first_half_sorted_indices[j];
                query_keys[i] = auxiliary_keys[idx];
                query_ptrs[i] = idx;
            }
        } else {
            for (int i = 0; i < n_query_keys; ++i) {
                long j = random_indices[i] * query_step + query_start_idx;
                long idx = j / 2;
                assert(idx < halfN / 2);
                query_keys[i] = allKeys[idx];
                query_ptrs[i] = allPtrs[idx];
            }
        }

        free_auxiliary_variables();
    }

    void insertion_data_load(const string &data_path, const string &random_indices_path,
                             const string &first_half_sorted_indices_path,
                             const string &second_half_random_indices_path) {
        allKeys = make_unique<keyType []>(totalDataSize + 2);
        allPtrs = make_unique<recordPtr []>(totalDataSize + 2);

        data_utils::load_data_to_uptr(random_indices_path.c_str(), data_utils::INT64_FLAG);
        random_indices = std::move(data_utils::array_int64);

        data_utils::load_data_to_uptr(first_half_sorted_indices_path.c_str(), data_utils::INT64_FLAG);
        first_half_sorted_indices = std::move(data_utils::array_int64);

        data_utils::load_data_to_uptr(second_half_random_indices_path.c_str(), data_utils::INT64_FLAG);
        second_half_random_indices = std::move(data_utils::array_int64);

        totalDataSize = data_utils::load_data_to_uptr(data_path.c_str(), data_utils::KEY_FLAG);
        auxiliary_keys = std::move(data_utils::array_key);
        halfN = totalDataSize / 2;

        keyType upper_bound = auxiliary_keys[totalDataSize - 1] + 1;
        allKeys[totalDataSize] = upper_bound;

        for (long i = 0; i < halfN; ++i) {
            long idx = first_half_sorted_indices[i];
            allKeys[i] = auxiliary_keys[idx];
            allPtrs[i] = idx;
        }
        allKeys[halfN] = upper_bound;
        allPtrs[halfN] = -1;
        for (long i = 0; i < halfN; ++i) {
            long idx = second_half_random_indices[i];
            allKeys[halfN + 1 + i] = auxiliary_keys[idx];
            allPtrs[halfN + 1 + i] = idx;
        }

        allKeys[totalDataSize + 1] = upper_bound;
        allPtrs[totalDataSize + 1] = -1;

        for (int i = 0; i < n_query_keys; ++i) {
            long j = random_indices[i] * query_step + query_start_idx;
            query_keys[i] = allKeys[j];
            query_ptrs[i] = allPtrs[j];
        }

        free_auxiliary_variables();
    }

    void deletion_data_load(const string &data_path, const string &random_indices_path,
                             const string &first_half_sorted_indices_path,
                             const string &second_half_random_indices_path) {


        data_utils::load_data_to_uptr(random_indices_path.c_str(), data_utils::INT64_FLAG);
        random_indices = std::move(data_utils::array_int64);

        data_utils::load_data_to_uptr(first_half_sorted_indices_path.c_str(), data_utils::INT64_FLAG);
        first_half_sorted_indices = std::move(data_utils::array_int64);

        data_utils::load_data_to_uptr(second_half_random_indices_path.c_str(), data_utils::INT64_FLAG);
        second_half_random_indices = std::move(data_utils::array_int64);

        totalDataSize = data_utils::load_data_to_uptr(data_path.c_str(), data_utils::KEY_FLAG);
        allKeys = std::move(data_utils::array_key);
        halfN = totalDataSize / 2;

        keyType upper_bound = allKeys[totalDataSize - 1] + 1;
        allKeys[totalDataSize] = upper_bound;

        allPtrs = make_unique<recordPtr []>(totalDataSize + 2);


        for (long i = 0; i < totalDataSize; ++i) {
            allPtrs[i] = i;
        }

        eraseKeys = make_unique<keyType []>(halfN);
        erasePtrs = make_unique<recordPtr []>(halfN);
        for (long i = 0; i < halfN; ++i) {
            long idx = second_half_random_indices[i];
            eraseKeys[i] = allKeys[idx];
            erasePtrs[i] = allPtrs[idx];
        }


        for (int i = 0; i < n_query_keys; ++i) {
            long j = random_indices[i] * query_step + query_start_idx;
            long idx = first_half_sorted_indices[j];
            query_keys[i] = allKeys[idx];
            query_ptrs[i] = allPtrs[idx];
        }
        free_auxiliary_variables();
    }

    void two_distribution_data_load(const string &data_path, const string &data_2_path, const string &random_indices_path,
                             const string &first_half_sorted_indices_path) {
//        keys_in_diff_dist = new keyType [halfN+1];


        data_utils::load_data_to_uptr(random_indices_path.c_str(), data_utils::INT64_FLAG);
        random_indices = std::move(data_utils::array_int64);
        data_utils::load_data_to_uptr(first_half_sorted_indices_path.c_str(), data_utils::INT64_FLAG);
        first_half_sorted_indices = std::move(data_utils::array_int64);

        totalDataSize = data_utils::load_data_to_uptr(data_path.c_str(), data_utils::KEY_FLAG);
        allKeys = std::move(data_utils::array_key);
        halfN = totalDataSize / 2;

        keyType upper_bound = allKeys[totalDataSize - 1] + 1;
        allKeys[totalDataSize] = upper_bound;

        allPtrs = make_unique<recordPtr []>(totalDataSize + 2);
        allPtrs[totalDataSize] = -1;

        data_utils::load_data_to_uptr(data_2_path.c_str(), data_utils::KEY_FLAG);
        keys_in_diff_dist = std::move(data_utils::array_key);

        ptrs_in_diff_dist = make_unique<recordPtr []>(halfN + 1);
        for (long i = 0; i < halfN; ++i) {
            ptrs_in_diff_dist[i] = totalDataSize + i;
        }

        for (int i = 0; i < n_query_keys; ++i) {
            long j = random_indices[i] * query_step + query_start_idx;
            long idx = first_half_sorted_indices[j];
            query_keys[i] = allKeys[idx];
            query_ptrs[i] = allPtrs[idx];
        }

        free_auxiliary_variables();
    }
}
#endif //DILI_LOAD_DATA_H
