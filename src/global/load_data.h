
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <cassert>
#include "global.h"
#include "utils.h"
using namespace std;

#ifndef DILI_LOAD_DATA_H
#define DILI_LOAD_DATA_H


long *allKeys = NULL;
long *allPayloads = NULL;
long *keys_in_diff_dist = NULL;
long *payloads_in_diff_dist = NULL;
long *eraseKeys = NULL;
long *erasePayloads = NULL;
long *random_indices = NULL;
long *first_half_sorted_indices = NULL;
long *first_n_quaters_sorted_indices = NULL;
long *second_half_random_indices = NULL;

long *query_keys = NULL;
long *query_ranges = NULL;
long *payloads_of_ranges = NULL;
long *query_lens = NULL;
long *range_query_payloads = NULL;
long *query_payloads = NULL;
long *auxiliary_keys = NULL;

/*
long *query_range_lens = NULL;
long *query_ranges = NULL;
*/



namespace load_data {

    void allocate() {
//        cout << "calling load_data::allocate()" << endl;
        query_keys = new long[n_query_keys];
        query_payloads = new long[n_query_keys];
        random_indices = new long[n_query_keys];
        first_half_sorted_indices = new long[halfN];
        allKeys = new long[totalDataSize + 2];
        allPayloads = new long[totalDataSize + 2];
        second_half_random_indices = new long[halfN];
    }

    void save_lookups(const char *data_path) {
        FILE *fp = NULL;

        if (NULL == (fp = fopen(data_path, "wb"))) {
            printf("%s cannot be created.\n", data_path);
            exit(1);
        }

        fwrite(&n_query_keys, sizeof(long), 1, fp);

        long *lookups = new long[2 * n_query_keys];
        for (long i = 0; i < n_query_keys; ++i) {
            lookups[2*i] = query_keys[i];
            lookups[2*i+1] = query_payloads[i];
        }

        long rc = fwrite(lookups, sizeof(long), 2 * n_query_keys, fp);
        assert(rc == 2 * n_query_keys);
        fclose(fp);
        delete[] lookups;
    }

    void free_auxiliary_variables() {
        if (random_indices) {
            delete[] random_indices;
            random_indices = NULL;
        }
        if (first_half_sorted_indices) {
            delete[] first_half_sorted_indices;
            first_half_sorted_indices = NULL;
        }
        if (first_n_quaters_sorted_indices) {
            delete[] first_n_quaters_sorted_indices;
            first_n_quaters_sorted_indices = NULL;
        }
        if (second_half_random_indices) {
            delete[] second_half_random_indices;
            second_half_random_indices = NULL;
        }
        if (auxiliary_keys) {
            delete[] auxiliary_keys;
            auxiliary_keys = NULL;
        }
    }

    void free_space_except_query_data() {
        free_auxiliary_variables();
        if (allKeys) {
            delete[] allKeys;
            allKeys = NULL;
        }
        if (allPayloads) {
            delete[] allPayloads;
            allPayloads = NULL;
        }
        if (eraseKeys) {
            delete[] eraseKeys;
            eraseKeys = NULL;
        }
        if (erasePayloads) {
            delete[] erasePayloads;
            erasePayloads = NULL;
        }
        if (query_ranges) {
            delete[] query_ranges;
            query_ranges = NULL;
        }
        if (query_lens) {
            delete[] query_lens;
            query_lens = NULL;
        }
        if (range_query_payloads) {
            delete[] range_query_payloads;
            range_query_payloads = NULL;
        }
        if (payloads_of_ranges) {
            delete[] payloads_of_ranges;
            payloads_of_ranges = NULL;
        }
        if (keys_in_diff_dist) {
            delete[] keys_in_diff_dist;
            keys_in_diff_dist = NULL;
        }
        if (payloads_in_diff_dist) {
            delete[] payloads_in_diff_dist;
            payloads_in_diff_dist = NULL;
        }
        /*
        if (query_range_lens) {
            delete[] query_range_lens;
            query_range_lens = NULL;
        }
        if (query_ranges) {
            delete[] query_ranges;
            query_ranges = NULL;
        }
        */
    }



    void free_space() {
        free_space_except_query_data();
        if (query_keys) {
            delete[] query_keys;
            query_keys = NULL;
        }
        if (query_payloads) {
            delete[] query_payloads;
            query_payloads = NULL;
        }
    }


    void bulk_loading_data_load(const string &data_path, const string &random_indices_path,
                                const string &first_half_sorted_indices_path) {
        load_data_int(random_indices_path.c_str(), random_indices);
        load_data_int(first_half_sorted_indices_path.c_str(), first_half_sorted_indices);
        load_data_int(data_path.c_str(), allKeys);
        long upper_bound = allKeys[totalDataSize - 1] + 1;
        allKeys[totalDataSize] = upper_bound;

        for (long i = 0; i < totalDataSize; ++i) {
            allPayloads[i] = i;
        }
        for (int i = 0; i < n_query_keys; ++i) {
            long j = random_indices[i] * query_step + query_start_idx;
            long idx = first_half_sorted_indices[j];
            query_keys[i] = allKeys[idx];
            query_payloads[i] = allPayloads[idx];
        }

        free_auxiliary_variables();
    }


    void range_query_data_load(const string &data_path, const string &range_query_start_idxes_path, const string &range_query_lens_path) {
        query_ranges = new long[2 * n_query_keys];
//        payloads_of_ranges = new long[2 * n_query_keys];
        query_lens = new long[2 * n_query_keys];

        if (query_keys) {
            delete[] query_keys;
            query_keys = NULL;
        }
        if (query_payloads) {
            delete[] query_payloads;
            query_payloads = NULL;
        }

        load_data_int(range_query_start_idxes_path.c_str(), random_indices);
        load_data_int(range_query_lens_path.c_str(), query_lens);
        load_data_int(data_path.c_str(), allKeys);
        long upper_bound = allKeys[totalDataSize - 1] + 1;
        allKeys[totalDataSize] = upper_bound;

        for (long i = 0; i < totalDataSize; ++i) {
            allPayloads[i] = i;
        }
        for (int i = 0; i < n_query_keys; ++i) {
            long j = random_indices[i];
            long k = j + query_lens[i];
            query_ranges[2 * i] = allKeys[j];
            query_ranges[2 * i + 1] = allKeys[k];

//            payloads_of_ranges[2 * i] = allPayloads[j];
//            payloads_of_ranges[2 * i + 1] = allPayloads[k];
        }

//        free_auxiliary_variables();
        range_query_payloads = new long[totalDataSize];
    }

    void different_card_data_load(const string &data_path, const string &random_indices_path, const string &first_half_sorted_indices_path,
                             const string &first_n_quaters_sorted_indices_path, int n_quarters) {
        assert(n_quarters >= 1 && n_quarters < 4);
        auxiliary_keys = new long[totalDataSize + 2];
        first_n_quaters_sorted_indices = new long[totalDataSize];
        load_data_int(random_indices_path.c_str(), random_indices);
        load_data_int(first_half_sorted_indices_path.c_str(), first_half_sorted_indices);
        long _n = load_data_int(first_n_quaters_sorted_indices_path.c_str(), first_n_quaters_sorted_indices);
        load_data_int(data_path.c_str(), auxiliary_keys);
        long upper_bound = auxiliary_keys[totalDataSize - 1] + 1;
        auxiliary_keys[totalDataSize] = upper_bound;

        long maxIdx = n_quarters * halfN / 2;
        assert(_n == maxIdx);

        for (long i = 0; i < maxIdx; ++i) {
            long idx = first_n_quaters_sorted_indices[i];
            allKeys[i] = auxiliary_keys[idx];
            allPayloads[i] = idx;
        }
        allKeys[maxIdx] = upper_bound;
        allPayloads[maxIdx] = -1;


        if (n_quarters >= 2) {
            for (int i = 0; i < n_query_keys; ++i) {
                long j = random_indices[i] * query_step + query_start_idx;
                long idx = first_half_sorted_indices[j];
                query_keys[i] = auxiliary_keys[idx];
                query_payloads[i] = idx;
            }
        } else {
            for (int i = 0; i < n_query_keys; ++i) {
                long j = random_indices[i] * query_step + query_start_idx;
                long idx = j / 2;
                assert(idx < halfN / 2);
                query_keys[i] = allKeys[idx];
                query_payloads[i] = allPayloads[idx];
            }
        }

        free_auxiliary_variables();
    }

    void insertion_data_load(const string &data_path, const string &random_indices_path,
                             const string &first_half_sorted_indices_path,
                             const string &second_half_random_indices_path) {
        auxiliary_keys = new long[totalDataSize + 2];
        load_data_int(random_indices_path.c_str(), random_indices);
        load_data_int(first_half_sorted_indices_path.c_str(), first_half_sorted_indices);
        load_data_int(second_half_random_indices_path.c_str(), second_half_random_indices);
        load_data_int(data_path.c_str(), auxiliary_keys);
        long upper_bound = auxiliary_keys[totalDataSize - 1] + 1;
        allKeys[totalDataSize] = upper_bound;

        for (long i = 0; i < halfN; ++i) {
            long idx = first_half_sorted_indices[i];
            allKeys[i] = auxiliary_keys[idx];
            allPayloads[i] = idx;
        }
        allKeys[halfN] = upper_bound;
        allPayloads[halfN] = -1;
        for (long i = 0; i < halfN; ++i) {
            long idx = second_half_random_indices[i];
            allKeys[halfN + 1 + i] = auxiliary_keys[idx];
            allPayloads[halfN + 1 + i] = idx;
        }

        allKeys[totalDataSize + 1] = upper_bound;
        allPayloads[totalDataSize + 1] = -1;

        for (int i = 0; i < n_query_keys; ++i) {
            long j = random_indices[i] * query_step + query_start_idx;
            query_keys[i] = allKeys[j];
            query_payloads[i] = allPayloads[j];
        }

        free_auxiliary_variables();
    }

    void deletion_data_load(const string &data_path, const string &random_indices_path,
                             const string &first_half_sorted_indices_path,
                             const string &second_half_random_indices_path) {
        eraseKeys = new long[halfN];
        erasePayloads = new long[halfN];
        load_data_int(random_indices_path.c_str(), random_indices);
        load_data_int(first_half_sorted_indices_path.c_str(), first_half_sorted_indices);
        load_data_int(second_half_random_indices_path.c_str(), second_half_random_indices);
        load_data_int(data_path.c_str(), allKeys);
        long upper_bound = allKeys[totalDataSize - 1] + 1;
        allKeys[totalDataSize] = upper_bound;

        for (long i = 0; i < totalDataSize; ++i) {
            allPayloads[i] = i;
        }
        for (long i = 0; i < halfN; ++i) {
            long idx = second_half_random_indices[i];
            eraseKeys[i] = allKeys[idx];
            erasePayloads[i] = allPayloads[idx];
        }


        for (int i = 0; i < n_query_keys; ++i) {
            long j = random_indices[i] * query_step + query_start_idx;
            long idx = first_half_sorted_indices[j];
            query_keys[i] = allKeys[idx];
            query_payloads[i] = allPayloads[idx];
        }
        free_auxiliary_variables();
    }

    void two_distribution_data_load(const string &data_path, const string &data_2_path, const string &random_indices_path,
                             const string &first_half_sorted_indices_path,
                             const string &second_half_random_indices_path) {
        keys_in_diff_dist = new long[halfN+1];
        payloads_in_diff_dist = new long[halfN+1];
        load_data_int(random_indices_path.c_str(), random_indices);
        load_data_int(first_half_sorted_indices_path.c_str(), first_half_sorted_indices);
        load_data_int(data_path.c_str(), allKeys);
        long upper_bound = allKeys[totalDataSize - 1] + 1;
        allKeys[totalDataSize] = upper_bound;
        allPayloads[totalDataSize] = -1;

        load_data_int(data_2_path.c_str(), keys_in_diff_dist);
        for (long i = 0; i < halfN; ++i) {
            payloads_in_diff_dist[i] = totalDataSize + i;
        }

        for (int i = 0; i < n_query_keys; ++i) {
            long j = random_indices[i] * query_step + query_start_idx;
            long idx = first_half_sorted_indices[j];
            query_keys[i] = allKeys[idx];
            query_payloads[i] = allPayloads[idx];
        }

        free_auxiliary_variables();
    }

//    void two_distribution_data_load(const string &data_path, const string &data_2_path, const string &random_indices_path,
//                                    const string &first_half_sorted_indices_path,
//                                    const string &second_half_random_indices_path) {
//        auxiliary_keys = new long[totalDataSize + 2];
//        load_data_int(random_indices_path.c_str(), random_indices);
//        load_data_int(first_half_sorted_indices_path.c_str(), first_half_sorted_indices);
//        load_data_int(second_half_random_indices_path.c_str(), second_half_random_indices);
//        load_data_int(data_path.c_str(), auxiliary_keys);
//        long upper_bound = auxiliary_keys[totalDataSize - 1] + 1;
//        auxiliary_keys[totalDataSize] = upper_bound;
//
//        for (long i = 0; i < halfN; ++i) {
//            long idx = first_half_sorted_indices[i];
//            allKeys[i] = auxiliary_keys[idx];
//            allPayloads[i] = idx;
//        }
//        allKeys[halfN] = upper_bound;
//        allPayloads[halfN] = -1;
//
//        load_data_int(data_2_path.c_str(), allKeys + halfN + 1);
//        for (long i = 0; i < halfN; ++i) {
//            allPayloads[halfN + 1 + i] = totalDataSize + i;
//        }
//
//        allKeys[totalDataSize + 1] = upper_bound;
//        allPayloads[totalDataSize + 1] = -1;
//
//        for (int i = 0; i < n_query_keys; ++i) {
//            long j = random_indices[i] * query_step + query_start_idx;
//            query_keys[i] = allKeys[j];
//            query_payloads[i] = allPayloads[j];
//        }
//
//        free_auxiliary_variables();
//    }
}
#endif //DILI_LOAD_DATA_H
