#include "../global/global_typedef.h"
#include <string>
#include <algorithm>
#include <iostream>
#include <vector>
#include <cassert>

#ifndef UTILS_DATA_UTILS_H
#define UTILS_DATA_UTILS_H

namespace data_utils {
    enum DATA_FLAG {
        INT32_FLAG = 0,
        UINT32_FLAG = 1,
        INT64_FLAG = 2,
        UINT64_FLAG = 3,
        FLOAT_FLAG = 4,
        DOUBLE_FLAG = 5,
        KEY_FLAG = 6,
        RECORDPTR_FLAG = 7,
        PAIR_FLAG = 8
    };

    extern int32Array array_int32;
    extern uint32Array array_uint32;
    extern int64Array array_int64;
    extern uint64Array array_uint64;
    extern floatArray array_float;
    extern doubleArray array_double;
    extern keyArray array_key;
    extern recordPtrArray array_recordPtr;

    long load_data_to_uptr(const char *path, DATA_FLAG flag=DATA_FLAG::INT64_FLAG);
    long load_data(const char *path, void *data, DATA_FLAG flag=DATA_FLAG::INT64_FLAG);
    long load_data_pair(const char *data_path, std::pair<long, long> *_data);

    template<typename T>
    long load_vec_data(const char *data_path, std::vector<T> &vec_data) {


        FILE *fp = NULL;

        if (NULL == (fp = fopen(data_path, "rb"))) {
            printf("%s cannot be opened.\n", data_path);
            exit(1);
        }

        long act_data_size = 0;
        long rc = fread(&act_data_size,sizeof(long), 1, fp);
        T *data = new T[act_data_size];

        rc = fread(data, sizeof(T), act_data_size, fp);
        assert(rc == act_data_size);
        fclose(fp);

        vec_data.clear();
        std::copy(data, data + act_data_size, std::back_inserter(vec_data));

        delete[] data;

        return act_data_size;
    }

    void save_data_pair(const char *data_path, std::pair<long, long> *_data, long data_size);

    void save_data(const char *data_path, void *data, long data_size, DATA_FLAG flag=DATA_FLAG::INT64_FLAG);


    template<typename T>
    void save_vec_data(const char *data_path, const std::vector<T> &vec_data) {
        long data_size = vec_data.size();
        T *data = new T[data_size];
        std::copy(vec_data.begin(), vec_data.end(), data);
//        save_data(data_path, data, data_size)

        FILE *fp = NULL;

        if (NULL == (fp = fopen(data_path, "wb"))) {
            printf("%s cannot be created.\n", data_path);
            exit(1);
        }

        fwrite(&data_size, sizeof(long), 1, fp);

        long rc = fwrite(data, sizeof(T), data_size, fp);
        assert(rc == data_size);
        fclose(fp);

        delete[] data;
    }

    template<typename T1, typename T2>
    void save_pair_vec_data(const char *data_path, const std::vector< std::pair<T1, T2> > &vec_data) {
        long data_size = vec_data.size();
        std::unique_ptr<T1 []> data_1 = std::make_unique<T1 []>(data_size);
        std::unique_ptr<T2 []> data_2 = std::make_unique<T2 []>(data_size);
        for (size_t i = 0; i < vec_data.size(); ++i) {
            data_1[i] = vec_data[i].first;
            data_2[i] = vec_data[i].second;
        }

        FILE *fp = NULL;

        if (NULL == (fp = fopen(data_path, "wb"))) {
            printf("%s cannot be created.\n", data_path);
            exit(1);
        }

        fwrite(&data_size, sizeof(long), 1, fp);

        long rc1 = fwrite(data_1.get(), sizeof(T1), data_size, fp);
        long rc2 = fwrite(data_2.get(), sizeof(T2), data_size, fp);
        assert(rc1 == data_size);
        assert(rc2 == data_size);
        fclose(fp);
    }

    template<typename T1, typename T2>
    void load_pair_vec_data(const char *data_path, std::vector< std::pair<T1, T2> > &vec_data) {
        FILE *fp = NULL;

        if (NULL == (fp = fopen(data_path, "rb"))) {
            printf("%s cannot be opened.\n", data_path);
            exit(1);
        }

        long act_data_size = 0;
        long rc = fread(&act_data_size,sizeof(long), 1, fp);
        std::unique_ptr<T1 []> data_1 = std::make_unique<T1 []>(act_data_size);
        std::unique_ptr<T2 []> data_2 = std::make_unique<T2 []>(act_data_size);

        rc = fread(data_1.get(), sizeof(T1), act_data_size, fp);
        assert(rc == act_data_size);
        rc = fread(data_2.get(), sizeof(T2), act_data_size, fp);
        assert(rc == act_data_size);
        fclose(fp);

        vec_data.clear();
        for (long i = 0; i < act_data_size; ++i) {
            vec_data.emplace_back(std::make_pair(data_1[i], data_2[i]));
        }
    }


    void save_mat_data_double(const char *data_path, const d_matrix &d_mat);


    void scale_sorted_data(long *sorted_data, long data_size, long min_val, long max_val, bool allow_repeat=false);

    template<typename T>
    T cal_min_diff(T *data, long n) {
        T min_diff = data[1] - data[0];
        for (size_t i = 2; i < n; ++i) {
            T diff = data[i] - data[i - 1];
            if (diff < min_diff) {
                min_diff = diff;
            }
            if (min_diff <= 0) {
                std::cout << "i = " << i << ", data[i] = " << data[i] << std::endl;
                return min_diff;
            }
        }
        return min_diff;
    }

    template<typename T>
    T max_value(const std::vector<T> &data) {
        T maxv = data[0];
        for (const auto &v : data) {
            if (v > maxv) {
                maxv = v;
            }
        }
        return maxv;
    }
    template<typename T>
    long max_idx(const std::vector<T> &data) {
        T maxv = data[0];
        long idx = 0;
        for (long i = 0; i < data.size(); ++i) {
            const T &v = data[i];
            if (v > maxv) {
                maxv = v;
                idx = i;
            }
        }
        return idx;
    }

    template<typename T>
    T min_value(const std::vector<T> &data) {
        T minv = data[0];
        for (const auto &v : data) {
            if (v > minv) {
                minv = v;
            }
        }
        return minv;
    }

    template<typename T>
    long min_idx(const std::vector<T> &data) {
        T minv = data[0];
        long idx = 0;
        for (long i = 0; i < data.size(); ++i) {
            const T &v = data[i];
            if (v < minv) {
                minv = v;
                idx = i;
            }
        }
        return idx;
    }

    template<typename T>
    void swap(std::vector<T> &data, long idx1, long idx2) {
        T v = data[idx1];
        data[idx1] = data[idx2];
        data[idx2] = v;
    }


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

    void check(long *keys, long len);

    void generate_data_descs(long *data, long n, const std::string &desc_path, const std::string &meta_info);
}


#endif // UTILS_DATA_UTILS_H
