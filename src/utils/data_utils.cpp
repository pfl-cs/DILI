#include "data_utils.h"
#include <cstdlib>
#include <cstdio>
#include <cassert>
#include <sys/types.h>
#include <sys/stat.h>
#include <iostream>
#include <string>
#include <fstream>
//#include <filesystem>
using namespace std;

namespace data_utils {
    extern int32Array array_int32 = nullptr;
    extern uint32Array array_uint32 = nullptr;
    extern int64Array array_int64 = nullptr;
    extern uint64Array array_uint64 = nullptr;
    extern floatArray array_float = nullptr;
    extern doubleArray array_double = nullptr;
    extern keyArray array_key = nullptr;
    extern recordPtrArray array_recordPtr = nullptr;

    long load_data_to_uptr(const char *data_path, DATA_FLAG flag) {
        FILE *fp = NULL;

        if (NULL == (fp = fopen(data_path, "rb"))) {
            printf("%s cannot be opened.\n", data_path);
            exit(1);
        }

        long act_data_size = 0;
        long rc = fread(&act_data_size, sizeof(long), 1, fp);
//    cout << "act_data_size = " << act_data_size << endl;


        if (flag == INT32_FLAG) {
            auto _size = sizeof(int);
            int *data = new int[act_data_size + 2];
            rc = fread(data, _size, act_data_size, fp);
            array_int32 = int32Array(data);
        } else if (flag == UINT32_FLAG) {
            auto _size = sizeof(unsigned int);
            unsigned int *data = new unsigned int[act_data_size + 2];
            rc = fread(data, _size, act_data_size, fp);
            array_uint32 = uint32Array(data);
        } else if (flag == INT64_FLAG) {
            auto _size = sizeof(long);
            long *data = new long[act_data_size + 2];
            rc = fread(data, _size, act_data_size, fp);
            array_int64 = int64Array(data);
        } else if (flag == UINT64_FLAG) {
            auto _size = sizeof(unsigned long);
            unsigned long *data = new unsigned long[act_data_size + 2];
            rc = fread(data, _size, act_data_size, fp);
            array_uint64 = uint64Array(data);
        } else if (flag == FLOAT_FLAG) {
            auto _size = sizeof(float);
            float *data = new float[act_data_size + 2];
            rc = fread(data, _size, act_data_size, fp);
            array_float = floatArray(data);
        } else if (flag == DOUBLE_FLAG) {
            auto _size = sizeof(double);
            double *data = new double[act_data_size + 2];
            rc = fread(data, _size, act_data_size, fp);
            array_double = doubleArray(data);
        } else if (flag == KEY_FLAG) {
            auto _size = sizeof(keyType);
            array_key = make_unique<keyType []>(act_data_size + 2);
            rc = fread(array_key.get(), _size, act_data_size, fp);

//            keyType *data = new keyType[act_data_size + 2];
//            rc = fread(data, _size, act_data_size, fp);
//            cout << "-------------" << endl;
//            array_key = keyArray(data);

            for (int i = 0; i < 5; ++i) {
                cout << array_key[i] << " ";
            }
            cout << endl;
            assert(array_key);
            cout << "-------------" << endl;
        } else if (flag == RECORDPTR_FLAG) {
            auto _size = sizeof(recordPtr);
            recordPtr *data = new recordPtr[act_data_size + 2];
            rc = fread(data, _size, act_data_size, fp);
            array_recordPtr = recordPtrArray(data);
        }

        assert(rc == act_data_size);
        fclose(fp);

        return act_data_size;
    }


    long load_data(const char *data_path, void *data, DATA_FLAG flag) {
        FILE *fp = NULL;

        if (NULL == (fp = fopen(data_path, "rb"))) {
            printf("%s cannot be opened.\n", data_path);
            exit(1);
        }

        long act_data_size = 0;
        long rc = fread(&act_data_size, sizeof(long), 1, fp);
//    cout << "act_data_size = " << act_data_size << endl;

        auto _size = sizeof(int);

        switch (flag) {
            case INT32_FLAG:
                _size = sizeof(int);
                break;
            case UINT32_FLAG:
                _size = sizeof(unsigned int);
                break;
            case INT64_FLAG:
                _size = sizeof(long);
                break;
            case UINT64_FLAG:
                _size = sizeof(unsigned long);
                break;
            case FLOAT_FLAG:
                _size = sizeof(float);
                break;
            case DOUBLE_FLAG:
                _size = sizeof(double);
                break;
            default:
                break;
        }

        rc = fread(data, _size, act_data_size, fp);
        assert(rc == act_data_size);
        fclose(fp);

        return act_data_size;
    }

    long load_data_pair(const char *data_path, std::pair<long, long> *_data) {
        FILE *fp = NULL;

        if (NULL == (fp = fopen(data_path, "rb"))) {
            printf("%s cannot be opened.\n", data_path);
            exit(1);
        }

        long act_data_size = 0;
        long rc = fread(&act_data_size,sizeof(long), 1, fp);
//    cout << "act_data_size = " << act_data_size << endl;
        long *data = new long[2 * act_data_size];
        rc = fread(data, sizeof(long), act_data_size * 2, fp);
        assert(rc == act_data_size * 2);
        fclose(fp);
        for (long i = 0; i < act_data_size; ++i) {
            _data[i].first = data[2 * i];
            _data[i].second = data[2 * i + 1];
        }

        delete[] data;

        return act_data_size;
    }


    void save_data(const char *data_path, void *data, long data_size, DATA_FLAG flag) {
        FILE *fp = NULL;

        if (NULL == (fp = fopen(data_path, "wb"))) {
            printf("%s cannot be created.\n", data_path);
            exit(1);
        }
        auto _size = sizeof(int);

        switch (flag) {
            case INT32_FLAG:
                _size = sizeof(int);
                break;
            case UINT32_FLAG:
                _size = sizeof(unsigned int);
                break;
            case INT64_FLAG:
                _size = sizeof(long);
                break;
            case UINT64_FLAG:
                _size = sizeof(unsigned long);
                break;
            case FLOAT_FLAG:
                _size = sizeof(float);
                break;
            case DOUBLE_FLAG:
                _size = sizeof(double);
                break;
            default:
                break;
        }


        fwrite(&data_size, sizeof(long), 1, fp);

        long rc = fwrite(data, _size, data_size, fp);
        assert(rc == data_size);
        fclose(fp);
    }


    void save_data_pair(const char *data_path, std::pair<long, long> *_data, long data_size) {
        FILE *fp = NULL;

        if (NULL == (fp = fopen(data_path, "wb"))) {
            printf("%s cannot be created.\n", data_path);
            exit(1);
        }

        fwrite(&data_size, sizeof(long), 1, fp);
        long *data = new long[data_size * 2];
        for (long i = 0; i < data_size; ++i) {
            data[2 * i] = _data[i].first;
            data[2 * i + 1] = _data[i].second;
        }

        long rc = fwrite(data, sizeof(long), data_size * 2, fp);
        assert(rc == data_size * 2);
        fclose(fp);
        delete[] data;
    }


    void save_mat_data_double(const char *data_path, const d_matrix &d_mat) {
        long n_rows = d_mat.size();
        if (n_rows > 0) {
            const doubleVec &first_row = d_mat[0];
            long n_cols = first_row.size();
            long data_size = n_rows * n_cols;
            double *d_data = new double[data_size];

            long last_idx = 0;
            for (long i = 0; i < n_rows; ++i) {
                const doubleVec &d_vec = d_mat[i];
                std::copy(d_vec.begin(), d_vec.end(), d_data + last_idx);
                last_idx += d_vec.size();
            }

            FILE *fp = NULL;
            if (NULL == (fp = fopen(data_path, "wb"))) {
                printf("%s cannot be created.\n", data_path);
                exit(1);
            }

            fwrite(&n_rows, sizeof(long), 1, fp);
            fwrite(&n_cols, sizeof(long), 1, fp);

            long rc = fwrite(d_data, sizeof(double), data_size, fp);
            assert(rc == data_size);
            fclose(fp);

            delete[] d_data;
        }

    }

    void scale_sorted_data(long *sorted_data, long data_size, long min_val, long max_val, bool allow_repeat) {
        if (!allow_repeat) {
            long offset = 0;
            for (long i = 1; i < data_size; ++i) {
                long x = sorted_data[i] + offset;
                long diff = x - sorted_data[i - 1];
                if (diff <= 0) {
                    assert(diff == 0);
                    offset += 1;
                }
                sorted_data[i] += offset;
            }
        }
        long data_min = sorted_data[0];
        long data_max = sorted_data[data_size - 1];
        double ratio = 1.0 * (max_val - min_val) / (data_max - data_min);

        if (ratio > 1) {
            for (long i = 1; i < data_size; ++i) {
                sorted_data[i] = static_cast<long>(ratio * (sorted_data[i] - data_min));
            }
            sorted_data[0] = min_val;
            sorted_data[data_size - 1] = max_val;
        }
        check(sorted_data, data_size);
    }


    void check(long *keys, long len) {
        for (int i = 0; i < len - 1; ++i) {
            if (keys[i] >= keys[i+1]) {
                cout << "i = " << i << ", keys[i] = " << keys[i] << ", keys[i+1] = " << keys[i+1] << endl;
            }
            assert(keys[i] < keys[i+1]);
        }
    }


    void generate_data_descs(long *data, long n, const std::string &desc_path, const std::string &meta_info) {
        long min_diff = data[1] - data[0];
        long max_diff = min_diff;
        if (min_diff <= 0) {
            cout << "i = 1, min_diff = " << min_diff << ", data[i - 1] = " << data[0] << ", data[i] = " << data[1]
                 << endl;
        }
        assert(min_diff > 0);
        for (long i = 2; i < n; ++i) {
            long diff = data[i] - data[i - 1];
            if (min_diff <= 0) {
                cout << "i = " << i << ", diff = " << min_diff << ", data[i - 1] = " << data[i - 1] << ", data[i] = "
                     << data[i] << endl;
            }
            assert(diff > 0);
            if (diff > max_diff) {
                max_diff = diff;
            }
            if (diff < min_diff) {
                min_diff = diff;
            }
        }
        ofstream fout(desc_path);
        if (meta_info.size() > 1) {
            fout << meta_info << endl;
        }
        fout << "data_size = " << n << endl
             << "min_key = " << data[0] << endl
             << "max_key = " << data[n - 1] << endl
             << "min_diff = " << min_diff << endl
             << "max_diff = " << max_diff << endl;

        fout.close();
    }
}
