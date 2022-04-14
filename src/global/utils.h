#include "global.h"
#include <string>
#include <algorithm>
#include <iostream>
#include <unistd.h>
#include <sys/ioctl.h>
#include <linux/perf_event.h>
#include <asm/unistd.h>

#ifndef DILI_UTILS_H
#define DILI_UTILS_H

long binarySearch(long *data, long l, long r, unsigned long x);


bool file_exists(const char *path);

float vec_vec_mul(float v0[], float v1[], float bias, int n);
float vec_vec_mul_and_relu(float v0[], float v1[], float bias, int n);
//double vec_vec_mul(d_vector v0, d_vector v1, double bias, int n);
//double vec_vec_mul_and_relu(d_vector v0, d_vector v1, double bias, int n);
long load_data_int(const char *path, long *int_data);
long load_data_uint(const char *data_path, unsigned long *uint_data);
long load_data_pair(const char *data_path, std::pair<long, long> *_data);
void load_data_int_vec_from_txt(const std::string &path, longVec &int_data_vec);
long load_data_int_vec(const char *path, longVec &int_data_vec);

long load_data_float(const char *path, float *float_data);

long load_data_double(const char *path, double *data);
long load_data_double_vec(const char *path, doubleVec &d_vec);

void save_data_int(const char *data_path, long *int_data, long data_size);
void save_data_pair(const char *data_path, std::pair<long, long> *_data, long data_size);
void save_data_int_vec(const char *data_path, const longVec &int_data_vec);

void save_data_double(const char *data_path, double *data, long data_size);
void save_data_double_vec(const char *data_path, const doubleVec &d_vec);

void save_data_double_mat(const char *data_path, const d_matrix &d_mat);

int path_status(const std::string &path_name);
void create_dir(const std::string &dir);

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

void check(long *keys, int len);


// trim from start (in place)
static inline void ltrim(std::string &s) {
    s.erase(s.begin(), std::find_if(s.begin(), s.end(), [](unsigned char ch) {
        return !std::isspace(ch);
    }));
}

// trim from end (in place)
static inline void rtrim(std::string &s) {
    s.erase(std::find_if(s.rbegin(), s.rend(), [](unsigned char ch) {
        return !std::isspace(ch);
    }).base(), s.end());
}

// trim from both ends (in place)
static inline void trim(std::string &s) {
    ltrim(s);
    rtrim(s);
}

// trim from start (copying)
static inline std::string ltrim_copy(std::string s) {
    ltrim(s);
    return s;
}

// trim from end (copying)
static inline std::string rtrim_copy(std::string s) {
    rtrim(s);
    return s;
}

// trim from both ends (copying)
static inline std::string trim_copy(std::string s) {
    trim(s);
    return s;
}

long perf_event_open(struct perf_event_attr *hw_event, pid_t pid,
                     int cpu, int group_fd, unsigned long flags);

void generate_data_descs(long *data, long n, const std::string &desc_path, const std::string &meta_info);

void read_lines(const std::string &path, std::vector<std::string> &lines);
#endif // DILI_UTILS_H
