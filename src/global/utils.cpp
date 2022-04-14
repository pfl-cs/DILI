#include "utils.h"
#include <cstdlib>
#include <cstdio>
#include <cassert>
#include <sys/types.h>
#include <sys/stat.h>
#include "global.h"
#include <iostream>
#include <string>
#include <fstream>
//#include <filesystem>
using namespace std;

bool file_exists(const char *path) {
    FILE *fp = NULL;

    if (NULL == (fp = fopen(path, "rb"))) {
        return false;
    }

    fclose(fp);
    return true;
}

float vec_vec_mul(float v0[], float v1[], float bias, int n) {
    int i;
    float sum = bias;
    for (i = 0; i < n; ++i) {
        sum += v0[i] * v1[i];
    }
    return sum;
}

float vec_vec_mul_and_relu(float v0[], float v1[], float bias, int n) {
    int i;
    float sum = bias;
    for (i = 0; i < n; ++i) {
        sum += v0[i] * v1[i];
    }
    return sum > 0 ? sum : 0;
}


long load_data_int(const char *data_path, long *int_data) {
    FILE *fp = NULL;

    if (NULL == (fp = fopen(data_path, "rb"))) {
        printf("%s cannot be opened.\n", data_path);
        exit(1);
    }

    long act_data_size = 0;
    long rc = fread(&act_data_size,sizeof(long), 1, fp);
//    cout << "act_data_size = " << act_data_size << endl;

    rc = fread(int_data, sizeof(long), act_data_size, fp);
    assert(rc == act_data_size);
    fclose(fp);

    return act_data_size;
}

long load_data_uint(const char *data_path, unsigned long *uint_data) {
    FILE *fp = NULL;

    if (NULL == (fp = fopen(data_path, "rb"))) {
        printf("%s cannot be opened.\n", data_path);
        exit(1);
    }

    long act_data_size = 0;
    long rc = fread(&act_data_size,sizeof(long), 1, fp);
//    cout << "act_data_size = " << act_data_size << endl;

    rc = fread(uint_data, sizeof(unsigned long), act_data_size, fp);
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

long load_data_int_vec(const char *data_path, longVec &int_data_vec) {
    FILE *fp = NULL;

    if (NULL == (fp = fopen(data_path, "rb"))) {
        printf("%s cannot be opened.\n", data_path);
        exit(1);
    }

    long act_data_size = 0;
    long rc = fread(&act_data_size,sizeof(long), 1, fp);
    long *int_data = new long[act_data_size];

    rc = fread(int_data, sizeof(long), act_data_size, fp);
    assert(rc == act_data_size);
    fclose(fp);

    int_data_vec.clear();
    std::copy(int_data, int_data + act_data_size, std::back_inserter(int_data_vec));

    delete[] int_data;

    return act_data_size;
}

void load_data_int_vec_from_txt(const string &path, longVec &int_data_vec) {
    ifstream fin(path, ios::in);
    if (!fin.is_open()) {
        cout << "Cannot open " << path << "." << endl;
    } else {
        while (!fin.eof()) {
            long n;
            fin >> n;
            int_data_vec.push_back(n);
        }
    }
    fin.close();
}



long load_data_float(const char *data_path, float *float_data) {
    FILE *fp = NULL;

    if (NULL == (fp = fopen(data_path, "rb"))) {
        printf("%s cannot be opened.\n", data_path);
        exit(1);
    }

    long act_data_size = 0;
    long rc = fread(&act_data_size,sizeof(long), 1, fp);

    rc = fread(float_data,sizeof(float), act_data_size, fp);
    assert(rc == act_data_size);
    fclose(fp);

    return act_data_size;
}

long load_data_double(const char *data_path, double *data) {
    FILE *fp = NULL;

    if (NULL == (fp = fopen(data_path, "rb"))) {
        printf("%s cannot be opened.\n", data_path);
        exit(1);
    }

    long act_data_size = 0;
    long rc = fread(&act_data_size, sizeof(long), 1, fp);
//    cout << "act_data_size = " << act_data_size << endl;

    rc = fread(data, sizeof(double), act_data_size, fp);
    assert(rc == act_data_size);
    fclose(fp);

    return act_data_size;
}


long load_data_double_vec(const char *data_path, doubleVec &d_vec) {
    FILE *fp = NULL;

    if (NULL == (fp = fopen(data_path, "rb"))) {
        printf("%s cannot be opened.\n", data_path);
        exit(1);
    }

    long act_data_size = 0;
    long rc = fread(&act_data_size,sizeof(long), 1, fp);
    double *d_data = new double[act_data_size];

    rc = fread(d_data, sizeof(double), act_data_size, fp);
    assert(rc == act_data_size);
    fclose(fp);

    d_vec.clear();
    std::copy(d_data, d_data + act_data_size, std::back_inserter(d_vec));

    delete[] d_data;

    return act_data_size;
}


void save_data_int(const char *data_path, long *int_data, long data_size) {
    FILE *fp = NULL;

    if (NULL == (fp = fopen(data_path, "wb"))) {
        printf("%s cannot be created.\n", data_path);
        exit(1);
    }

    fwrite(&data_size, sizeof(long), 1, fp);

    long rc = fwrite(int_data, sizeof(long), data_size, fp);
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


void save_data_int_vec(const char *data_path, const longVec &int_data_vec) {
    long data_size = int_data_vec.size();
    long *int_data = new long[data_size];
    std::copy(int_data_vec.begin(), int_data_vec.end(), int_data);

    FILE *fp = NULL;

    if (NULL == (fp = fopen(data_path, "wb"))) {
        printf("%s cannot be created.\n", data_path);
        exit(1);
    }

    fwrite(&data_size, sizeof(long), 1, fp);

    long rc = fwrite(int_data, sizeof(long), data_size, fp);
    assert(rc == data_size);
    fclose(fp);

    delete[] int_data;
}


void save_data_double(const char *data_path, double *d_data, long data_size) {
    FILE *fp = NULL;

    if (NULL == (fp = fopen(data_path, "wb"))) {
        printf("%s cannot be created.\n", data_path);
        exit(1);
    }

    fwrite(&data_size, sizeof(long), 1, fp);

    long rc = fwrite(d_data, sizeof(double), data_size, fp);
    assert(rc == data_size);
    fclose(fp);
}

void save_data_double_vec(const char *data_path, const doubleVec &d_vec) {
    long data_size = d_vec.size();
    double *d_data = new double[data_size];
    std::copy(d_vec.begin(), d_vec.end(), d_data);

    FILE *fp = NULL;

    if (NULL == (fp = fopen(data_path, "wb"))) {
        printf("%s cannot be created.\n", data_path);
        exit(1);
    }

    fwrite(&data_size, sizeof(long), 1, fp);

    long rc = fwrite(d_data, sizeof(double), data_size, fp);
    assert(rc == data_size);
    fclose(fp);

    delete[] d_data;
}

void save_data_double_mat(const char *data_path, const d_matrix &d_mat) {
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

long binarySearch(long *data, long l, long r, long x) {
    if (r >= l) {
        long mid = l + (r - l) / 2;

        // If the element is present at the middle
        // itself
        if (data[mid] == x)
            return mid;

        if (data[mid] > x)
            return binarySearch(data, l, mid - 1, x);

        return binarySearch(data, mid + 1, r, x);
    }
    return -1;
}


int path_status(const string &path_name) {
    struct stat info;
    if (stat(path_name.c_str(), &info) != 0) {
        return 0; // path not exists
    }
    else if (info.st_mode & S_IFDIR) {
        return 1; // path is a directory
    }
    else {
        return 2; // path is not a directory
    }
}

void create_dir(const string &dir) {
    string cmd = "mkdir -p " + dir;
    system(cmd.c_str());
//    if (__cplusplus < 201703L) {
//        string cmd = "mkdir -p " + dir;
//        system(cmd.c_str());
//    } else {
//        namespace fs = std::filesystem;
//        fs::create_directory(dir);
////        cout << "not supported in current." << endl;
//    }
}


void check(long *keys, int len) {
    for (int i = 0; i < len - 1; ++i) {
        if (keys[i] >= keys[i+1]) {
            cout << "i = " << i << ", keys[i] = " << keys[i] << ", keys[i+1] = " << keys[i+1] << endl;
        }
        assert(keys[i] < keys[i+1]);
    }
}


long perf_event_open(struct perf_event_attr *hw_event, pid_t pid,
                     int cpu, int group_fd, unsigned long flags) {
    int ret;

    ret = syscall(__NR_perf_event_open, hw_event, pid, cpu,
                  group_fd, flags);
    return ret;
}

void generate_data_descs(long *data, long n, const std::string &desc_path, const std::string &meta_info) {
    long min_diff = data[1] - data[0];
    long max_diff = min_diff;
    if (min_diff <= 0) {
        cout << "i = 1, min_diff = " << min_diff << ", data[i - 1] = " << data[0] << ", data[i] = " << data[1] << endl;
    }
    assert(min_diff > 0);
    for (long i = 2; i < n; ++i) {
        long diff = data[i] - data[i-1];
        if (min_diff <= 0) {
            cout << "i = " << i << ", diff = " << min_diff << ", data[i - 1] = " << data[i - 1] << ", data[i] = " << data[i] << endl;
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
         << "max_key = " << data[n-1] << endl
         << "min_diff = " << min_diff << endl
         << "max_diff = " << max_diff << endl;

    fout.close();
}

void read_lines(const std::string &path, std::vector<std::string> &lines) {
    ifstream fin(path, ios::in);
    string line;
    while(getline(fin, line)) {
        lines.push_back(line);
    }
    fin.close();
}