#include "global/utils.h"
#include "global/global.h"
#include <iostream>
#include <cmath>
#include <cassert>
#include <algorithm>
#include <functional>
#include <utility>
#include <random>
#include <set>
#include <filesystem>
#include "dili/DILI.h"
using namespace std;

#define DATA_SIZE 10000000
#define INSERT_SIZE 1000000
#define DELETE_SIZE 1000000
#define QUERY_SIZE 1000000
KEY_TYPE *keys = NULL; // #define KEY_TYPE long
PAYLOAD_TYPE *payloads = NULL; // #define PAYLOAD_TYPE long
KEY_TYPE *keys_to_insert = NULL;
PAYLOAD_TYPE *payloads_to_insert = NULL;
KEY_TYPE *keys_to_delete = NULL;

pair<KEY_TYPE, PAYLOAD_TYPE> *queries;

void data_alloc() {
    keys = new KEY_TYPE[DATA_SIZE + INSERT_SIZE + DELETE_SIZE];
    payloads = new PAYLOAD_TYPE[DATA_SIZE + INSERT_SIZE];
    queries = new pair<KEY_TYPE, PAYLOAD_TYPE>[QUERY_SIZE];
    keys_to_insert = keys + DATA_SIZE;
    payloads_to_insert = payloads + DATA_SIZE;
    keys_to_delete = keys + DATA_SIZE + 1 + INSERT_SIZE;
}

void data_free() {
    if (keys) { delete[] keys; keys = NULL; }
    if (payloads) { delete[] payloads; payloads = NULL; }
    keys_to_insert = NULL; payloads_to_insert = NULL; keys_to_delete = NULL;
    if (queries) { delete[] queries; queries = NULL; }
}

bool data_exists(const string &data_dir) {
    auto keys_path = filesystem::path(data_dir) / "keys.dat";
    auto payloads_path = filesystem::path(data_dir) / "payloads.dat";
    auto queries_path = filesystem::path(data_dir) / "queries.dat";
    return file_exists(keys_path.c_str()) && file_exists(payloads_path.c_str()) && file_exists(queries_path.c_str());
}

void data_save(const string &data_dir) {
    auto keys_path = filesystem::path(data_dir) / "keys.dat";
    save_data_int(keys_path.c_str(), keys, DATA_SIZE + INSERT_SIZE + DELETE_SIZE);
    auto payloads_path = filesystem::path(data_dir) / "payloads.dat";
    save_data_int(payloads_path.c_str(), payloads, DATA_SIZE + INSERT_SIZE);
    auto queries_path = filesystem::path(data_dir) / "queries.dat";
    save_data_pair(queries_path.c_str(), queries, QUERY_SIZE);
}

void data_load(const string &data_dir) {
    auto keys_path = filesystem::path(data_dir) / "keys.dat";
    long size = load_data_int(keys_path.c_str(), keys);
    assert(size == DATA_SIZE + INSERT_SIZE + DELETE_SIZE);
    auto payloads_path = filesystem::path(data_dir) / "payloads.dat";
    size = load_data_int(payloads_path.c_str(), payloads);
    assert(size == DATA_SIZE + INSERT_SIZE + DELETE_SIZE);

    auto queries_path = filesystem::path(data_dir) / "queries.dat";
    size = load_data_pair(queries_path.c_str(), queries);
    assert(size == QUERY_SIZE);
}

void data_sampling() {
    std::random_device rd;
    std::mt19937 engine(rd());
    std::uniform_int_distribution<KEY_TYPE> dist(0, DATA_SIZE * 100);
    set<KEY_TYPE> key_set;
    while (key_set.size() < DATA_SIZE + INSERT_SIZE) {
        long key = dist(engine);
        key_set.insert(key);
    }
    std::copy(key_set.begin(), key_set.end(), keys);
    unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
    shuffle(keys, keys + DATA_SIZE + INSERT_SIZE, std::default_random_engine(seed));

    std::sort(keys, keys + DATA_SIZE);
    for (int i = 0; i < DATA_SIZE + INSERT_SIZE; ++i) {
        payloads[i] = i;
    }

    std::uniform_int_distribution<KEY_TYPE> query_and_delete_dist(0, DATA_SIZE-1);
    int n = QUERY_SIZE + DELETE_SIZE * 9 / 10;
    long *idxes = new long[n];
    set<KEY_TYPE> idx_set;
    while (idx_set.size() < n) {
        long idx = query_and_delete_dist(engine);
        idx_set.insert(idx);
    }
    std::copy(idx_set.begin(), idx_set.end(), idxes);
    shuffle(idxes, idxes + n, std::default_random_engine(seed));
    for (int i = 0; i < QUERY_SIZE; ++i) {
        queries[i].first = keys[idxes[i]];
        queries[i].second = payloads[idxes[i]];
    }
    for (int i = 0; i < DELETE_SIZE * 9 / 10; ++i) {
        keys_to_delete[i] = keys[idxes[QUERY_SIZE + i]];
    }


    std::uniform_int_distribution<KEY_TYPE> delete_second_part_dist(0, INSERT_SIZE - 1);
    idx_set.clear();
    while (idx_set.size() < DELETE_SIZE / 10) {
        long idx = delete_second_part_dist(engine);
        idx_set.insert(idx);
    }
    std::copy(idx_set.begin(), idx_set.end(), idxes);
    for (int i = 0; i < DELETE_SIZE / 10; ++i) {
        keys_to_delete[i + DELETE_SIZE * 9 / 10] = keys_to_insert[idxes[i]];
    }
    delete[] idxes;
}


int main(int argc, char *argv[]) {
    string data_dir = "data";
    int status = path_status(data_dir);
    assert(status != 2); // path is not a directory;
    if (status == 0) { // path does not exist
        create_dir(data_dir);
    }

    data_alloc();

    if (data_exists(data_dir)) {
        data_load(data_dir);
    } else {
        data_sampling();
        data_save(data_dir);
    }

    string mirror_dir = "models/buTree";
    status = path_status(mirror_dir);
    assert(status != 2);
    if (status == 0) {
        create_dir(mirror_dir);
    }


    DILI dili;
    dili.set_mirror_dir(mirror_dir);
    dili.bulk_load(keys, payloads, DATA_SIZE);
    cout << "0.here is OK." << endl;

    for (int i = 0; i < QUERY_SIZE; ++i) {
        PAYLOAD_TYPE pred = dili.search(queries[i].first);
        assert (pred == queries[i].second);
    }
    cout << "1.here is OK." << endl;

    for (int i = 0; i < INSERT_SIZE; ++i) {
        bool succeed = dili.insert(keys_to_insert[i], payloads_to_insert[i]);
        assert(succeed);
    }
    cout << "2.here is OK." << endl;

    for (int i = 0; i < QUERY_SIZE; ++i) {
        PAYLOAD_TYPE pred = dili.search(queries[i].first);
        assert (pred == queries[i].second);
    }
    cout << "3.here is OK." << endl;

    for (int i = 0; i < INSERT_SIZE; ++i) {
        bool succeed = dili.erase(keys_to_delete[i]);
        assert(succeed);
    }
    cout << "4.here is OK." << endl;
    for (int i = 0; i < QUERY_SIZE; ++i) {
        PAYLOAD_TYPE pred = dili.search(queries[i].first);
        assert (pred == queries[i].second);
    }
    cout << "5.here is OK." << endl;

    data_free();
    dili.clear();
    return 0;
}
