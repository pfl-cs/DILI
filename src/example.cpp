#include "global/global.h"
#include "utils/data_utils.h"
#include "utils/file_utils.h"
#include <iostream>
#include <cmath>
#include <cassert>
#include <algorithm>
#include <functional>
#include <utility>
#include <random>
#include <set>
#include <vector>
#include <filesystem>
#include "dili/DILI.h"
using namespace std;

struct Record {
    keyType key;
    recordContent content; // simulation

    Record(keyType _key, recordContent _content): key(_key), content(_content) {}
    void setNull() { content = -1; }
};

namespace database {
    vector<Record> underlying_data;
    recordPtr insert_to_database(const Record &record) {
        underlying_data.emplace_back(record);
        // use the index of the data array as the simulation of the record pointer
        return static_cast<recordPtr>(underlying_data.size() - 1);
    }

    void delete_from_database(const recordPtr &ptr) {
        underlying_data[ptr].setNull();
    }

    void clear() {
        underlying_data.clear();
    }
}

#define DATA_SIZE 10000000l
#define INSERT_SIZE 1000000l
#define DELETE_SIZE 1000000l
#define QUERY_SIZE 1000000l


typedef std::unique_ptr<Record []> RecordArray;
//vector< pair<keyType, recordPtr> > query_infos;
vector<keyType> query_keys;
vector<recordPtr> query_ptrs; // used to check if index return correct results, its size == query_keys.size()
vector<Record> records_to_insert;
vector<keyType> keys_to_delete;

vector< pair<keyType, recordPtr> > bulk_load_data;

void data_init() {
    std::random_device rd;
    std::mt19937 engine(rd());
    unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();

    std::uniform_int_distribution<keyType> dist(0, DATA_SIZE * 100);
    set<keyType> key_set;
    vector<keyType> keys;
    auto _n = DATA_SIZE + INSERT_SIZE;
    while (key_set.size() < _n) {
        keyType key = dist(engine);
        key_set.insert(key);
    }
    std::copy(key_set.begin(), key_set.end(), std::back_inserter(keys));

    shuffle(keys.begin(), keys.end(), std::default_random_engine(seed));

    data_utils::swap(keys, 0, data_utils::min_idx(keys));
    data_utils::swap(keys, DATA_SIZE - 1, data_utils::max_idx(keys));

    std::sort(keys.begin(), keys.begin() + DATA_SIZE);
    vector<Record> records;

    recordPtr last_ptr = -1;
    for (long i = 0; i < DATA_SIZE; ++i) {
        keyType key = keys[i];
        if (i > 0) {
            keyType last_key = keys[i - 1];
            assert(last_key < key);
        }
        recordContent content = abs(keys[i] + i); // simulation
        recordPtr ptr = database::insert_to_database(Record(keys[i], content));
        assert (ptr > last_ptr);
        last_ptr = ptr;
        bulk_load_data.emplace_back(make_pair(key, ptr));
//        records.emplace_back(Record(keys[i], content));
    }
    for (long i = DATA_SIZE; i < DATA_SIZE + INSERT_SIZE; ++i) {
        recordContent content = abs(keys[i] + i); // simulation
        records_to_insert.emplace_back(Record(keys[i], content));
    }
}

void test_data_sampling() {
    std::random_device rd;
    std::mt19937 engine(rd());
    unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();

    std::uniform_int_distribution<long> idx_dist(0, DATA_SIZE-1);
    size_t n = QUERY_SIZE + DELETE_SIZE * 9 / 10;
    set<long> idx_set;
    while (idx_set.size() < n) {
        long idx = idx_dist(engine);
        idx_set.insert(idx);
    }
    vector<long> idxes;
    std::copy(idx_set.begin(), idx_set.end(), std::back_inserter(idxes));
    shuffle(idxes.begin(), idxes.end(), std::default_random_engine(seed));

    for (long i = 0; i < QUERY_SIZE; ++i) {
        recordPtr ptr = static_cast<recordPtr>(idxes[i]);
        keyType query_key = database::underlying_data[ptr].key;
        query_keys.emplace_back(query_key);
        query_ptrs.emplace_back(ptr);
    }
    for (long i = 0; i < DELETE_SIZE * 9 / 10; ++i) {
        recordPtr ptr = static_cast<recordPtr>(idxes[QUERY_SIZE + i]);
        keyType key_to_delete = database::underlying_data[ptr].key;
        keys_to_delete.emplace_back(key_to_delete);
    }

    std::uniform_int_distribution<long> delete_idx_second_part_dist(0, INSERT_SIZE - 1);
    idx_set.clear();
    while (idx_set.size() < DELETE_SIZE / 10) {
        long idx = delete_idx_second_part_dist(engine);
        idx_set.insert(idx);
    }
    idxes.clear();
    std::copy(idx_set.begin(), idx_set.end(), std::back_inserter(idxes));

    for (long i = 0; i < DELETE_SIZE / 10; ++i) {
        keys_to_delete.emplace_back(records_to_insert[idxes[i]].key);
    }
}

bool data_exists(const string &data_dir) {
    auto records_path = filesystem::path(data_dir) / "records.dat";
    return file_utils::file_exists(records_path.c_str());
}

void data_save(const string &data_dir) {
    vector< pair<keyType, recordContent> > record_pairs;
    for (size_t i = 0; i < database::underlying_data.size(); ++i) {
        const Record &record = database::underlying_data[i];
        record_pairs.emplace_back(make_pair(record.key, record.content));
    }
    for (size_t i = 0; i < records_to_insert.size(); ++i) {
        const Record &record = records_to_insert[i];
        record_pairs.emplace_back(make_pair(record.key, record.content));
    }
    auto records_path = filesystem::path(data_dir) / "records.dat";
    data_utils::save_pair_vec_data(records_path.c_str(), record_pairs);
}

void data_load(const string &data_dir) {
    vector< pair<keyType, recordContent> > record_pairs;
    auto records_path = filesystem::path(data_dir) / "records.dat";
    data_utils::load_pair_vec_data(records_path.c_str(), record_pairs);
    database::clear();
    records_to_insert.clear();
    assert (record_pairs.size() == DATA_SIZE + INSERT_SIZE);
    for (long i = 0; i < DATA_SIZE; ++i) {
        const auto &p = record_pairs[i];
        bulk_load_data.emplace_back(make_pair(p.first, static_cast<recordPtr>(i)));
        database::insert_to_database(Record(p.first, p.second));
    }
    for (long i = DATA_SIZE; i < DATA_SIZE + INSERT_SIZE; ++i) {
        const auto &p = record_pairs[i];
        records_to_insert.emplace_back(Record(p.first, p.second));
    }
}


int main(int argc, char *argv[]) {
    cout << "Sampling data......." << endl;
    string data_dir = "data";
    int status = file_utils::path_status(data_dir);
    assert(status != 2); // path is not a directory;
    if (status == 0) { // path does not exist
        file_utils::detect_and_create_dir(data_dir);
    }
    if (data_exists(data_dir)) {
        data_load(data_dir);
    } else {
        data_init();
        data_save(data_dir);
    }
    test_data_sampling();

    string mirror_dir = "data/buTree";
    status = file_utils::path_status(mirror_dir);
    assert(status != 2);
    if (status == 0) {
        file_utils::detect_and_create_dir(mirror_dir);
    }

    DILI dili;
    dili.set_mirror_dir(mirror_dir);
    dili.bulk_load(bulk_load_data); // or uint64_t construction_time = dili.Build(bulk_load_data);

    cout << "Bulk loading test......";
    for (long i = 0; i < QUERY_SIZE; ++i) {
        recordPtr pred = dili.search(query_keys[i]);
        assert (pred == query_ptrs[i]);
    }
    cout << "finished." << endl;

    cout << "Insertion test......";
    for (long i = 0; i < INSERT_SIZE; ++i) {
        const Record &record = records_to_insert[i];
        recordPtr ptr = database::insert_to_database(records_to_insert[i]);
        assert(ptr > 0);
        bool succeed = dili.insert(make_pair(record.key, ptr));
        assert(succeed);
    }
    for (long i = 0; i < QUERY_SIZE; ++i) {
        recordPtr pred = dili.search(query_keys[i]);
        assert (pred == query_ptrs[i]);
    }

    cout << "finished." << endl;

    cout << "Deletion test......";
    for (int i = 0; i < DELETE_SIZE; ++i) {
        keyType delete_key = keys_to_delete[i];
        recordPtr ptr = dili.delete_key(delete_key);
        assert (ptr >= 0);
        database::delete_from_database(ptr);
    }

    for (long i = 0; i < QUERY_SIZE; ++i) {
        recordPtr pred = dili.search(query_keys[i]);
        assert (pred == query_ptrs[i]);
    }
    cout << "finished." << endl;

    dili.clear();
    return 0;
}
