#include "diliNode.h"
#include <iostream>
using namespace std;

namespace dili_auxiliary {
    std::vector<fan2Leaf*> empty_fan2leaves;
    std::vector<diliNode*> empty_fan2nodes;
    std::vector<diliNode*> empty_nodes;
    long *retrain_keys = NULL;
    long *retrain_payloads = NULL;

    void init_insert_aux_vars() {
        if (retrain_keys) { delete[] retrain_keys; }
        if (retrain_payloads) { delete[] retrain_payloads; }
        retrain_keys = new long[LEAF_MAX_CAPACIY];
        retrain_payloads = new long[LEAF_MAX_CAPACIY];
    }

    void free_insert_aux_vars() {
        if (retrain_keys) {
            delete[] retrain_keys;
            retrain_keys = NULL;
        }
        if (retrain_payloads) {
            delete[] retrain_payloads;
            retrain_payloads = NULL;
        }
    }
}

