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
        empty_fan2nodes.clear();
        empty_nodes.clear();
        empty_fan2leaves.clear();

        for (int i = 0; i < 1e7; ++i) {
//#ifndef ALLOW_FAN2_NODE
            empty_fan2leaves.push_back(new fan2Leaf());
//#else
            diliNode *_n2 = new diliNode(false);
            _n2->init(2);
            empty_fan2nodes.push_back(_n2);
//#endif
            diliNode *_n = new diliNode(false);
            _n->init(3);
            empty_nodes.push_back(_n);
        }
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
        empty_fan2leaves.clear();
        empty_nodes.clear();
        empty_fan2nodes.clear();
//        while (!empty_fan2leaves.empty()) {
//            empty_fan2leaves.pop_back();
//        }
//        while (!empty_nodes.empty()) {
//            empty_nodes.pop_back();
//        }
    }
}

