#include <cstdio>
#include <cstdlib>
#include <cstring>

#ifndef DILI_FAN2LEAF_H
#define DILI_FAN2LEAF_H


struct fan2Leaf {
    keyType k1;
    recordPtr p1;
    keyType k2;
    recordPtr p2;
    fan2Leaf() {}
    fan2Leaf(const keyType &_k1, const recordPtr &_p1, const keyType &_k2, const recordPtr &_p2): k1(_k1), p1(_p1), k2(_k2), p2(_p2) {}
    inline void set(const keyType &_k1, const recordPtr &_p1, const keyType &_k2, const recordPtr &_p2) {
        k1 = _k1;
        p1 = _p1;
        k2 = _k2;
        p2 = _p2;
    }

    void save(FILE *fp) {
        fwrite(&k1, sizeof(keyType),1, fp);
        fwrite(&p1, sizeof(recordPtr),1, fp);
        fwrite(&k2, sizeof(keyType),1, fp);
        fwrite(&p2, sizeof(recordPtr),1, fp);
    }

    void load(FILE *fp) {
        fread(&k1, sizeof(keyType),1, fp);
        fread(&p1, sizeof(recordPtr),1, fp);
        fread(&k2, sizeof(keyType),1, fp);
        fread(&p2, sizeof(recordPtr),1, fp);
    }
};

#endif //DILI_FAN2LEAF_H
