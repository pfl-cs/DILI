#include <cstdio>
#include <cstdlib>
#include <cstring>

#ifndef DILI_FAN2LEAF_H
#define DILI_FAN2LEAF_H


struct fan2Leaf {
    long k1;
    long p1;
    long k2;
    long p2;
    fan2Leaf() {}
    fan2Leaf(const long &_k1, const long &_p1, const long &_k2, const long &_p2): k1(_k1), p1(_p1), k2(_k2), p2(_p2) {}
    inline void set(const long &_k1, const long &_p1, const long &_k2, const long &_p2) {
        k1 = _k1;
        p1 = _p1;
        k2 = _k2;
        p2 = _p2;
    }

    void save(FILE *fp) {
        fwrite(&k1, sizeof(long),1, fp);
        fwrite(&p1, sizeof(long),1, fp);
        fwrite(&k2, sizeof(long),1, fp);
        fwrite(&p2, sizeof(long),1, fp);
    }

    void load(FILE *fp) {
        fread(&k1, sizeof(long),1, fp);
        fread(&p1, sizeof(long),1, fp);
        fread(&k2, sizeof(long),1, fp);
        fread(&p2, sizeof(long),1, fp);
    }
};

#endif //DILI_FAN2LEAF_H
