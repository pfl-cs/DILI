#ifndef DILI_BINARY_SEARCH_H
#define DILI_BINARY_SEARCH_H


// return i with data[i] <= x < data[i + 1]
template <class T, class K>
inline int binary_search_type_1(T *data, const K &x, int l, int r)  {
    while (l < r) {
        int mid = (l + r) >> 1;
        if (data[mid] <= x) {
            l = mid + 1;
        } else {
            r = mid;
        }
    }

    return l - 1;
}

#endif //DILI_BINARY_SEARCH_H
