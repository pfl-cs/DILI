#include <iostream>
#include <vector>
#include <memory>

#ifndef GLOBAL_TYPEDEF_H
#define GLOBAL_TYPEDEF_H

#define KEY_TYPE long
#define PAYLOAD_TYPE long

typedef std::vector<int> intVec;
typedef std::vector<double> doubleVec;
typedef std::vector<long> longVec;
typedef std::vector<double> d_vector;
typedef std::vector< std::vector<double> > d_matrix;
typedef std::vector< std::vector<long> > l_matrix;
typedef std::vector< std::vector<int> > i_matrix;
typedef long recordPtr;
typedef long keyType;
typedef long recordContent;

typedef std::unique_ptr<int []> int32Array;
typedef std::unique_ptr<unsigned int []> uint32Array;
typedef std::unique_ptr<long []> int64Array;
typedef std::unique_ptr<unsigned long[]> uint64Array;
typedef std::unique_ptr<float []> floatArray;
typedef std::unique_ptr<double []> doubleArray;

typedef std::unique_ptr<keyType []> keyArray;
typedef std::unique_ptr<recordPtr []> recordPtrArray;

//typedef std::shared_ptr<keyType []> sharedKeyArray;
//typedef std::shared_ptr<recordPtr []> sharedRecordPtrArray;
//typedef std::shared_ptr<double []> sharedDoubleArray;

#endif //GLOBAL_TYPEDEF_H
