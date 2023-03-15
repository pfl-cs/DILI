#include "global.h"

//extern const double CACHE_MISS_LOSS = 1;
//extern const double EXPANDING_RATIO = -1;
//extern const int totalDataSize = 200000000;
//extern const long halfN = totalDataSize / 2;
//extern const long step = 100;
//extern const long start_idx = step / 2;
//extern const int one_in_n = 10;
//extern const long n_query_keys = halfN/step;
//extern const double RATIO = 10000;
//
//extern const double R1 = 25.0 / 130.0;
//extern const double R2 = 17.0 / 130.0;
//extern const double R3 = 5.0 / 130.0;


const double CACHE_MISS_LOSS = 1;
long totalDataSize = 200000000l;
long halfN = totalDataSize / 2;
const long n_query_keys = 1000000l;
long query_step = 100;
long query_start_idx = query_step / 2;

long num_adjust_stats = 0;


const double R1 = 25.0 / 130.0;
const double R2 = 17.0 / 130.0;
const double R3 = 5.0 / 130.0;


const int fanThreashold = 8192;
//const int minFan = 2;
//const int LEAF_MAX_CAPACIY = 8192;


//const int one_in_n = 10;
//const double RATIO = 10000;
//const int Delta = 1024;
//const int nThreads = 31;
//const int maxFan = 2048;
//const int minFanforSplit = 16;


double RHO = 0.1;
int buMinFan = 16;
double max_expanding_ratio = 6;
double retrain_threshold = 2;

