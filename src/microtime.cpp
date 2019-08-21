#include "lu.h"
#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#define MICRO_IN_SEC 1000000.00

double microtime(){
        int tv_sec,tv_usec;
        double time;
        struct timeval tv;
//         struct timezone tz;
        gettimeofday(&tv,NULL);

        return tv.tv_sec+tv.tv_usec/MICRO_IN_SEC;
}

