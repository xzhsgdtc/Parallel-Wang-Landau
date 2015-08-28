/* 
 * File:   Random.cpp
 * Author: jerry
 * 
 * Created on August 20, 2012, 1:52 PM
 */

#include "Random.h"
#include <cstdlib>
#include <cstdio>


Random::Random(const long s){
    seed = s;
    rng = gsl_rng_alloc(gsl_rng_ranlxd1);
    gsl_rng_set(rng, seed);
}

Random::Random(const Random& r){
    seed = r.seed;
    rng = gsl_rng_clone(r.rng);
}

Random& Random::operator=(const Random& orig){
    if(this == &orig) return *this;
    
    if(rng != NULL) gsl_rng_memcpy(rng, orig.rng);
    else rng = gsl_rng_clone(orig. rng);

    return *this;
}

Random::~Random(){
    gsl_rng_free(rng);
}

// 0-size
long Random::nextLong(long size){
    if(size == 0) return size;
    return gsl_rng_uniform_int(rng, size);
}

long Random::nextLong(long from , long to){
    long size = to-from;
    if(size == 0) return from;
    return from + gsl_rng_uniform_int(rng, size);
}

double Random::nextDouble(double from, double to){
    double size = to-from;
    return from + gsl_rng_uniform(rng)*size;
}

void Random::backup(std::string bkfname){
    FILE* f;
    f = fopen(bkfname.c_str(), "wb");
    gsl_rng_fwrite(f, rng);
    fclose(f);
}

void Random::recover(std::string bkfname){
    FILE* f;
    f = fopen(bkfname.c_str(), "rb");
    gsl_rng_fread(f, rng);
    fclose(f);
}
