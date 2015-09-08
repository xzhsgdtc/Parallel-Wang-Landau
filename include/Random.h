/* *********************************************************
 * :::::::::::::::::::::::::::::::::::::::::::::::::::::::::
 *
 *  This file contains the definition of class Random.
 *  Which provide convient access to randome 
 *  number generator
 *
 *   @File:   Random.h/.cpp
 *   @Author: Guangjie Shi (Jerry)
 *
 *   @Version 1.0: Aug. 20, 2012
 *      Define random number generator as global function
 *
 *   @Version 2.0: Oct. 31, 2013
 *      Changed random number generator to become 
 *      a class instead of global functions
 *
 *   @Version 3.0: May  28, 2014
 *      Use random number generator from gsl (gsl_rng_ranlxd1) 
 *      instead of marsaglia
 *
 * :::::::::::::::::::::::::::::::::::::::::::::::::::::::::
 */


#ifndef RANDOM_H
#define	RANDOM_H

#include <ctime>
#include <iostream>
#include <fstream>
#include <gsl/gsl_rng.h>

/**
 * @class Random @file Random.h/.cpp
 * 
 * Random number generator
 *
 */

class Random {
    
  public:
    /**
     * Constructor
     */
    Random(const long s = 1234);

    /**
     * Copy Constructor
     */
    Random(const Random& r);

    /**
     * Destructor
     */
    virtual ~Random();

    /**
     * overload operator = for Random class
     */
    Random& operator=(const Random& r);

    /**
     * @return uniformly distributed double value between [from,to), 
     *         default range is [0,1)
     */
    double nextDouble(double from=0, double to=1);
    
    /**
     * @param to upper limit
     * @return uniformly distributed long value between [0,to);
     */
    long nextLong(long to=100);

    /**
     * @param from lower limit
     * @param to upper limit
     * @return uniformly distributed long value between [from,to);
     */
    long nextLong(long from , long to);

    /**
     * backup to binary file
     */
    void backup(std::string bkfname);

    /**
     * recover from binary file
     */
    void recover(std::string bkfname);

  private:
    long seed;
    gsl_rng* rng;

};

#endif	/* RANDOM_H */

