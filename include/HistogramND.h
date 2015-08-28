/* *********************************************************
 * :::::::::::::::::::::::::::::::::::::::::::::::::::::::::
 *  
 *   This file contains the template class for HistogramND
 *   (with helper class Bin). In histogram, the energy will 
 *   be divided longo bins and the energy range is [low, high)
 *   
 *   @File:   HistogramND.h/.cpp
 *   @Author: Guangjie Shi (Jerry)
 *
 *   @Version 1.0: Aug. 17, 2012
 *      Provide multiple dimension histogram 
 * 
 * :::::::::::::::::::::::::::::::::::::::::::::::::::::::::
 */

#ifndef HISTOGRAM_H
#define	HISTOGRAM_H
#include <iostream>
#include <fstream>
#include <sstream>

const double HISTOGRAM_DOUBLE_TOLERANCE = 1e-7;   // tolerance for double equaling

enum FlatnessType{
    DEFAULT_FLAT = 0,
    STRICT_FLAT = 1,
    ZHOU_PRE = 2
};

enum HistogramUpdateType{
    DEFAULT_UPDATE = 0
};


class HistogramND {
  
private:
    class Bin {
      public:        
        Bin();
        Bin(const Bin& one);
        ~Bin();
        void init(long d);
        void set(long d, double s, double e);
        Bin& operator=(Bin& one);
        
        double* min;    // lower bound of energy range
        double* max;    // upper bound of energy range
        double dos;     // density of states in logrithm
        long count;     // count of histogram
        long dim;
    };

public:

    /**
     * Constructor
     */
    HistogramND(long d=1);
    
    /**
     * Copy Constructor
     * @param one
     */
    HistogramND(const HistogramND& one);
    
    /**
     * Destructor
     */
    ~HistogramND();
    
    /**
     * Overwrite operator = for class HistogramND
     * @param one
     * @return 
     */
    HistogramND& operator=(const HistogramND& one);

    /**
     * Specify the range of given dimension
     * @param dim: which dimension
     * @param start: lower bound
     * @param end: upper bound
     * @param num: number of bins
     * @param flag: whether is integer type
     */
    void set(long dim, double start, double end, long num, bool flag=false);
     
    
    
    /**
     * build histogram into bins, using minE,maxE,size inside the class
     */
    void build();
        
    /**
     * To check whether the flatness criterion has been satisfied
     * @param fcriterion the flatness criterion
     * @param modfactor the current modification factor
     * @return true if histogram is flat, false otherwise
     */
    bool flat(const double& fcriterion, const double& modfactor);
        
    /**
     * Reset the count of histogram to zero
     */
    void reset();
        
    /**
     * Given energy, update corresponding bin, and return bin index
     * @param energy
     * @return index of bin that has been updated 
     */
    long update(const double* observables, double mf);
    
    /**
     * Given observables, find out the index of bin it/they corresponding to.
     * @param observables observables of the system
     */
    long bIndex(const double* observables) const;

    /**
     * Given the observables ( ususlly energy ), return the g(E)
     * @param observables 
     * @return density of state, or -1 if the observables are out of ranges
     */
    double dosAt(const double* observables);
    
    /**
     * Given a output stream, print the histogram
     * @param out
     */
    void print(std::ostream& out = std::cout);

    /**
     * Read in histogram from file. Note that this function only can read the file format
     * print out by the same program
     * @param in input file stream
     */
    void read(std::ifstream& fin);

    /**
     * Read in histogram from file. Note that this function only can read the file format
     * print out by the same program
     * @param fname input file name 
     */
    void read(std::string fname);
    
    /**
     * calculate the normalized density of states
     */
    void dosNorm();
    
    /**
     * This method is designed to write all informations of HistogramND
     * class to a binary file, in case of some accidence.
     */
    void backup(std::ofstream& fout);
    
    /**
     * This method is designed to recovery all information of HistogramND
     * class from a binary file.
     */
    void recover(std::ifstream& fin);
    
    
     /**
     *  Some methods help to access some 
     *  datas inside histogram 
     */
    inline long size(long d)     const   {return mBinsEachDim[d];}
    inline long size()           const   {return mNumberOfBins;}
    inline long dim()            const   {return mDim;}
    inline double min(long d=0)  const   {return mMinEachDim[d];}
    inline double max(long d=0)  const   {return mMaxEachDim[d];}
    inline double step(long d=0) const   {return mStepEachDim[d];}

    /**
     * Copy all histogram data into a array
     * @param logge array for copying histogram data
     */
    void dosCopyAll(double* logge);

    /**
     * Assign values to each bins, given a double array
     * @param logge array contain histogram data
     */
    void dosAssignAll(double* logge);

    void setOutput(std::ostream& out=std::cout);

    /**
     * Set the update type for histogram
     * @param ut HistogramUpdateType
     */
    inline void setUpdateType(HistogramUpdateType ut){
        mUpdateType = ut;
    }
    
    /**
     * Set the flat type for histogram
     * @param ut HistogramUpdateType
     */
    inline void setFlatType(FlatnessType ft){
        mFlatnessType = ft;
    }

private:
    

    FlatnessType mFlatnessType;
    HistogramUpdateType mUpdateType;
    

    std::ostream* mOut;

    // ::::::::::::::::::::::::::::::::::::::::
    // boolean
    bool* mDiscreteFlag;
    
    // ::::::::::::::::::::::::::::::::::::::::
    // long
    long mDim;
    long mNumberOfBins;    //total number of bins
    long* mBinsEachDim;        // number of bins for each dimension
    long* mWeightEachDim;
    long* mTempIndex;
    
    // ::::::::::::::::::::::::::::::::::::::::
    // double
    double* mStepEachDim;      // energy step for each bin
    double* mMinEachDim;       // minimum energy in histogram
    double* mMaxEachDim;       // maximum energy in histogram
   
    
    // ::::::::::::::::::::::::::::::::::::::::
    // self-defined
    Bin* mBins;
    
    // ::::::::::::::::::::::::::::::::::::::::
    // private methods
    
    /**
     * find the index
     * @param d: dimension
     * @param quantity: value of quantity in dimension-d
     * @return index of this quantity in this dimension
     */
    long findIndex(long d, double quantity) const;
};

#endif	/* HISTOGRAM_H */

