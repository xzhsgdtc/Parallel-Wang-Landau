/* *********************************************************
 * :::::::::::::::::::::::::::::::::::::::::::::::::::::::::
 *  
 *   This file define the framework of Wang-Landau simulation,
 *   as well as parallel Wang-Landau simulation.
 *
 *   @File:   WLFrame.h
 *   @Author: Guangjie Shi (Jerry)
 *
 *   @Version 1.0: Aug. 16, 2012
 *      
 *   @Version 2.0: Nov. 11, 2013
 * 
 * :::::::::::::::::::::::::::::::::::::::::::::::::::::::::
 */


#ifndef WANGLANDAU_H
#define	WANGLANDAU_H

#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>
#include <cfloat>
#include <cstdlib>
#include <mpi.h>
#include "HistogramND.h"
#include "Random.h"
#include "Global.h"



const long   DEFAULT_RANDOM_SEED                = 7784;
const double DEFAULT_FLATNESS_CRITERION         = 0.8;
const double DEFAULT_INIT_MODIFICATION_FACTOR   = 1.0;
const double DEFAULT_MODIFICATION_DIVIDER       = 2.0;
const double DEFAULT_MODIFICATION_THRESHOLD     = 1e-8;

const std::string DEFAULT_CONF_PRINT_FILENAME_BASE = "_ConfPrint.mol2";
double alpha = 0.97;

class Replica {
public:
    long mID; // id for tracking purpose
    long* mLongData;
    double* mDoubleData;

    Replica():mID(0), mLongData(NULL), mDoubleData(NULL){}
    Replica(int n, int m){
        mLongData = new long[n];
        mDoubleData = new double[m];
    }
    ~Replica(){
        if(mLongData != NULL) delete [] mLongData;
        if(mDoubleData != NULL) delete [] mDoubleData;
    }
};




/** 
 * @class WLFrame Framework for Wang-Landau sampling.
 *
 * This is a template class, users need to define 
 * the model classes themselves and those model classes
 * should support certain methods to be called by 
 * class WLFrame<T>.
 *
 */

template <class ModelType>
class WLFrame {
    
public:
   
    /**
     * Constructor of WLFrame
     * @param filename input filename
     * @param id id of current process, in order to distingguish 
     *        between different processes in parallel simulation.
     */
    WLFrame(std::string filename, long id=0);

    /**
     * Copy constructor of WLFrame
     * @param orig another WLFrame object 
     */
    WLFrame(const WLFrame& orig);

    /**
     * Destructor
     */
    ~WLFrame();
    
    /**
     * Input routine for reading parameters
     * @param in input file stream object
     */
    void input(std::istream& in);

    /**
     * Input routine for reading parameters
     * @param filename input file name in string object
     */
    void input(const std::string& filename);

    /**
     * Initialized whole system. It will call init() routine in model class
     * and build up the histogram according to specified setting in input file
     */
    void init();
    
    /**
     * Call doMCMove routine in model class and return value accordingly.
     *
     * @return true if Monte Carlo move has been done, false if the move is rejected
     */
    void doMCMove();

    /**
     * Call undoMCMove routine in model class
     */
    void undoMCMove();

    /**
     * Based on backup observables and current observables,
     * decide the probability that previous Monte Carlo move be accepted.
     * 
     * @return probability of accpeting previous Monte Carlo move, 
     *         or negative number if current observable is out of range
     */
    double probAccept();
    
    /**
     * Call dosCopyAll(...) routine in HistogramND class,
     * to copy all the dos data into array logge
     *
     * @param logge a 1D double array (or pointer) to copy all log(dos) data in histogram
     *
     */
    void dosCopyAll(double* logge);

    /**
     * Call dosAssignAll(...) routine in HistogramND class,
     * to assign dos data from array logge into HistogramND object
     *
     * @param logge a 1D double array (or pointer) contain dos data in logrithm scale
     */
    void dosAssignAll(double* logge);
    
    /**
     * Calculate the dos ratio of two given sets of observables.
     * @param obs_up observables set in the nominator
     * @param obs_down observables set in the denominator
     * 
     * @return ratio of two given sets of observables, 
     *         or -1 if at least one of these observables are out of range
     */
    double dosRatio(const double* obs_up, const double* obs_down);
    double dosRatioInLog(const double* obs_up, const double* obs_down);
    
    /**
     * Collect data after simulation is satisfied the stopping criterion
     * ???????????????????????????????
     * to be implemented
     */
    void collectData();
    
    /**
     * Print out information of WLFrame
     * @param out
     */
    void print(std::ostream& out=std::cout, PrintType type = PRINT_WL_INFORMATION);
    
    
    void backup();
    bool recover();
    
    
    // :::::::::::::::::::::::::::::::::::::::::::: //
    //                                              //
    //          inline helper functions             //
    //                                              //
    // :::::::::::::::::::::::::::::::::::::::::::: //
    
    /**
     * Helper that help to generate error message
     *
     * @param func function in which error is encountered
     * @param msg error message
     */
    inline void errorMsg(std::string func, std::string msg){
        *mOut << "Error @" << func << ", " << msg << std::endl;
    }
    
    /**
     * Decided whether previous Monte Carlo move is accepted or not.
     *
     * @return true if accepted, otherwise false.
     */
    inline bool isMoveAccepted(){
        if(! mModel -> isMoveValid()) {
            return false;
        }
        else{
            return  (mRandom -> nextDouble() > probAccept() ? false : true) ;   
        }
    }

    /**
     * Set output stream
     * @param out output stream object, default value is std::cout
     */
    inline void setOutput(std::ostream& out = std::cout) { 
        mOut = &out; 
    }
    
    /**
     * Update the modification factor
     */
    inline void updateModificationFactor(){
        mCurrentIteration++;
        mCurrentMCMoves = 0;
        mModificationFactor /= mModificationDivider;
    }

    /**
     * Judge whether current iteration is satisfied the stopping criterion
     * @return true if ok to stop, 
     *         false if check histogram interval isn't satisfied or it is not ok to stop.
     */
    bool isIterationComplete(){
        return ( this -> isReadyCheckHistogram() ? this -> isHistogramFlat() : false ) ; 
    }
    
    /**
     * Judge whether whole simulation is ready to stop
     *
     * @return true if stopping criterion is satisfied, otherwise false
     */
    inline bool isSimulationComplete(){
        return ( mModificationFactor <  mModificationThreshold );
    }  
    
    /**
     * Judge whether the check histogram interval is satisfied
     *
     * @return true if satisfied, otherwise false
     */
    inline bool isReadyCheckHistogram(){
        return (mCurrentMCMoves != 0 && mCurrentMCMoves % mHistogramCheckInterval == 0);
    }

    /**
     * Judge whether the output interval is satisfied
     *
     * @return true if satisfied, otherwise false
     */
    inline bool isReadyPrint(){
        return (mCurrentMCMoves % mOutputInterval == 0);
    }

    /**
     * Judge whether the flatness criterion of histogram is satisfied
     *
     * @return true if satisfied, otherwise false
     */
    inline bool isHistogramFlat(){
        return mHistogram -> flat (mFlatnessCriterion, mModificationFactor);
    }

    /**
     * @return return current interation count
     */
    inline long iteration(){
        return mCurrentIteration;
    }

    /**
     * @return histogram check interval
     */
    inline long histogramCheckInterval(){
        return mHistogramCheckInterval;
    } 

    /**
     * @return current Monte Carlo moves count, 
     *         value should be clear when one iteration is complete
     */
    inline long currMoves(){
        return mCurrentMCMoves;
    }

    /**
     * @return total Monte Carlo moves count
     */
    inline long totalMoves(){
        return mTotalMCMoves;
    }

    /**
     * @return model size, particularly the total number of monomers in the system
     */
    inline long modelSize(){
        return mModel -> size(); 
    }
    
    /**
     * @return the total number of bins inside histogram
     */
    inline long histogramSize(){
        return mHistogram -> size();
    }

   
    /**
     * Given the energy, update histogram along with density of state
     * @param energy
     */
    inline void updateHistogram(){
        long i = mHistogram -> update(mModel->observables(), mModificationFactor);
        if(mConfPrintFlag){
            if( mConfPrintCount[i] < mConfPrintNumber ){
               print(mConfPrintStream, PRINT_MODEL_MOL2);
               mConfPrintCount[i]++ ;
            }
        }
    }

    /**
     * Reset histogram to zero, it should be called after completeness of an interation. 
     */
    inline void resetHistogram(){
        mHistogram -> reset();
    }
    
    /**
     * Adding prefix for given base string.
     * @param base
     * @return resultant string, prefix usually look like P00000{id}
     */
    inline std::string addPrefix(std::string base){
        return mPrefix + base;
    }

    inline std::string backupFilename(){
        return mPrefix+"PWLSim.BKUP";
    }

    inline double randomDouble(double from=0, double to=1){
        return mRandom -> nextDouble(from, to);
    }

    inline long randomLong(double from, double to){
        return mRandom -> nextLong(from, to);
    }

private:
    
    /**
     * Modify model to satisfy the energy range
     * @param mMinLimitsE
     * @param mMaxLimitsE
     */
    bool modifyModel();
    
    MPI_Datatype ReplicaType;
    std::ostream* mOut;
    std::ofstream mConfPrintStream;
    HistogramND* mHistogram;
    ModelType* mModel;
    Random* mRandom;
    Replica* mReplica;
    FlatnessType mHistogramFlatType;
    HistogramUpdateType mHistogramUpdateType;

    bool* mHistogramType;    // histogram type for every dimension
    long* mHistogramSize;    // histogram size for every dimension
    long* mConfPrintCount;  // count down the number of configurations that need to be print
    double* mMinLimits;    // mMinLimitsimum for every dimension
    double* mMaxLimits;    // mMaxLimitsimum for every
    double* mGlobalMinLimits;   // global mMinLimitsimum
    double* mGlobalMaxLimits;   // global mMaxLimitsimum


    bool mParallelSimFlag;  // flag represents whether it is parallel simulation
    bool mRecoverFlag;      // flag represents whether it is at recover phase
    bool mConfPrintFlag;    // configuration print flag, print specified number of structures in every bin of histogram
    long mConfPrintNumber;  // number of structures that will be print out for each bin of histogram
    long mRandomSeed;       // random number mRandomSeed
    long mID;               // my procees id in parallel simulation
    long mReplicaExchangeProposeCount;
    long mReplicaExchangeAcceptCount;


    long mCurrentIteration;    // Wang-Landau iteration count
    long mHistogramCheckInterval; // histogram checking interval
    long mReplicaExchangeInterval; // replica exchange interval
    long mCompleteCheckInterval; // PWL complete check interval
    long mOutputInterval; // observables tracking interval
    long mBackupInterval; // back up interval
    long mHistogramDim;    // histogram dimension
    long mSweep;            // monte carlo sweep, usually equal to number of monomers
    unsigned long mCurrentMCMoves;  // current Monte Carlo moves
    unsigned long mTotalMCMoves;
    
    
    double mFlatnessCriterion;  // flatness criterion
    double mModificationFactor;   // modification factor
    double mModificationDivider;   // modification factor divider
    double mModificationThreshold; // final modification factor
    
    std::string mModelInputFile; // input file name of model
    std::string mPrefix;


   
    
    
    // :::::::::::::::::::::::::::::::::::::::::::: //
    //                                              //
    //      Parallel Wang-Landau routines           //
    //                                              //
    // :::::::::::::::::::::::::::::::::::::::::::: //
       
public:

    inline bool PWL_isReadyReplicaExchange(){
        return (mTotalMCMoves != 0 && mTotalMCMoves % mReplicaExchangeInterval == 0);
    }
    inline bool PWL_isReadyCheckComplete(){
        return (mTotalMCMoves != 0 && mTotalMCMoves % mCompleteCheckInterval == 0);
    }

    bool PWL_isHistogramFlat(MPI_Comm* mpi_intra_win_comm);
    bool PWL_isSimulationComplete();
    bool PWL_proposeReplicaExchange(int num_procs_per_win, int pid_inter_win_comm, int inter_win_comm_id, MPI_Comm* mpi_inter_win_comm);
    void PWL_dosMerge(MPI_Comm* mpi_intra_win_comm, int num_procs_per_win);
    void PWL_packupReplica(Replica& pack);
    void PWL_unpackReplica(Replica& pack);
    void PWL_replicaRegisterAndInit();
    
    /**
     *
     * @dim : dimension of windows
     * @size : number of windows in each direction
     * @id : window id of this process in each direction 
     * @overlap: overlap of each direction
     */
    void PWL_init(int dim, int* size, int* id, double* overlap);

    
};

//----------------------Below are implementation------------------------------

template<class ModelType>
WLFrame<ModelType>::WLFrame(std::string filename, long id){
    
    mID                                 = id;
    
    mRecoverFlag                        = false;
    mParallelSimFlag                    = false;
    mConfPrintFlag                      = false;
    
    mConfPrintCount                     = NULL;
    mRandom                             = NULL;
    mModel                              = NULL;
    mHistogram                          = NULL;
    mReplica                            = NULL;
    mGlobalMinLimits                    = NULL;
    mGlobalMaxLimits                    = NULL;
    mMinLimits                          = NULL;
    mMaxLimits                          = NULL;
    mHistogramSize                      = NULL;
    mHistogramType                      = NULL;
    setOutput(std::cout);

    mRandomSeed                         = DEFAULT_RANDOM_SEED;
    
    mFlatnessCriterion                  = DEFAULT_FLATNESS_CRITERION       ;
    mModificationFactor                 = DEFAULT_INIT_MODIFICATION_FACTOR ;
    mModificationDivider                = DEFAULT_MODIFICATION_DIVIDER     ;
    mModificationThreshold              = DEFAULT_MODIFICATION_THRESHOLD   ;

    mHistogramUpdateType                = DEFAULT_UPDATE;
    mHistogramFlatType                  = DEFAULT_FLAT;
    mHistogramDim                       = 0;
    mCurrentIteration                   = 0;
    mHistogramCheckInterval             = 0; 
    mReplicaExchangeInterval            = 0; 
    mOutputInterval                     = 0;
    mCompleteCheckInterval              = 0;
    mBackupInterval                     = 0;
    mCurrentMCMoves                     = 0;
    mTotalMCMoves                       = 0; 
    mReplicaExchangeProposeCount        = 0;
    mReplicaExchangeAcceptCount         = 0;
    mSweep                              = 0;
    mConfPrintNumber                    = 0;

    mPrefix                             = "P" + intToString(mID,7);
    input(filename);
}

template <class ModelType>
WLFrame<ModelType>::WLFrame(const WLFrame& orig) {

    setOutput(*orig.mOut);
    mID = orig.mID;
    mRecoverFlag = orig.mRecoverFlag;
    mParallelSimFlag = orig.mParallelSimFlag;
    mConfPrintFlag = orig.mConfPrintFlag;
    mConfPrintNumber = orig.mConfPrintNumber;
    mRandomSeed = orig.mRandomSeed;
    mFlatnessCriterion = orig.mFlatnessCriterion;
    mHistogramDim = orig.mHistogramDim;
    mHistogramUpdateType = orig.mHistogramUpdateType;
    mHistogramFlatType = orig.mHistogramFlatType;
    mModificationFactor = orig.mModificationFactor;
    mModificationThreshold = orig.mModificationThreshold;
    mModificationDivider = orig.mModificationDivider;
    mSweep = orig.mSweep;
    mCurrentIteration = orig.mCurrentIteration;
    mCurrentMCMoves = orig.mCurrentMCMoves;
    mTotalMCMoves = orig.mTotalMCMoves;
    mBackupInterval = orig.mBackupInterval;
    mOutputInterval = orig.mOutputInterval;
    mCompleteCheckInterval = orig.mCompleteCheckInterval;
    mReplicaExchangeInterval = orig.mReplicaExchangeInterval;
    mHistogramCheckInterval = orig.mHistogramCheckInterval;
    mReplicaExchangeProposeCount = orig.mReplicaExchangeProposeCount;
    mReplicaExchangeAcceptCount = orig.mReplicaExchangeAcceptCount;
    mModelInputFile = orig.mModelInputFile;
    mPrefix = orig.mPrefix; 

    mGlobalMinLimits = new double[mHistogramDim];
    mGlobalMaxLimits = new double[mHistogramDim];
    mMinLimits = new double[mHistogramDim];
    mMaxLimits = new double[mHistogramDim];
    mHistogramSize = new long[mHistogramDim];
    mHistogramType = new bool[mHistogramDim];

    for(int i=0; i<mHistogramDim; ++i){
        if(orig.mGlobalMinLimits != NULL)   mGlobalMinLimits[i] = orig.mGlobalMinLimits[i];
        if(orig.mGlobalMaxLimits != NULL)   mGlobalMaxLimits[i] = orig.mGlobalMaxLimits[i];
        if(orig.mMinLimits       != NULL)   mMinLimits[i] = orig.mMinLimits[i];
        if(orig.mMaxLimits       != NULL)   mMaxLimits[i] = orig.mMaxLimits[i];
        if(orig.mHistogramSize   != NULL)   mHistogramSize[i] = orig.mHistogramSize[i];
        if(orig.mHistogramType   != NULL)   mHistogramType[i] = orig.mHistogramType[i];
    }

    mRandom = new Random();
    *mRandom = *orig.mRandom;
    
    mHistogram = new HistogramND(mHistogramDim);
    *mHistogram = *orig.mHistogram;
    
    mModel = new ModelType();
    *mModel = *orig.mModel;
    
    mReplica = NULL;
    if(mParallelSimFlag){ 
        if( orig.mReplica != NULL){
            int num_long = 2 + 2* mModel -> size();
            int num_double = 1 +  mModel -> dim() * mModel -> size();
            PWL_replicaRegisterAndInit();
            for(long i=0; i<num_long; ++i)
                mReplica -> mLongData[i] = orig.mReplica -> mLongData[i];
            for(long i=0; i<num_double; ++i)
                mReplica -> mDoubleData[i] = orig.mReplica -> mDoubleData[i];
        }
    }

    if (orig.mConfPrintFlag){
        long int size = mHistogram->size();
        mConfPrintCount = new long[size];
        for(int i=0; i < size; ++i){
            mConfPrintCount[i] = orig.mConfPrintCount[i];
        }
        getOutputStream(mConfPrintStream, mPrefix + DEFAULT_CONF_PRINT_FILENAME_BASE, true);
    }else{
        mConfPrintCount = NULL;
        mConfPrintStream = NULL;
    }
        
}

template <class ModelType>
WLFrame<ModelType>::~WLFrame() {
    if(  mRandom             != NULL )  delete mRandom         ;
    if(  mModel              != NULL )  delete mModel          ;
    if(  mHistogram          != NULL )  delete mHistogram      ;
    if(  mReplica            != NULL )  delete mReplica        ;
    if(  mGlobalMinLimits    != NULL )  delete []   mGlobalMinLimits;
    if(  mGlobalMaxLimits    != NULL )  delete []   mGlobalMaxLimits;
    if(  mMinLimits          != NULL )  delete []   mMinLimits      ;
    if(  mMaxLimits          != NULL )  delete []   mMaxLimits      ;
    if(  mHistogramSize      != NULL )  delete []   mHistogramSize  ;
    if(  mHistogramType      != NULL )  delete []   mHistogramType  ;
    if ( mConfPrintCount     != NULL )  delete []   mConfPrintCount ;
    if ( mConfPrintStream.is_open()  )  mConfPrintStream.close();
}

template<class ModelType>
bool WLFrame<ModelType>::modifyModel(){
    // using simulated annealing 

    long count(0);
    long coolInterval = 10*mSweep;
    double T(10), thresholdT(0.01);
    double de(0), e(0);

    while( !mModel -> satisfyLimits() && T > thresholdT){
        mModel->doMCMove();
        de = mModel -> energyDifference();
        e = mModel -> energy();
        if(isnan(e)){
            mModel -> undoMCMove();
        }else{
            if ( mRandom -> nextDouble() > exp( -de / T ) ){
                mModel -> undoMCMove();
            }
        }
        ++count;
        if(count % coolInterval == 0) {
            if (T > thresholdT){
                T *= alpha;
                *mOut << "          * T = " << T << "; "; mModel -> print(*mOut, PRINT_MODEL_OBSERVABLES_ONLY);            
            }

        }
    }

    bool result(true);
    if (mModel -> satisfyLimits()){
        result = true;
        mModel -> resetStatistics();
        *mOut  << "          * finish modifying, observables: "; mModel -> print(*mOut, PRINT_MODEL_OBSERVABLES_ONLY); 
    } else{
        result = false;
        *mOut  << "          * modify failed, need restart, current observables: "; mModel -> print(*mOut, PRINT_MODEL_OBSERVABLES_ONLY); 
    }
    return result;
}


template <class ModelType>
void WLFrame<ModelType>::doMCMove(){
    // check whether backup criterion satisfied
    if( mTotalMCMoves !=0 && (mTotalMCMoves % mBackupInterval ) == 0) { 
        if (! mRecoverFlag)
                backup();
        else 
            mRecoverFlag = false;
    }

    // do MC moves
    mCurrentMCMoves++;
    mTotalMCMoves++;
    
    mModel -> doMCMove();
}

template <class ModelType>
void WLFrame<ModelType>::undoMCMove(){
    mModel -> undoMCMove();
}


template <class ModelType>
double WLFrame<ModelType>::probAccept(){
    double geold, genew, poss;
    geold = mHistogram->dosAt(mModel->prevObservables());
    genew = mHistogram->dosAt(mModel->observables());
    if(genew < 0) return -1;
    if(geold < 0){
        errorMsg("probAccept", "backup observables cannot excess the energy limit !");
        exit(1);
    }
    poss = (mModel->cFactor())*exp(geold - genew);
    return (1 > poss ? poss : 1);
}

template <class ModelType>
void WLFrame<ModelType>::dosCopyAll(double* logge){
    mHistogram->dosCopyAll(logge);
}

template <class ModelType>
void WLFrame<ModelType>::dosAssignAll(double* logge){
    mHistogram->dosAssignAll(logge);
}

template <class ModelType>
double WLFrame<ModelType>::dosRatio(const double* obs_up, const double* obs_down){
    double ge_up, ge_down, poss;
    ge_up = mHistogram->dosAt(obs_up);
    ge_down = mHistogram->dosAt(obs_down);
    
    if(ge_up < 0 || ge_down < 0){
        *mOut << "Error: observables are out of the range of histogram !" << std::endl;
        exit(1);
    }
    poss = exp(ge_up - ge_down);
    return poss;
}

template <class ModelType>
double WLFrame<ModelType>::dosRatioInLog(const double* obs_up, const double* obs_down){
    double ge_up, ge_down, log_poss;
    ge_up = mHistogram->dosAt(obs_up);
    ge_down = mHistogram->dosAt(obs_down);
    
    if(ge_up < 0 || ge_down < 0){
        *mOut << "Error: observables are out of the range of histogram !" << std::endl;
        exit(1);
    }
    log_poss = ge_up - ge_down;
    return log_poss;
}



template <class ModelType>
void WLFrame<ModelType>::collectData(){
    // to be implemented
    
    std::ofstream fout[2];
    getOutputStream(fout[0],mPrefix+"_dos.dat");
    getOutputStream(fout[1],mPrefix+"_norm_dos.dat");
    mHistogram-> print(fout[0]);
    mHistogram-> dosNorm();
    mHistogram-> print(fout[1]);
    fout[0].close();
    fout[1].close();
    
}

template<class ModelType>
void WLFrame< ModelType >::input(std::istream& in){
    std::string smMinLimits, smMaxLimits, smHistogramSize, smHistogramType; // temporary storage of histogram information
    std::string line;
    while(getline(in, line)){
        if(line[0] == '/') continue;
        std::istringstream ins(line);
        std::string word;
        ins >> word;
        if(word.compare("RANDOMSEED")                       == 0){
            ins >> mRandomSeed;
        } else if(word.compare("FLATNESSCRITERION")         == 0){
            ins >> mFlatnessCriterion;
        } else if(word.compare("MODIFICATIONFACTOR")        == 0){
            ins >> mModificationFactor;
        } else if(word.compare("MODIFICATIONDIVIDER")       == 0){
            ins >> mModificationDivider;
        } else if(word.compare("MODIFICATIONTHRESHOLD")     == 0){
            ins >> mModificationThreshold;
        } else if(word.compare("MINIMUMLIMITS")             == 0){
            smMinLimits = line;
        } else if(word.compare("MAXIMUMLIMITS")             == 0){
            smMaxLimits = line;
        } else if(word.compare("HISTOGRAMTYPES")            == 0){
            smHistogramType = line;
        } else if(word.compare("HISTOGRAMUPDATETYPE")       == 0){
            long temp;
            ins >> temp;
            mHistogramUpdateType = (HistogramUpdateType)temp;
        } else if(word.compare("HISTOGRAMFLATTYPE")         == 0){
            long temp;
           ins >> temp;
           mHistogramFlatType = (FlatnessType) temp;
        } else if(word.compare("HISTOGRAMDIMENSION")        == 0){
            ins >> mHistogramDim;
        } else if(word.compare("HISTOGRAMSIZE")             == 0){
            smHistogramSize = line;
        } else if(word.compare("HISTOGRAMCHECKINTERVAL")    == 0){
            ins >> mHistogramCheckInterval;
        } else if(word.compare("MODELINPUTFILE")            == 0){
            ins >> mModelInputFile;
        } else if(word.compare("BACKUPINTERVAL")            == 0){
            ins >> mBackupInterval;
        } else if(word.compare("REPLICAEXCHANGEINTERVAL")   == 0){
            ins >> mReplicaExchangeInterval;
        } else if(word.compare("PARALLELSIMULATIONFLAG")    == 0){
            ins >> mParallelSimFlag;
        } else if(word.compare("COMPLETECHECKINTERVAL")     ==0){
            ins >> mCompleteCheckInterval;
        } else if(word.compare("OUTPUTINTERVAL")            == 0){
            ins >> mOutputInterval;
        } else if(word.compare("CONFPRINTFLAG")            == 0){
            ins >> mConfPrintFlag >> mConfPrintNumber;
        }
    }
    
       
    // ::::::::::::::::::::::::::::::::::::::
    // get parameters for building histogram
    std::istringstream sin[4];
    sin[0].str(smMinLimits);
    sin[1].str(smMaxLimits);
    sin[2].str(smHistogramSize);
    sin[3].str(smHistogramType);
    std::string waste;
    for(long i=0; i<4; ++i)
        sin[i] >> waste;
    mGlobalMinLimits = new double[mHistogramDim];
    mGlobalMaxLimits = new double[mHistogramDim];
    mMinLimits = new double[mHistogramDim];
    mMaxLimits = new double[mHistogramDim];
    mHistogramSize = new long[mHistogramDim];
    mHistogramType = new bool[mHistogramDim];
    for(long i=0; i<mHistogramDim; ++i){
        sin[0] >> mMinLimits[i];
        sin[1] >> mMaxLimits[i];
        sin[2] >> mHistogramSize[i];
        sin[3] >> mHistogramType[i];
        
        mGlobalMinLimits[i] = mMinLimits[i];
        mGlobalMaxLimits[i] = mMaxLimits[i];
    }
}

template<class ModelType>
void WLFrame<ModelType>::input(const std::string& filename){
    std::ifstream fin;
    getInputStream(fin,filename);
    input(fin);
    fin.close();
}


template<class ModelType>
void WLFrame<ModelType>::init(){
    
    if ( recover() ) return;
     
    *mOut << "\n"
          << "\n# Wang-Landau Simulation Normal Start: \n" 
          << "  >> initialize random number generator ... " << std::endl;
    // Initialize Random Number Generator
    Random* temp_ran = new Random(mRandomSeed);
    for(int i=0; i< mID; ++i)   temp_ran -> nextLong(0,10000000);
    mRandomSeed = temp_ran -> nextLong(0,10000000);
    delete temp_ran;
    
    mRandom     = new Random(mRandomSeed);
    
    *mOut << "     * random seed: " << std::setw(10) << mRandomSeed << std::endl; 
     
    *mOut << "  >> build up histogram ... " << std::endl;
    // Build Up the Histogram
    mHistogram = new HistogramND(mHistogramDim);
    mHistogram->setOutput(*mOut);
    mHistogram->setUpdateType(mHistogramUpdateType);
    mHistogram->setFlatType(mHistogramFlatType);
     
    for (int i = 0; i < mHistogramDim; ++i) {
        mHistogram -> set(i, mMinLimits[i], mMaxLimits[i], mHistogramSize[i], mHistogramType[i]);
    }
    mHistogram -> build();
    
    if( mConfPrintFlag ){
        *mOut << "  >> configuration print flag is ON ! \n"
              << "       * number of configurations that is going to print for each bin is " << mConfPrintNumber << std::endl;
        long size = mHistogram -> size();
        mConfPrintCount = new long[size];
        for(int i=0; i<size; ++i) mConfPrintCount[i] = 0;
        getOutputStream(mConfPrintStream, mPrefix + DEFAULT_CONF_PRINT_FILENAME_BASE, true);
    } else{
        *mOut << "  >> configuration print flag is OFF !" << std::endl;
    }
        
    *mOut << "  >> initialize model ... " << std::endl;
    // Initialize Model
    mModel = new ModelType();
    mModel -> input(mModelInputFile);
    mModel -> setOutput(*mOut);
    mModel -> setConstraint(mHistogram);
    mModel -> setRandom(*mRandom);
    mModel -> init(addPrefix("_init_conf.mol2"));

    // Processing Other Parameters
    mSweep = mModel -> size();
    mHistogramCheckInterval             *= mSweep;
    mReplicaExchangeInterval            *= mSweep; 
    mOutputInterval                     *= mSweep;
    mCompleteCheckInterval              *= mSweep;
    mBackupInterval                     *= mSweep;

    mCurrentIteration                   = 1;
    mCurrentMCMoves                     = 0;
    mTotalMCMoves                       = 0; 
    mReplicaExchangeProposeCount        = 0;
    mReplicaExchangeAcceptCount         = 0;
    
    // Modify model to fit constraints 
    long count_fail = 0;
    *mOut << "  >> modifying model to fit constraints ... " << std::endl;
    while( !modifyModel() ){
        //modify failed, need to re-initialize model and random number generator
        
        count_fail ++;
        if (count_fail % 10 == 0 && alpha < 0.999){
           alpha *= 1.01; 
           if (alpha > 0.999){
               alpha = 0.999;
           }
           *mOut << "     * increase the factor alpha to " << std::setprecision(10) << alpha << std::endl;
        }

        mRandomSeed = mRandom->nextLong(0,10000000);
        
        delete mRandom;
        delete mModel;
        *mOut << "  >> re-initialize random number generator ... \n" 
              << "     * re-generated random seed: " << std::setw(10) << mRandomSeed << std::endl; 
        mRandom     = new Random(mRandomSeed);
    
        *mOut << "  >> re-initialize model ... " << std::endl;
        mModel = new ModelType();
        mModel -> input(mModelInputFile);
        mModel -> setOutput(*mOut);
        mModel -> setConstraint(mHistogram);
        mModel -> setRandom(*mRandom);
        mModel -> init(addPrefix("_init_conf.mol2"));
    }
     
        
    // Print Initial Information
    const long NOF = 4; // NOF = Number Of Files
    std::string files[NOF]={  mPrefix+"_model_info.dat",
                              mPrefix+"_cell_info.dat",
                              mPrefix+"_conf_info.mol2",
                              mPrefix+"_histogram_info.dat"};
                          
    std::ofstream fout[NOF];
     
    for(long i=0; i<NOF; ++i){
        getOutputStream(fout[i],files[i]);
    }
     
    print(fout[0], PRINT_WL_INFORMATION);
    print(fout[1], PRINT_MODEL_CELL);
    print(fout[2], PRINT_MODEL_MOL2);
    print(fout[3], PRINT_WL_HISTOGRAM);
     
     
    for(long i=0; i<NOF; ++i)  fout[i].close();

    *mOut << "  >> finish initialization !" << std::endl;
}


template<class ModelType>
void WLFrame<ModelType>::print(std::ostream& out, PrintType type){
    switch(type){

        case PRINT_WL_INFORMATION:
            {
                out << "### This file contains the information of Wang-Landau sampling.\n"
                    << "### Which should contain the parameters of Wang-Landau method,\n"
                    << "### as well as model and histogram informations.\n\n"
                    << "***********************************************************\n" 
                    << "*                                                         *\n" 
                    << "*               WANG-LANDAU INFORMATION                   *\n" 
                    << "*                                                         *\n" 
                    << "***********************************************************\n"
                    << std::setw(30) << "MyID "                     << std::setw(20) << mID                          <<"\n" 
                    << std::setw(30) << "RandomSeed "               << std::setw(20) << mRandomSeed                  <<"\n" 
                    << std::setw(30) << "FlatnessCriterion "        << std::setw(20) << mFlatnessCriterion           <<"\n" 
                    << std::setw(30) << "ModificationFactor "       << std::setw(20) << mModificationFactor          <<"\n" 
                    << std::setw(30) << "ModificationDivider "      << std::setw(20) << mModificationDivider         <<"\n" 
                    << std::setw(30) << "ModificationThreshold "    << std::setw(20) << mModificationThreshold       <<"\n" 
                    << std::setw(30) << "His. Checking Interval "   << std::setw(20) << mHistogramCheckInterval   <<"\n" 
                    << std::setw(30) << "Output Interval "          << std::setw(20) << mOutputInterval              <<"\n" 
                    << std::setw(30) << "Replica Exchange Interval" << std::setw(20) << mReplicaExchangeInterval     <<"\n" 
                    << std::setw(30) << "Complete Check Interval "  << std::setw(20) << mCompleteCheckInterval       <<"\n" 
                    << std::setw(30) << "Backup Interval "          << std::setw(20) << mBackupInterval              <<"\n" 
                    << std::setw(30) << "Current MC Moves "         << std::setw(20) << mCurrentMCMoves              <<"\n" 
                    << std::setw(30) << "Total MC Moves "           << std::setw(20) << mTotalMCMoves                <<"\n"  
                    << std::endl;

                
                // print out histogram information
                for(long i=0; i<mHistogram->dim(); ++i){
                    out << std::setw(30) << "Dimension "  + intToString(i)   << " : "                               << "\n" 
                        << std::setw(30) << "Min. Limits"                    << std::setw(20)<< mHistogram->min(i) << "\n" 
                        << std::setw(30) << "Max. Limits"                    << std::setw(20)<< mHistogram->max(i) << "\n" 
                        << std::setw(30) << "Histogram Size"                 << std::setw(20)<< mHistogram->size(i)       << "\n";
                }

                // print out model information
                mModel -> print(out, PRINT_MODEL_INFORMATION);    
            }
            break;

        case PRINT_WL_HISTOGRAM:
            mHistogram -> print(out);
            break;

        case PRINT_WL_TRACKING_INFORMATION:
            {    
                out << std::setw(15) << mTotalMCMoves/mSweep   << " "
                    << std::setw(15) << mCurrentMCMoves/mSweep << " ";
                if(mParallelSimFlag)
                    out << std::setw(7) << mReplica -> mID << " ";
                out << std::setw(7) << mID << " "
                    << std::setw(10) << std::setprecision(8) << (mReplicaExchangeProposeCount == 0? 0 : (double) (mReplicaExchangeAcceptCount)/mReplicaExchangeProposeCount) << " ";
                mModel -> print(out, PRINT_MODEL_OBSERVABLES);
            }
            break;

        case PRINT_MODEL_INFORMATION:
        case PRINT_MODEL_CELL:
        case PRINT_MODEL_XYZ:
        case PRINT_MODEL_MOL2:
        case PRINT_MODEL_STATISTICS:
        case PRINT_MODEL_LIPIDWATER_ARRAY:
        case PRINT_MODEL_OBSERVABLES:
        case PRINT_MODEL_OBSERVABLES_ONLY:
        case PRINT_MAX_MOVE_DISTANCE:
            mModel -> print(out, type);
            break;

    }
}

template<class ModelType>
void WLFrame< ModelType>::backup(){
      
    std::string backupfile = backupFilename();
    std::ofstream fout(backupfile.c_str(), std::ios::binary);
    
    if(!fout.is_open()){
        errorMsg("backup","cannot open file: WLSim.BKUP!" );
        exit(1);
    }

    long sod = sizeof(double);
    long sol = sizeof(long);
    
    *mOut << "### system backup @iteration " << mCurrentIteration
          << "   @MCMoves  " << mCurrentMCMoves << " / " << mTotalMCMoves << std::endl;

    fout.write(reinterpret_cast<char*>(&mParallelSimFlag)               , sizeof(bool));
    fout.write(reinterpret_cast<char*>(&mRandomSeed)                    , sol);
    fout.write(reinterpret_cast<char*>(&mID)                            , sol);
    fout.write(reinterpret_cast<char*>(&mReplicaExchangeProposeCount)   , sol);
    fout.write(reinterpret_cast<char*>(&mReplicaExchangeAcceptCount)    , sol);
    fout.write(reinterpret_cast<char*>(&mCurrentIteration)              , sol);
    fout.write(reinterpret_cast<char*>(&mHistogramCheckInterval)        , sol);
    fout.write(reinterpret_cast<char*>(&mReplicaExchangeInterval)       , sol);
    fout.write(reinterpret_cast<char*>(&mCompleteCheckInterval)         , sol);
    fout.write(reinterpret_cast<char*>(&mOutputInterval)                , sol);
    fout.write(reinterpret_cast<char*>(&mBackupInterval)                , sol);
    fout.write(reinterpret_cast<char*>(&mHistogramDim)                  , sol);
    fout.write(reinterpret_cast<char*>(&mSweep)                         , sol);
    fout.write(reinterpret_cast<char*>(&mCurrentMCMoves)                , sizeof(unsigned long));
    fout.write(reinterpret_cast<char*>(&mTotalMCMoves)                  , sizeof(unsigned long));
    fout.write(reinterpret_cast<char*>(&mFlatnessCriterion)             , sod);
    fout.write(reinterpret_cast<char*>(&mModificationFactor)            , sod);
    fout.write(reinterpret_cast<char*>(&mModificationDivider)           , sod);
    fout.write(reinterpret_cast<char*>(&mModificationThreshold)         , sod);

    fout.write(reinterpret_cast<char*>(mHistogramType)                  , mHistogramDim*sizeof(bool));
    fout.write(reinterpret_cast<char*>(mHistogramSize)                  , mHistogramDim*sol);
    fout.write(reinterpret_cast<char*>(mMinLimits)                      , mHistogramDim*sod);
    fout.write(reinterpret_cast<char*>(mMaxLimits)                      , mHistogramDim*sod);
    fout.write(reinterpret_cast<char*>(mGlobalMinLimits)                , mHistogramDim*sod);
    fout.write(reinterpret_cast<char*>(mGlobalMaxLimits)                , mHistogramDim*sod);

    mRandom     -> backup(addPrefix("rng.bkup"));
    mHistogram  -> backup(fout);
    mModel      -> backup(fout);

    if(mParallelSimFlag){
        int num_long = 2 + 2* mModel -> size();
        int num_double = 1 +  mModel -> dim() * mModel -> size();
        
        fout.write(reinterpret_cast<char*>(&mReplica->mID)               ,sol);
        fout.write(reinterpret_cast<char*>(mReplica->mLongData)         ,sol*num_long);
        fout.write(reinterpret_cast<char*>(mReplica->mDoubleData)       , sod*num_double);
    }


    fout.write(reinterpret_cast<char*>(&mConfPrintNumber)               , sol );
    fout.write(reinterpret_cast<char*>(&mConfPrintFlag)                 , sizeof(bool));
    if(mConfPrintFlag){
        fout.write(reinterpret_cast<char*>(mConfPrintCount)                 , mHistogram->size() * sol);
    }
    fout.close();
}

template<class ModelType>
bool WLFrame< ModelType>::recover() {

    std::string  backupfile = backupFilename();
    std::ifstream fin(backupfile.c_str(), std::ios::binary);

    mRecoverFlag = false;
    if (fin.is_open()) {
        mRecoverFlag = true;


        *mOut << "\n##############\n"
              << "Wang-Landau Simulation recovering ... : \n" 
              << "  >> recover from file: " << backupfile << std::endl;
    
        long sod = sizeof (double);
        long sol = sizeof (long);

        fin.read(reinterpret_cast<char*>(&mParallelSimFlag)               , sizeof(bool));
        fin.read(reinterpret_cast<char*>(&mRandomSeed)                    , sol);
        fin.read(reinterpret_cast<char*>(&mID)                            , sol);
        fin.read(reinterpret_cast<char*>(&mReplicaExchangeProposeCount)   , sol);
        fin.read(reinterpret_cast<char*>(&mReplicaExchangeAcceptCount)    , sol);
        fin.read(reinterpret_cast<char*>(&mCurrentIteration)              , sol);
        fin.read(reinterpret_cast<char*>(&mHistogramCheckInterval)        , sol);
        fin.read(reinterpret_cast<char*>(&mReplicaExchangeInterval)       , sol);
        fin.read(reinterpret_cast<char*>(&mCompleteCheckInterval)         , sol);
        fin.read(reinterpret_cast<char*>(&mOutputInterval)                , sol);
        fin.read(reinterpret_cast<char*>(&mBackupInterval)                , sol);
        fin.read(reinterpret_cast<char*>(&mHistogramDim)                  , sol);
        fin.read(reinterpret_cast<char*>(&mSweep)                         , sol);
        fin.read(reinterpret_cast<char*>(&mCurrentMCMoves)                , sizeof(unsigned long));
        fin.read(reinterpret_cast<char*>(&mTotalMCMoves)                  , sizeof(unsigned long));
        fin.read(reinterpret_cast<char*>(&mFlatnessCriterion)             , sod);
        fin.read(reinterpret_cast<char*>(&mModificationFactor)            , sod);
        fin.read(reinterpret_cast<char*>(&mModificationDivider)           , sod);
        fin.read(reinterpret_cast<char*>(&mModificationThreshold)         , sod);


        mMinLimits = new double[mHistogramDim];
        mMaxLimits = new double[mHistogramDim];
        mGlobalMinLimits = new double[mHistogramDim];
        mGlobalMaxLimits = new double[mHistogramDim];
        mHistogramSize = new long[mHistogramDim];
        mHistogramType = new bool[mHistogramDim];

        fin.read(reinterpret_cast<char*>(mHistogramType)                  , mHistogramDim*sizeof(bool));
        fin.read(reinterpret_cast<char*>(mHistogramSize)                  , mHistogramDim*sol);
        fin.read(reinterpret_cast<char*>(mMinLimits)                      , mHistogramDim*sod);
        fin.read(reinterpret_cast<char*>(mMaxLimits)                      , mHistogramDim*sod);
        fin.read(reinterpret_cast<char*>(mGlobalMinLimits)                , mHistogramDim*sod);
        fin.read(reinterpret_cast<char*>(mGlobalMaxLimits)                , mHistogramDim*sod);

        
        mRandom     = new Random(mRandomSeed); 
        mHistogram  = new HistogramND(mHistogramDim);
        mModel      = new ModelType();
        
        mHistogram  ->  setOutput(*mOut);
        mModel      ->  setOutput(*mOut);

        mRandom     ->  recover(addPrefix("rng.bkup"));

        mModel      ->  setOutput(*mOut);
        mModel      ->  setRandom(*mRandom);
        mHistogram  ->  setOutput(*mOut);

        mHistogram  ->  recover(fin);
        mModel      ->  recover(fin);
       
        mModel -> setConstraint(mHistogram);
        
        if(mParallelSimFlag){
            int num_long = 2 + 2* mModel -> size();
            int num_double = 1 +  mModel -> dim() * mModel -> size();
            
            PWL_replicaRegisterAndInit();
            fin.read(reinterpret_cast<char*>(&mReplica->mID)               , sol);
            fin.read(reinterpret_cast<char*>(mReplica->mLongData)         , sol*num_long);
            fin.read(reinterpret_cast<char*>(mReplica->mDoubleData)       , sod*num_double);
        }

        fin.read(reinterpret_cast<char*>(&mConfPrintNumber)               , sol );
        fin.read(reinterpret_cast<char*>(&mConfPrintFlag)                 , sizeof(bool));

        if( mConfPrintFlag ){
            *mOut << "  >> configuration print flag is ON ! \n"
                  << "       * number of configurations that is going to print for each bin is " << mConfPrintNumber << std::endl;
            mConfPrintCount = new long[mHistogram->size()];
            fin.read(reinterpret_cast<char*>(mConfPrintCount)                 , mHistogram->size() * sol);
            getOutputStream(mConfPrintStream, mPrefix + DEFAULT_CONF_PRINT_FILENAME_BASE, true);
        } else{
            *mOut << "  >> configuration print flag is OFF !" << std::endl;
        }

        fin.close();

    }
    return mRecoverFlag;
}


/* 
 * **********************************************************
 * :::::::::::::::::::::::::::::::::::::::::::::::::::::::::
 *              
 *              Parallel WL Related Methods below
 * 
 * :::::::::::::::::::::::::::::::::::::::::::::::::::::::::
 */


template <class ModelType>
bool WLFrame<ModelType>::PWL_isHistogramFlat(MPI_Comm* mpi_intra_win_comm){
    int myflat, allflat;
    this -> isHistogramFlat() ? myflat = 1 : myflat = 0;
    int rc = MPI_Allreduce(&myflat, &allflat, 1, MPI_INT, MPI_PROD, *mpi_intra_win_comm);
    if (rc != MPI_SUCCESS) {
        errorMsg("PWL_isHistogramFlat", "MPI_Allreduce(...) was not successfully performed ! ");
        MPI_Abort(MPI_COMM_WORLD, rc);
    }
    return ( allflat == 0 ? false : true );
}

template <class ModelType>
bool WLFrame<ModelType>::PWL_isSimulationComplete(){
    int mystop, allstop;
    this -> isSimulationComplete() ? mystop = 1 : mystop = 0;
    int rc = MPI_Allreduce(&mystop, &allstop, 1, MPI_INT, MPI_PROD, MPI_COMM_WORLD);
    if (rc != MPI_SUCCESS) {
        errorMsg("PWL_isSimulationComplete", "MPI_Allreduce(...) was not successfully performed ! ");
        MPI_Abort(MPI_COMM_WORLD, rc);
    }
    
    return ( allstop == 0 ? false : true );
}

template <class ModelType>
bool WLFrame<ModelType>::PWL_proposeReplicaExchange(  int num_procs_per_win, 
                                                      int pid_inter_win_comm, 
                                                      int inter_win_comm_id, 
                                                      MPI_Comm* mpi_inter_win_comm ){
    /**
     * Steps to implement replica exchange
     *  1. find out pairs for everybody
     *  2. exchange information about observable to calculate "myratio"
     *  3. exchange "myratio" to decided whether replica exchange should be performed or not
     *  4. if yes, do exchange, otherwise do nothing
     */

    // check whether backup criterion satisfied
    if( mTotalMCMoves !=0 && (mTotalMCMoves % mBackupInterval ) == 0) { 
        if (! mRecoverFlag)
                backup();
        else 
            mRecoverFlag = false;
    }

    mCurrentMCMoves++;
    mTotalMCMoves++;

    mReplicaExchangeProposeCount ++;
    MPI_Status status;
    
    //:::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    // look for swap partner
    int* pairs = new int[2 * num_procs_per_win];
    int swap_partner = -1;
    
    if (pid_inter_win_comm == 0 && inter_win_comm_id != -1) {
        int choose_from = num_procs_per_win;
        int select;
        int* libre = new int[num_procs_per_win];

        for (int i = 0; i < num_procs_per_win; ++i) libre[i] = num_procs_per_win + i;

        // idea: processes from the lower box choose someone from the higher box at random
        for (int i = 0; i < num_procs_per_win; ++i) {
            select = mRandom -> nextLong(0, choose_from);
            pairs[i] = libre[select];
            pairs[libre[select]] = i;
            choose_from--;
            for (int j = select; j < choose_from; ++j) {
                libre[j] = libre[j + 1];
            }
        }
        delete [] libre;
    }
    
    // scatter swap partner information to everyone
    if(inter_win_comm_id != -1){
        MPI_Scatter(pairs, 1, MPI_INT, &swap_partner, 1, MPI_INT, 0, mpi_inter_win_comm[inter_win_comm_id]);
    }
    delete [] pairs;

    int change(0), isInMyLimit(0), isInOtherLimit(0);
    if(inter_win_comm_id != -1){
        //:::::::::::::::::::::::::::::::::::::::::::::::::::::::::
        // check probability to swap
        
        double* observables = NULL;
        int num_observables = mModel->numberOfObservables();
        double myfrac, otherfrac, accept_prob;
        observables = new double[num_observables];
        mModel -> copyObservables(observables); 
        
        MPI_Sendrecv_replace(observables, num_observables, MPI_DOUBLE, swap_partner, 1, swap_partner, 1, mpi_inter_win_comm[inter_win_comm_id], &status);

        // check whether the observable is inside the limits of both threads
        if(mHistogram -> isInLimit(observables)) {
            isInMyLimit = 1;
            myfrac = dosRatioInLog(mModel->observables(), observables);
        }
        else isInMyLimit = 0;

        if(pid_inter_win_comm < swap_partner){
            MPI_Recv(&isInOtherLimit, 1, MPI_INT, swap_partner, 4, mpi_inter_win_comm[inter_win_comm_id], &status);
        }else{
            MPI_Send(&isInMyLimit, 1, MPI_INT, swap_partner, 4, mpi_inter_win_comm[inter_win_comm_id]);
        }
        
        if(pid_inter_win_comm < swap_partner){
            MPI_Recv(&otherfrac, 1, MPI_DOUBLE, swap_partner, 2, mpi_inter_win_comm[inter_win_comm_id], &status);
            if(isInMyLimit == 1 && isInOtherLimit == 1){
                accept_prob = myfrac + otherfrac;
                accept_prob = accept_prob > 0 ? 1:exp(accept_prob);
                if(mRandom -> nextDouble() < accept_prob){
                    change =1;
                }else{
                    change =0;
                }
            }else change = 0;

            MPI_Send(&change, 1, MPI_INT, swap_partner, 3, mpi_inter_win_comm[inter_win_comm_id]);
        }else{
            MPI_Send(&myfrac, 1, MPI_DOUBLE, swap_partner, 2, mpi_inter_win_comm[inter_win_comm_id]);
            MPI_Recv(&change, 1, MPI_INT, swap_partner, 3, mpi_inter_win_comm[inter_win_comm_id], &status);
        }
        delete [] observables; 
        if(change == 1){
            mReplicaExchangeAcceptCount++;
            PWL_packupReplica(*mReplica);

            int rc=MPI_Sendrecv_replace(mReplica, 1, ReplicaType, swap_partner, 1, swap_partner, 1, mpi_inter_win_comm[inter_win_comm_id],&status);
            if(rc != MPI_SUCCESS){
                // something wrong, need recovery, to be implemented
                MPI_Abort(MPI_COMM_WORLD, rc);
            }
            PWL_unpackReplica(*mReplica);
        }
    }

    return (bool)change;

}

template <class ModelType>
void WLFrame<ModelType>::PWL_dosMerge(MPI_Comm* mpi_intra_win_comm, int num_procs_per_win){
    double *logge, *logge_buffer;
    int totalBins = this -> histogramSize();
    logge = new double[totalBins];
    logge_buffer = new double[totalBins];
    this -> dosCopyAll(logge);
    int status = MPI_Allreduce(logge, logge_buffer, totalBins, MPI_DOUBLE, MPI_SUM, *mpi_intra_win_comm);
    if (status != MPI_SUCCESS) {
        errorMsg("PWL_dosMerge", "MPI_Allreduce(...) was not successfully performed ! ");
        MPI_Abort(MPI_COMM_WORLD, status);
    }
    for (int i = 0; i < totalBins; ++i) {
        logge_buffer[i] /= num_procs_per_win;
    }
    this -> dosAssignAll(logge_buffer);
    delete [] logge;
    delete [] logge_buffer;
}

template <class ModelType>
void WLFrame<ModelType>::PWL_packupReplica(Replica& replica){
    mModel -> packup(replica.mLongData, replica.mDoubleData);
}

template <class ModelType>
void WLFrame<ModelType>::PWL_unpackReplica(Replica& replica){
    mModel -> unpack(replica.mLongData, replica.mDoubleData);
}

template <class ModelType>
void WLFrame<ModelType>::PWL_replicaRegisterAndInit(){
     //:::::::::::::::::::::::::::::::::::::::;::::::::::::::
    // initialize package and register package type for MPI 
    int num_long = 2 + 2* mModel -> size();
    int num_double = 1 +  mModel -> dim() * mModel -> size();
    *mOut << "      ->  register PACKAGE for MPI... \n"
          << "          * No. of MPI_LONG: " << num_long        << "\n"
          << "          * No. of MPI_Double: " << num_double    << std::endl;
           
    mReplica = new Replica(num_long, num_double);
    mReplica -> mID = mID;

    int num_block = 3;
    int blocklen[3] = {1, num_long, num_double};
    MPI_Datatype old_types[3] = {MPI_LONG, MPI_LONG, MPI_DOUBLE}; 
    //MPI_Aint displacement[3]  = {0, sizeof(long), sizeof(long)*(num_long+1)};
    MPI_Aint displacement[3];
    MPI_Aint start_address, address;
    MPI_Address(& (mReplica -> mID), &start_address);
    displacement[0]=0;
    MPI_Address( mReplica -> mLongData, &address);
    displacement[1]= address - start_address;
    MPI_Address( mReplica -> mDoubleData, &address);
    displacement[2]= address - start_address;
    
    MPI_Type_struct(num_block, blocklen, displacement, old_types ,&ReplicaType);
    MPI_Type_commit(&ReplicaType);

}

template <class ModelType>
void WLFrame<ModelType>::PWL_init(int dim, int* num_wins, int* id, double* overlap){

    mParallelSimFlag = true;
    //:::::::::::::::::::::::::::::::::::::::::::
    // modify range of observables based on id
    double range;
    double width;
    double* step = new double[dim];
    for(int d=0; d<dim; ++d){
        if(mHistogramType[d]){
            step[d] = floor( (mMaxLimits[d]-mMinLimits[d]+1)/(double)mHistogramSize[d]+0.5 );
            range = step[d]*mHistogramSize[d];
            width = floor( range/(1 + (1-overlap[d])*(num_wins[d]-1)) + 0.5 );
            mMinLimits[d] = mMinLimits[d]+floor((1-overlap[d])*width+0.5)*id[d];
            mMaxLimits[d] = mMinLimits[d]+width-1;
            mHistogramSize[d] = floor( mHistogramSize[d] * 1.0/(1 + (1-overlap[d])*(num_wins[d]-1)) + 0.5 ); 
        }else{
            step[d] = (mMaxLimits[d]-mMinLimits[d])/(double)mHistogramSize[d];
            range = step[d]*mHistogramSize[d];
            width = range/(1 + (1-overlap[d])*(num_wins[d]-1));
            mMinLimits[d] = mMinLimits[d]+(1-overlap[d])*id[d]*width;
            mMaxLimits[d] = mMinLimits[d]+width;
            mHistogramSize[d] = floor(mHistogramSize[d]/(1 + (1-overlap[d])*(num_wins[d]-1))+0.5); 
        }
    }
    delete [] step;

    // regular initialization
    init();
    
    // if not recover, register package
    if(!mRecoverFlag)    PWL_replicaRegisterAndInit();

}

#endif	/* WANGLANDAU_H */

