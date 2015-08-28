/* *********************************************************
 * :::::::::::::::::::::::::::::::::::::::::::::::::::::::::
 *
 *   This file contains the definition of classes LipidModel
 *
 *   @File:   LipidModel.h/.cpp
 *   @Author: Guangjie Shi (Jerry)
 *
 *   @Version 1.0: Nov. 12, 2013
 *   @Version 1.1: Apr.  8, 2013
 *      In this version, some helper arrays have been added.
 *      They are mCellID2Coord, mCoord2CellID, mCellCornerCoord respectively.
 *      The purpose of this modification is to speed up the calculation.
 *      However, at the same time, some issue raised, in this version only
 *      3 dimension is considered. Therefore, one has to made certain 
 *      modification if s/he wants to run 2 dimension simulation
 *                
 *
 * :::::::::::::::::::::::::::::::::::::::::::::::::::::::::
 */


#ifndef LIPIDMODEL_H
#define	LIPIDMODEL_H

#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include "Monomer.h"
#include "Random.h"
#include "Global.h"
#include "HistogramND.h"
#include "RandomAccessNeighborList.h"


enum BoundaryConditionType{
    PeriodicBoundaryCondition=0
};

enum TrialMoveType{
    RANDOMMOVE = 0,
    INSERTIONMOVE,
    DELETIONMOVE,
    NONEMOVE, 
    RANDOMSHIFTMOVE,
    REPTATIONMOVE
};

const double DOUBLE_TOLERANCE = 1e-7;   // tolerance for double equaling
const long MAXBOND = 2;                 // maximum bonds one monomer could obtain
const long NUMBEROFSTATS = 10;           // number of statistics
const long NUMBEROFMOVETYPES = 3;       // number of trial move types, does not include NONEMOVE

const long   ADAWL_RECALCULATE_MOVE_DISTANCE_THRESHOLD = 10000;
const double OPTIMAL_ACCEPT_PROB = 0.4;
//const double ADAWL_CONSTANT_A = 0.82988;
//const double ADAWL_CONSTANT_B = 0.014625;
const double ADAWL_CONSTANT_A = 0.2;
const double ADAWL_CONSTANT_B = 0.5;



class LipidModel{

public:
    
    /**
     * Constructor
     */
    LipidModel();

    /**
     * Constructor
     * @param in input file stream
     */
    LipidModel(std::istream& in);

    /**
     * Copy Constructor
     */
    LipidModel(const LipidModel& orig);

    /**
     * Destructor
     */
    ~LipidModel();      
    
    /**
     * Get input from file
     * @param filename input file name
     */
    void input(std::string filename);

    /**
     * Get input from input file stream
     * @param in input file stream
     */
    void input(std::istream& in);

    /**
     * initialize whole system
     */
    void init(std::string initConfFile=NULL);
    
    /**
     * calculate the total energy of the system
     * @return the potential energy of whole system
     */
    double potential();
    

    /**
     * calculate potential between two monomers based on monomer types
     * @param m Monomer one
     * @param n Monomer two
     * @return potential energy between two monomers
     */
    double potential(const Monomer& m, const Monomer& n);
    
    /**
     * @param index
     * @param flag
     * @return the potential energy related to certain atom
     */
    double potential(const long& index, const bool& countWaterflag = false);
    

    /**
     * measure the system and assign according values to observalbes 
     */
    void doMeasure();
    
    /**
     * do Monte Carlo trial move
     * @return true if move successfully, false if move is rejected
     */
    void doMCMove();
     
    /**
     * undo previous Monte Carlo trial move
     */
    void undoMCMove();
    
    /**
     * overwrite operator =
     * @param one object of OLModel
     * @return 
     */
    LipidModel& operator=(const LipidModel& orig);
    
    
    /**
     * This method provide the function of recoding important information of this class
     * @param out
     */
    void print(std::ostream& out=std::cout, PrintType type = PRINT_MODEL_INFORMATION);


    /**
     * This method will read in the model configuration stored in the file
     * fileformat should be mol2 format that print out by this program
     * @param fin
     */
    void readModel(std::ifstream& fin);
    
    /**
     * @return correction factor for Monte Carlo trial move
     */
    double cFactor();

    /**
     * this function is primarily used for replica exchange
     */
    double cFactor(const double* obs_up, const double* obs_down);

    /**
     * @return total number of monomers in the system
     */
    inline long size(){
        return mNumberOfMonomers;
    }

    inline long numberOfObservables(){
        return mNumberOfObservables;
    }

    /**
     * @return dimension of model
     */
    inline long dim(){
        return mDim;
    }

    /**
     * @return energy of current system
     */
    inline double energy(){
        return mEnergy;
    }

    /**
     * @return energy difference caused by current MC move
     */
    inline double energyDifference(){
        return (mEnergy - mPrevEnergy);
    }
    
    /**
     * check model correctness
     */
    void checkModelIntegrity();

    /**
     * set model output method
     */
    void setOutput(std::ostream& out=std::cout);

    /**
     * Given a instance of HistogramND class, specify the constraints
     * of this model. 
     * @param his a pointer of HistogramND class
     */
    void setConstraint(const HistogramND* histogram = NULL);

    /**
     * set random number generator
     */
    void setRandom(Random& mRandom);

    /**
     * Packup all informations for replica exchange
     * @param dataLong long array for storing long data
     * @param dataDouble double array for storing double data
     */
    void packup(long* longData, double* doubleData);
    
    /**
     * Unpack from replica, and rebuild the whole system base on it. 
     * @param dataLong long array for storing long data
     * @param dataDouble double array for storing double data
     */
    void unpack(long* longData, double* doubleData);
     
    void backup(std::ofstream& fout);
    void recover(std::ifstream& fin);

    // clear all the information MC move involved.
    void resetStatistics();

    inline bool satisfyLimits(){
        return ( ( mMinEnergy <= mEnergy ) && 
                 ( mMaxEnergy >= mEnergy ) && 
                 (mNumberOfLipids >= mMinNumberOfLipids || mMinNumberOfLipids == -1) && 
                 (mNumberOfLipids <= mMaxNumberOfLipids || mMaxNumberOfLipids == -1) );
    }

    inline const double* observables(){
        return mObservables;
    }

    inline const double* prevObservables(){
        return mBackupObservables;
    }

    inline void copyObservables(double* obs){
        for(int i=0; i<mNumberOfObservables; ++i)
            obs[i] = mObservables[i];
    }


private:

    const HistogramND* mHistogram;

    Random* mRandom; 
    std::ostream* mOut;             // output stream
    Cell* mCells;                   // cell structure, used as storage of monomer index
    Monomer* mMonomers;             // monomer array
    //Monomer* mMovedMonomers;        // this array is help to record the monomer been moved in last MC move
    //                                // however, for single monomer RandomMove, this array is not used, instead mBackupMonomer is used
    //                                // this array does not require backup operation during the simulation
    MonomerIndex* mMonomerIndexes;  // helper that speed up the access of the mMonomers in boxes
    Monomer mBackupMonomer;         // backup monomer for random displacement move
    Monomer mBackupLipid[3];           // backup monomers of a lipid for random shift move
    TrialMoveType mMoveProposal;    // trail move proposal, helper for calculating the correction factor
    BoundaryConditionType mBoundaryCondition;


    bool mInitType;                 //Flag for init configuration, false for randomly, true for input from file
    bool mAdaptiveMoveFlag;         // Flag for adaptive maximum move distance
    bool mReptationFlag;            // True crawl toward head, otherwise false  
    long mNumberOfCells;            // number of cells in whole space
    long mNumberOfCellsPerLine;     // number of cells per line
    long mNumberOfCellsPerLine2;    // cellPerLine*cellPerLine
    long mNumberOfMonomers;         // number of Monomers in whole space
    long mOrder;                    // order of number of monomers, to help fix output format
    long mNumberOfObservables;       
    long mDim;                      // model dimension
    long mLastDeletedLipid[3];
    long mMaxNumberOfLipids;
    long mMinNumberOfLipids;
    long mNumberOfLipids;
    long mLastCreatedLipidHeadIDInLipidWater;
    long mNumberOfMovedMonomers;    // record how many monomers been moved in last MC move, not including single monomer RandomMove
                                    // this variable does not require backup operation during the simulation
    long mMovedBinIndex;
    long* mMoveCount;
    long* mMoveRejectCount;
    long* mLipidWater;              // helper that speed up the process of picking up certain type of mMonomers randomly 
                                    // [0, 3*numberOfLipids) are lipids, and rest are water mMonomers
    long** mBonds;                  // array that record information of bonds
    long** mCellID2Coord;           // map cell index to cartesian coordinates
    long*** mCoord2CellID;          // map coordiate of a cell to its cell id

    MathVector<double>** mCellCornerCoord;      // coordinates of lowest corner of cells  
                                                //          -------- 
                                                //        /|        /| 
                                                //       / |       / |
                                                //       ---------   |
                                                //      |   ------|- /
                                                //      | /       | /
                                                //      |/        |/  
                                                //   (*) --------- 

    double mLargestCutoff;

    double mReptationCFactor;
    double mMoveFraction[3];        // [0] random local displacement [1] insertion or deletion move
    double mMaxEnergy;
    double mMinEnergy;
    double mBondLength[3];          // involved in creation of new bond
    double mEnergy;                 // total energy
    double mPrevEnergy;             // energy before trial move
    double mLengthOfSimBox;         // length of simulation box
    double mDensity;                // density of system
    double mMinLengthOfCell;        // minimum cell size ( equal to the cutoff distance)
    double mLengthOfCell;           // cell size
    double mInitMaxMoveDistance;        // max distance the random move could make
    double* mMaxMoveDistance;        // max distance the random move could make
    double* mChangeOfCoord;         // coordinate changed by random move
    double* mObservables;           // generic variable for HistogramND
    double* mBackupObservables; 
    

    void errorMsg(std::string func, std::string msg);
    void buildCells();
    void buildCellHelpers();
    void buildMonomerIndexes();
    void initAdapTrialMove();
    void updateTrialMoveDistance(long i);
    void addBond(const long& a, const long& b);
    void delBond(const long& a, const long& b);
    bool canBeBonded(const Monomer& m, const Monomer& n);

    /**
     * Delete current bonds information, and rebuild bonds relied on 
     * current information of lipids;
     */
    void initBonds();

    /**
     * backup obseravables before trial moves
     */
    void backupMeasure();
    
    /**
     * test whether two given mMonomers are connected with a bond
     */
    bool isBonded(const long& a,const long& b);

    /**
     * randomly insert a lipid
     * @return true if sucessfully, false otherwise
     */
    bool insertLipid();
    
    /**
     * given three monomers and construct a lipid 
     * @head polar head
     * @middle hydrophobic middle 
     * @tail hydrophobic tail
     * @return true if sucessfully, false otherwise
     */
    void insertLipid(const long& head, const long& middle, const long& tail);
    void undoInsertLipid();
    
    /**
     * delete one lipid of the system, lipid will be randomly picked up if not index is specified
     * @param i index of lipid
     */
    void deleteLipid(long i=-1);
    void undoDeleteLipid();
    
    /**
     * Perform random local displacement move for given index, or random index if not specified
     * @param index index of monomer to be randomly moved, random value will be assigned if not specified
     */
    double doRandomMove(long index = -1);
    void undoRandomMove();
    
    /**
     * Perform random local displacement move for a given lipid, or random lipid if not specified
     * @param index index of lipid
     */
    double doRandomShiftMove(long index = -1);
    void undoRandomShiftMove();


    /**
     * Perform reptation move for lipid
     */
    bool doReptationMove(long index = -1);
    void undoReptationMove();

    /**
    * Given index and coordinates to move according monomer
    * The function does not deal with backup of moved monomer.
    * Thus one need to backup the monomers before calling this function, 
    * in order to efficiently undo the move.
    * @param index index of monomer to be moved
    * @param coord change of coordinates
    */
    void moveMonomer(const long& index, double* coord);

    /**
     * Given indexes and coordiantes to move a group of monomers,
     * This function simply just call moveMonomer(...) several times.
     * @param index array of indexes of monomer to move
     * @param size of the array
     * @param coord change of coordinates
     */
    void moveMonomerGroup(const long* index, const long& size, double* coord);

    
    /**
     * imigrate atom index, from old box to new box
     * @param index
     * @param coord
     */
    void migrate(const long& index, const long& oldB, const long& newB);

    /**
     * Given the index of an atom, find its neighbor list
     * @param id index of given atom
     * @return neighbor list of given atom
     */
    RandomAccessNeighborList findWaterNeighborList(const long& id);

    /**
     * To test whether a cell is close enough to a given monomer,
     * such that monomers that are contained in this cell could possibly interact
     * with the monomer.
     *
     * @param monomer_id index of given monomer
     * @param t_cell_id index of target cell
     * @param r distance to judge whether is close enough
     * @return true if they are close enough, otherwise false
     */
    bool isCellClose(const long& monomer_id, const long& t_cell_id, const double& r); 


    
    // :::::::::::::::::::::::::::::::::::::::::::: //
    //                                              //
    //          inline helper functions             //
    //                                              //
    // :::::::::::::::::::::::::::::::::::::::::::: //
    
    inline void applyBoundaryCondition(Monomer& a){
        a.mCoord -> periodic(0,mLengthOfSimBox);
    }
     
    inline long cellIndex(const long& ix, const long& iy, const long& iz=0){
        return  ix + iy * mNumberOfCellsPerLine + iz*mNumberOfCellsPerLine2;
    } 

    inline long whichCell(const double& x, const double& y, const double& z = 0){
       return mCoord2CellID[(long)floor(x / mLengthOfCell)][(long)floor(y / mLengthOfCell)][(long)floor(z / mLengthOfCell)];
       // return (long)floor(x / mLengthOfCell) + 
       //        (long)floor(y / mLengthOfCell)*mNumberOfCellsPerLine +
       //        (long)floor(z / mLengthOfCell)*mNumberOfCellsPerLine2;
    }
    
    inline long whichCell(const Monomer& a){
        return mCoord2CellID[(long)floor(a.mCoord->get(0) / mLengthOfCell)][(long)floor(a.mCoord->get(1) / mLengthOfCell)][(long)floor(a.mCoord->get(2) / mLengthOfCell)];
        //return (long)floor(a.mCoord->get(0) / mLengthOfCell) +  
        //       (long)floor(a.mCoord->get(1) / mLengthOfCell)*mNumberOfCellsPerLine +
        //       (long)floor(a.mCoord->get(2) / mLengthOfCell)*mNumberOfCellsPerLine2;
    }

    void swapLipidWaterArray(long i, long j){
        long mi = mLipidWater[i], mj = mLipidWater[j];
        mMonomers[mi].mLipidWaterID = j;
        mMonomers[mj].mLipidWaterID = i;
        mLipidWater[i] = mj;
        mLipidWater[j] = mi;
    }

    // :::::::::::::::::::::::::::::::::::::::::::::::: //
    //                                                  //
    //     Statistics:                                  //
    //         0: number of random move                 //
    //         1: number of rejected random move        //
    //         2: number of insertion move              //
    //         3: number of rejected insertion move     //
    //         4: number of deletion move               //
    //         5: number of rejected deletion move      //
    //                                                  //
    // :::::::::::::::::::::::::::::::::::::::::::::::: //
    
    long stat[NUMBEROFSTATS];

    // :::::::::::::::::::::::::::::::::::::::::::::::::::: //
    //                                                      //
    //          Parameters and routines for                 //
    //          calculation of potential energy             //
    //                                                      //
    // :::::::::::::::::::::::::::::::::::::::::::::::::::: //
    
    double potentialLJ(const Monomer& m1, const Monomer& m2);
    double potentialSCR(const Monomer& m1, const Monomer& m2);
    double potentialFENE(const Monomer& m1, const Monomer& m2);

    double LJ_cutoff[NUMBEROFMONOMERTYPES][NUMBEROFMONOMERTYPES];     
    double LJ_sigma[NUMBEROFMONOMERTYPES][NUMBEROFMONOMERTYPES];      
    double LJ_epsilon[NUMBEROFMONOMERTYPES][NUMBEROFMONOMERTYPES];    
    double LJ_sigma_6[NUMBEROFMONOMERTYPES][NUMBEROFMONOMERTYPES];    
    double LJ_shift[NUMBEROFMONOMERTYPES][NUMBEROFMONOMERTYPES];      
    double SCR_epsilon;     
    double SCR_sigma;       
    double SCR_sigma_9;     
    double SCR_cutoff;      
    double SCR_shift;       
    double FENE_K;          
    double FENE_R;          
    double FENE_R2;         
    double FENE_ro;         

};

#endif	/* LIPIDMODEL_H */

