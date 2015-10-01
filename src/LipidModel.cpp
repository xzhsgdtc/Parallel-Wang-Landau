#include "LipidModel.h"



LipidModel::LipidModel(){
    setOutput(std::cout);
    mCells              = NULL;
    mMonomers           = NULL;
    mMonomerIndexes     = NULL;
    mLipidWater         = NULL;
    mBonds              = NULL;
    mCellID2Coord       = NULL;
    mCoord2CellID       = NULL;
    mCellCornerCoord    = NULL;
    mChangeOfCoord      = NULL;
    mObservables        = NULL;
    mBackupObservables  = NULL;
    mRandom             = NULL;
    mMaxMoveDistance    = NULL;
    mHistogram          = NULL;
    mMoveCount          = NULL;
    mMoveRejectCount    = NULL;
    mMovedBinIndex      = -1;
    mBoundaryCondition = PeriodicBoundaryCondition;
    mAdaptiveMoveFlag = false;
    mInitType = false;
    mIsMoveValid = true;
    mLargestCutoff     = 0;
    mNeighborWaterAtInsertion = 0;
    mNeighborWaterAtDeletion = 0;
    
    mNumberOfLipids = mMinEnergy = mMaxEnergy = mMinNumberOfLipids = mMaxNumberOfLipids = 0;
    mNumberOfObservables = 1;
    mVolumeChangeRate = 0.01;
    mMoveFraction[0] = 1;
    mMoveFraction[1] = 1;
    mMoveFraction[2] = 1;
    for(int i=0; i<NUMBEROFSTATS; ++i){
        stat[i] = 0;
    }

    for(int i=0; i<NUMBEROFMONOMERTYPES; ++i){
        for(int j=0; j<NUMBEROFMONOMERTYPES; ++j){
            LJ_cutoff[i][j] = 0;
            FENE_WCA_cutoff[i][j] = 0;
            LJ_sigma[i][j] = 0;
            LJ_epsilon[i][j] = 0;
            LJ_sigma_6[i][j] = 0;
            LJ_shift[i][j] = 0;
        }
    }
    mOut = &std::cout; 
}

LipidModel::LipidModel(const LipidModel& orig){

    mCells              = NULL;
    mMonomers           = NULL;
    mMonomerIndexes     = NULL;
    mLipidWater         = NULL;
    mBonds              = NULL;
    mCellID2Coord       = NULL;
    mCellCornerCoord    = NULL;
    mCoord2CellID       = NULL;
    mChangeOfCoord      = NULL;
    mObservables        = NULL;
    mBackupObservables  = NULL;
    mRandom             = NULL;
    mMaxMoveDistance    = NULL;
    mHistogram          = orig.mHistogram;

    if(mHistogram != NULL){
        mMaxMoveDistance = new double[mHistogram->size()];
        mMoveCount = new long[mHistogram->size()];
        mMoveRejectCount = new long[mHistogram->size()];
        for(int i=0; i<mHistogram->size(); ++i){
            mMaxMoveDistance[i] = orig.mMaxMoveDistance[i];
            mMoveCount[i] = orig.mMoveCount[i];
            mMoveRejectCount[i] = orig.mMoveRejectCount[i];
        }
    }

    setOutput(*orig.mOut);
    setRandom(*orig.mRandom);
    
    mIsMoveValid = orig.mIsMoveValid;
    mLargestCutoff = orig.mLargestCutoff;
    mInitType = orig.mInitType;              
    mAdaptiveMoveFlag = orig.mAdaptiveMoveFlag;
    mNumberOfCells = orig.mNumberOfCells;         
    mNumberOfCellsPerLine = orig.mNumberOfCellsPerLine;  
    mNumberOfCellsPerLine2 = orig.mNumberOfCellsPerLine2; 
    mNumberOfMonomers   = orig.mNumberOfMonomers;      
    mOrder  = orig.mOrder;                 
    mNumberOfObservables = orig.mNumberOfObservables;   
    mDim    = orig.mDim;                   
    mMaxNumberOfLipids = orig.mMaxNumberOfLipids;
    mMinNumberOfLipids  = orig.mMinNumberOfLipids;
    mNumberOfLipids =   orig.mNumberOfLipids;
    mLastCreatedLipidHeadIDInLipidWater = orig.mLastCreatedLipidHeadIDInLipidWater;
    mNeighborWaterAtInsertion = orig.mNeighborWaterAtInsertion;
    mNeighborWaterAtDeletion = orig.mNeighborWaterAtDeletion;

    for(int i=0; i<3; ++i){
        mLastDeletedLipid[i] = orig.mLastDeletedLipid[i];
        mBondLength[i] = orig.mBondLength[i];
    }

    mBackupMonomer = orig.mBackupMonomer;
    mMoveProposal = orig.mMoveProposal;
    mBoundaryCondition = orig.mBoundaryCondition;

    for(int i=0; i<3; ++i) mMoveFraction[i] = orig.mMoveFraction[i];
    mMaxEnergy  = orig.mMaxEnergy;
    mMinEnergy  = orig.mMinEnergy;
    mEnergy = orig.mEnergy;         
    mPrevEnergy = orig.mPrevEnergy;     
    mLengthOfSimBox = orig.mLengthOfSimBox; 
    mHigherBoundOfSimBox = orig.mHigherBoundOfSimBox; 
    mLowerBoundOfSimBox = orig.mLowerBoundOfSimBox; 
    mDensity    = orig.mDensity;        
    mMinLengthOfCell    = orig.mMinLengthOfCell;
    mLengthOfCell   = orig.mLengthOfCell;   
    mInitMaxMoveDistance    = orig.mInitMaxMoveDistance;

    if(mNumberOfCells > 0){
        mCells = new Cell[mNumberOfCells];
        for(int i=0; i<mNumberOfCells; ++i) mCells[i] = orig.mCells[i];
        mCellID2Coord = new long*[mNumberOfCells];
        mCellCornerCoord = new MathVector<double>*[mNumberOfCells];
        for(int i=0; i<mNumberOfCells; ++i){
            mCellID2Coord[i] = new long[mDim];
            for(int j=0; j<mDim; ++j) {
                mCellID2Coord[i][j] = orig.mCellID2Coord[i][j];
            }
            mCellCornerCoord[i] = new MathVector<double>(mDim);
            *mCellCornerCoord[i] = *orig.mCellCornerCoord[i];
        }

    

        if(mNumberOfCellsPerLine == 0 || mNumberOfCellsPerLine == 1){
            mCoord2CellID = new long**[1];
            mCoord2CellID[0] = new long*[1];
            mCoord2CellID[0][0] = new long[1];
        }else if(mNumberOfCellsPerLine > 1){
            mCoord2CellID = new long**[mNumberOfCellsPerLine];
            for(int i=0; i<mNumberOfCellsPerLine; ++i){
                mCoord2CellID[i] = new long*[mNumberOfCellsPerLine];
                for(int j=0; j<mNumberOfCellsPerLine; ++j){
                    mCoord2CellID[i][j] = new long[mNumberOfCellsPerLine];
                    for(int k=0; k<mNumberOfCellsPerLine; ++k)
                        mCoord2CellID[i][j][k] = orig.mCoord2CellID[i][j][k];
                }
            }
        }
    }
    if(mNumberOfMonomers > 0){
        mMonomers = new Monomer[mNumberOfMonomers];
        mLipidWater = new long[mNumberOfMonomers];
        initBonds();
        for(int i=0; i<mNumberOfMonomers; ++i) {
            mMonomers[i] = orig.mMonomers[i];
            mLipidWater[i] = orig.mLipidWater[i];
            for(int j=0; j<MAXBOND;++j)
                mBonds[i][j] = orig.mBonds[i][j];
        }
    }
    if(mNumberOfCells > 0 && mNumberOfMonomers > 0) buildMonomerIndexes();
    
    if(mDim > 0){
        mChangeOfCoord = new double[mDim];
        for(int i=0; i < mDim; ++i) mChangeOfCoord[i] = orig.mChangeOfCoord[i];
    }

    if(mNumberOfObservables > 0){
        mObservables = new double[mNumberOfObservables];
        mBackupObservables = new double[mNumberOfObservables];
        for(int i=0; i<mNumberOfObservables; ++ i){
            mObservables[i] = orig.mObservables[i];
            mBackupObservables[i] = orig.mBackupObservables[i];
        }
    }

    for(int i=0; i<NUMBEROFSTATS;++i) stat[i] = orig.stat[i];

    // copy parameters for potential
    for(int i=0; i<NUMBEROFMONOMERTYPES; ++i){
        for(int j=0; j<NUMBEROFMONOMERTYPES; ++j){
            LJ_cutoff[i][j]   =orig.LJ_cutoff[i][j];  
            LJ_sigma[i][j]    =orig.LJ_sigma[i][j];  
            LJ_epsilon[i][j]  =orig.LJ_epsilon[i][j];
            LJ_sigma[i][j]    =orig.LJ_sigma[i][j];  
            LJ_shift[i][j]    =orig.LJ_shift[i][j];  
            FENE_WCA_cutoff[i][j] = orig.FENE_WCA_cutoff[i][j];
        }
    }
    SCR_epsilon     = orig.SCR_epsilon ;     
    SCR_sigma       = orig.SCR_sigma   ;       
    SCR_sigma_9     = orig.SCR_sigma_9 ;     
    SCR_cutoff      = orig.SCR_cutoff  ;      
    SCR_shift       = orig.SCR_shift   ;       
    FENE_K          = orig.FENE_K      ;          
    FENE_R          = orig.FENE_R      ;          
    FENE_R2         = orig.FENE_R2     ;         
    FENE_ro         = orig.FENE_ro     ;         
}


LipidModel& LipidModel::operator=(const LipidModel& orig){
    if(this == &orig) return *this;

    if ( mCells              != NULL) delete [] mCells              ;
    if ( mMonomers           != NULL) delete [] mMonomers           ;
    if ( mMonomerIndexes     != NULL) delete [] mMonomerIndexes     ;
    if ( mLipidWater         != NULL) delete [] mLipidWater         ;
    if ( mChangeOfCoord      != NULL) delete [] mChangeOfCoord      ;
    if ( mObservables        != NULL) delete [] mObservables        ;
    if ( mBackupObservables  != NULL) delete [] mBackupObservables  ;
    if ( mMoveCount          != NULL) delete [] mMoveCount          ;
    if ( mMoveRejectCount    != NULL) delete [] mMoveRejectCount    ;
    if ( mMaxMoveDistance    != NULL) delete [] mMaxMoveDistance    ;

    if( mBonds              != NULL ){
        for(int i=0; i<mNumberOfMonomers; ++i)
            delete mBonds[i];
        delete [] mBonds;
    }
    if( mCellID2Coord       != NULL ){
        for(int i=0; i<mNumberOfCells; ++i)
            delete mCellID2Coord[i];
        delete [] mCellID2Coord;
    }
    if(mCoord2CellID        != NULL ){
        for(int i=0; i<mNumberOfCellsPerLine; ++i){
            for(int j=0; j<mNumberOfCellsPerLine; ++j){
                delete [] mCoord2CellID[i][j];
            }
            delete [] mCoord2CellID[i];
        }
        delete [] mCoord2CellID;
    }
    if(mCellCornerCoord != NULL){
        for(int i=0; i<mNumberOfCells; ++i){
            delete mCellCornerCoord[i];
        }
        delete [] mCellCornerCoord;
    }

    mHistogram          = orig.mHistogram;

    if(mHistogram != NULL){
        mMaxMoveDistance = new double[mHistogram->size()];
        mMoveCount = new long[mHistogram->size()];
        mMoveRejectCount = new long[mHistogram->size()];
        for(int i=0; i<mHistogram->size(); ++i){
            mMaxMoveDistance[i] = orig.mMaxMoveDistance[i];
            mMoveCount[i] = orig.mMoveCount[i];
            mMoveRejectCount[i] = orig.mMoveRejectCount[i];
        }
    }

    setOutput(*orig.mOut);
    setRandom(*orig.mRandom);
    
    mIsMoveValid = orig.mIsMoveValid;
    mLargestCutoff = orig.mLargestCutoff;
    mInitType = orig.mInitType;              
    mAdaptiveMoveFlag = orig.mAdaptiveMoveFlag;
    mNumberOfCells = orig.mNumberOfCells;         
    mNumberOfCellsPerLine = orig.mNumberOfCellsPerLine;  
    mNumberOfCellsPerLine2 = orig.mNumberOfCellsPerLine2; 
    mNumberOfMonomers   = orig.mNumberOfMonomers;      
    mOrder  = orig.mOrder;                 
    mNumberOfObservables = orig.mNumberOfObservables;   
    mDim    = orig.mDim;                   
    mMaxNumberOfLipids = orig.mMaxNumberOfLipids;
    mMinNumberOfLipids  = orig.mMinNumberOfLipids;
    mNumberOfLipids =   orig.mNumberOfLipids;
    mLastCreatedLipidHeadIDInLipidWater = orig.mLastCreatedLipidHeadIDInLipidWater;
    mNeighborWaterAtInsertion = orig.mNeighborWaterAtInsertion;
    mNeighborWaterAtDeletion = orig.mNeighborWaterAtDeletion;

    for(int i=0; i<3; ++i){
        mLastDeletedLipid[i] = orig.mLastDeletedLipid[i];
        mBondLength[i] = orig.mBondLength[i];
    }

    mBackupMonomer = orig.mBackupMonomer;
    mMoveProposal = orig.mMoveProposal;
    mBoundaryCondition = orig.mBoundaryCondition;

    for(int i=0; i<3; ++i) mMoveFraction[i] = orig.mMoveFraction[i];
    mMaxEnergy  = orig.mMaxEnergy;
    mMinEnergy  = orig.mMinEnergy;
    mEnergy = orig.mEnergy;         
    mPrevEnergy = orig.mPrevEnergy;     
    mLengthOfSimBox = orig.mLengthOfSimBox; 
    mHigherBoundOfSimBox = orig.mHigherBoundOfSimBox; 
    mLowerBoundOfSimBox = orig.mLowerBoundOfSimBox; 
    mDensity    = orig.mDensity;        
    mMinLengthOfCell    = orig.mMinLengthOfCell;
    mLengthOfCell   = orig.mLengthOfCell;   
    mInitMaxMoveDistance    = orig.mInitMaxMoveDistance;

    if(mNumberOfCells > 0){
        mCells = new Cell[mNumberOfCells];
        for(int i=0; i<mNumberOfCells; ++i) mCells[i] = orig.mCells[i];
        mCellID2Coord = new long*[mNumberOfCells];
        mCellCornerCoord = new MathVector<double>*[mNumberOfCells];
        for(int i=0; i<mNumberOfCells; ++i){
            mCellID2Coord[i] = new long[mDim];
            for(int j=0; j<mDim; ++j) mCellID2Coord[i][j] = orig.mCellID2Coord[i][j];
            
            mCellCornerCoord[i] = new MathVector<double>(mDim);
            *mCellCornerCoord[i] = *orig.mCellCornerCoord[i];
        }

        if(mNumberOfCellsPerLine == 0 || mNumberOfCellsPerLine == 1){
            mCoord2CellID = new long**[1];
            mCoord2CellID[0] = new long*[1];
            mCoord2CellID[0][0] = new long[1];
        }else if(mNumberOfCellsPerLine > 1){
            mCoord2CellID = new long**[mNumberOfCellsPerLine];
            for(int i=0; i<mNumberOfCellsPerLine; ++i){
                mCoord2CellID[i] = new long*[mNumberOfCellsPerLine];
                for(int j=0; j<mNumberOfCellsPerLine; ++j){
                    mCoord2CellID[i][j] = new long[mNumberOfCellsPerLine];
                    for(int k=0; k<mNumberOfCellsPerLine; ++k)
                        mCoord2CellID[i][j][k] = orig.mCoord2CellID[i][j][k];
                }
            }
        }

    }
    if(mNumberOfMonomers > 0){
        mMonomers = new Monomer[mNumberOfMonomers];
        mLipidWater = new long[mNumberOfMonomers];
        initBonds();
        for(int i=0; i<mNumberOfMonomers; ++i) {
            mMonomers[i] = orig.mMonomers[i];
            mLipidWater[i] = orig.mLipidWater[i];
            for(int j=0; j<MAXBOND;++j)
                mBonds[i][j] = orig.mBonds[i][j];
        }
    }
    if(mNumberOfCells > 0 && mNumberOfMonomers > 0) buildMonomerIndexes();
    
    if(mDim > 0){
        mChangeOfCoord = new double[mDim];
        for(int i=0; i < mDim; ++i) mChangeOfCoord[i] = orig.mChangeOfCoord[i];
    }

    if(mNumberOfObservables > 0){
        mObservables = new double[mNumberOfObservables];
        mBackupObservables = new double[mNumberOfObservables];
        for(int i=0; i<mNumberOfObservables; ++ i){
            mObservables[i] = orig.mObservables[i];
            mBackupObservables[i] = orig.mBackupObservables[i];
        }
    }

    for(int i=0; i<NUMBEROFSTATS;++i) stat[i] = orig.stat[i];

    // copy parameters for potential
    for(int i=0; i<NUMBEROFMONOMERTYPES; ++i){
        for(int j=0; j<NUMBEROFMONOMERTYPES; ++j){
            LJ_cutoff[i][j]   =orig.LJ_cutoff[i][j];  
            LJ_sigma[i][j]    =orig.LJ_sigma[i][j];  
            LJ_epsilon[i][j]  =orig.LJ_epsilon[i][j];
            LJ_sigma[i][j]    =orig.LJ_sigma[i][j];  
            LJ_shift[i][j]    =orig.LJ_shift[i][j];  
            FENE_WCA_cutoff[i][j] = orig.FENE_WCA_cutoff[i][j];
        }
    }
    SCR_epsilon     = orig.SCR_epsilon ;     
    SCR_sigma       = orig.SCR_sigma   ;       
    SCR_sigma_9     = orig.SCR_sigma_9 ;     
    SCR_cutoff      = orig.SCR_cutoff  ;      
    SCR_shift       = orig.SCR_shift   ;       
    FENE_K          = orig.FENE_K      ;          
    FENE_R          = orig.FENE_R      ;          
    FENE_R2         = orig.FENE_R2     ;         
    FENE_ro         = orig.FENE_ro     ;         

    return *this;
}

LipidModel::~LipidModel(){
    
    if( mCells              != NULL ) delete [] mCells              ; 
    if( mMonomers           != NULL ) delete [] mMonomers           ; 
    if( mMonomerIndexes     != NULL ) delete [] mMonomerIndexes     ; 
    if( mLipidWater         != NULL ) delete [] mLipidWater         ; 
    if( mChangeOfCoord      != NULL ) delete [] mChangeOfCoord      ; 
    if( mObservables        != NULL ) delete [] mObservables        ; 
    if( mBackupObservables  != NULL ) delete [] mBackupObservables  ; 
    if( mMaxMoveDistance    != NULL ) delete [] mMaxMoveDistance    ; 
    if( mMoveCount          != NULL ) delete [] mMoveCount          ;
    if( mMoveRejectCount    != NULL ) delete [] mMoveRejectCount    ;
    if( mBonds              != NULL ){
        for(int i=0; i<mNumberOfMonomers; ++i)
            delete mBonds[i];
        delete [] mBonds;
    }
    if( mCellID2Coord       != NULL ){
        for(int i=0; i<mNumberOfCells; ++i)
            delete mCellID2Coord[i];
        delete [] mCellID2Coord;
    }
    if(mCoord2CellID        != NULL ){
        for(int i=0; i<mNumberOfCellsPerLine; ++i){
            for(int j=0; j<mNumberOfCellsPerLine; ++j){
                delete [] mCoord2CellID[i][j];
            }
            delete [] mCoord2CellID[i];
        }
        delete [] mCoord2CellID;
    }
    if(mCellCornerCoord != NULL){
        for(int i=0; i<mNumberOfCells; ++i){
            delete mCellCornerCoord[i];
        }
        delete [] mCellCornerCoord;
    }
}

void LipidModel::setOutput(std::ostream& out){
    mOut = &out;
}

void LipidModel::input(std::string filename){
    std::ifstream fin;
    getInputStream(fin,filename);
    input(fin);
    fin.close();
}

void LipidModel::input(std::istream& in){
    std::string line;
    while ( getline(in, line) ) {
        if (line[0] == '#') continue;
        std::istringstream ins;
        ins.str(line);
        std::string word;
        ins >> word;
        if (        word.compare("DIMENSION")               ==0 ){
            ins >> mDim;
        } else if ( word.compare("INITTYPE")                ==0 ){
            ins >> mInitType;
        } else if ( word.compare("NUMBEROFMONOMERS")        ==0 ){
            ins >> mNumberOfMonomers;
            long t = mNumberOfMonomers;
            mOrder = 0;
            do{
                t /= 10;
                mOrder ++;
            }while(t >= 1);
        } else if ( word.compare("NUMBEROFINITLIPIDS")        ==0 ){
            ins >> mNumberOfLipids;
        } else if ( word.compare("MINLENGTHOFCELL")         ==0 ){
            ins >> mMinLengthOfCell;
        } else if ( word.compare("DENSITY")                 ==0 ){
            ins >> mDensity;
        } else if ( word.compare( "MAXMOVEDISTANCE")        ==0 ){
            ins >> mInitMaxMoveDistance >> mAdaptiveMoveFlag;
        } else if ( word.compare("NUMBEROFOBSERVABLES")     ==0 ){
            ins >> mNumberOfObservables;
        } else if ( word.compare("BONDLENGTH")              == 0){
            ins >> mBondLength[0];
            mBondLength[1] = mBondLength[2] = mBondLength[0];
        } else if ( word.compare("BONDLENGTHRANGE")         == 0){
            double bondLengthRange;
            ins >> bondLengthRange;
            mBondLength[1] -= bondLengthRange;
            mBondLength[2] += bondLengthRange;
        } else if ( word.compare("LJ_EPSILON_1")             == 0){
            ins >> LJ_epsilon[0][0] 
                >> LJ_epsilon[0][1] 
                >> LJ_epsilon[0][2];
        } else if ( word.compare("LJ_EPSILON_2")             == 0){
            ins >> LJ_epsilon[1][0] 
                >> LJ_epsilon[1][1] 
                >> LJ_epsilon[1][2];
        } else if ( word.compare("LJ_EPSILON_3")             == 0){
            ins >> LJ_epsilon[2][0] 
                >> LJ_epsilon[2][1] 
                >> LJ_epsilon[2][2];
        } else if ( word.compare("LJ_SIGMA_1")               == 0){
            ins >> LJ_sigma[0][0] 
                >> LJ_sigma[0][1] 
                >> LJ_sigma[0][2];
        } else if ( word.compare("LJ_SIGMA_2")               == 0){
            ins >> LJ_sigma[1][0] 
                >> LJ_sigma[1][1] 
                >> LJ_sigma[1][2];
        } else if ( word.compare("LJ_SIGMA_3")               == 0){
            ins >> LJ_sigma[2][0] 
                >> LJ_sigma[2][1] 
                >> LJ_sigma[2][2];
        } else if ( word.compare("LJ_CUTOFF_1")              == 0){
            ins >> LJ_cutoff[0][0] 
                >> LJ_cutoff[0][1] 
                >> LJ_cutoff[0][2];
        } else if ( word.compare("LJ_CUTOFF_2")              == 0){
            ins >> LJ_cutoff[1][0] 
                >> LJ_cutoff[1][1] 
                >> LJ_cutoff[1][2];
        } else if ( word.compare("LJ_CUTOFF_3")              == 0){
            ins >> LJ_cutoff[2][0] 
                >> LJ_cutoff[2][1] 
                >> LJ_cutoff[2][2];
        } else if ( word.compare("RSC_EPSILON")              == 0){
            ins >> SCR_epsilon;
        } else if ( word.compare("RSC_SIGMA")                == 0){
            ins >> SCR_sigma;
        } else if ( word.compare("RSC_CUTOFF")               == 0){
            ins >> SCR_cutoff;
        } else if ( word.compare("FENE_K")                   == 0){
            ins >> FENE_K;
        } else if ( word.compare("FENE_R")                   == 0){
            ins >> FENE_R;
        } else if ( word.compare("FENE_RO")                  == 0){
            ins >> FENE_ro;
        } else if ( word.compare("MOVEFRACTION")             == 0){
            ins >> mMoveFraction[0] >> mMoveFraction[1] >> mMoveFraction[2];
        }
        ins.clear();
    }
    // find the smallest cutoff
    mLargestCutoff = SCR_cutoff;
    for(int i=0; i<NUMBEROFMONOMERTYPES; ++i){
        for(int j=0; j<NUMBEROFMONOMERTYPES; ++j){
            if(mLargestCutoff < LJ_cutoff[i][j]){
                mLargestCutoff = LJ_cutoff[i][j];
            }
            FENE_WCA_cutoff[i][j] = pow(2, 1.0/6) * LJ_sigma[i][j];
        }
    }
}

void LipidModel::init(std::string initConfFile){
    *mOut  << "     => start Iniliazation of Lipid Model ... \n"; 

    *mOut << "       -> pre-calculating potential parameters ... " << std::endl;
    // Lennard-Jones potential
    for(int i=0; i<NUMBEROFMONOMERTYPES; ++i){
        for(int j=0; j<NUMBEROFMONOMERTYPES; ++j){
            LJ_sigma_6[i][j] = pow(LJ_sigma[i][j], 6.0);
            LJ_shift[i][j] = 4 * LJ_epsilon[i][j] * LJ_sigma_6[i][j] * ( LJ_sigma_6[i][j] / pow(LJ_cutoff[i][j], 12.0) - 1 / pow(LJ_cutoff[i][j], 6.0) );
        }
    }
    
    // Soft Core Repulsive potential
    SCR_sigma_9 = pow(SCR_sigma, 9.0);
    SCR_shift = 4 * SCR_epsilon * SCR_sigma_9 / pow(SCR_cutoff, 9);

    // FENE potential
    FENE_R2 = FENE_R*FENE_R;

    *mOut << "       -> building cells ... " << std::endl;
    buildCells();
    
    initBonds();
    mObservables = new double[mNumberOfObservables];
    mBackupObservables = new double[mNumberOfObservables];
    mChangeOfCoord = new double[mDim];
    mMonomers = new Monomer[mNumberOfMonomers];
    mMonomerIndexes = new MonomerIndex[mNumberOfMonomers];
    mLipidWater = new long[mNumberOfMonomers];
    for(int i=0; i < mNumberOfMonomers; ++i){
        mLipidWater[i] = i;
    }

    long apl = ceil(pow((double)mNumberOfMonomers, 1.0/3)), box_index, atom_index(0);
    double x, y, z, xo, yo, zo, delta, hdelta, hbox;
    
    *mOut << "       -> adding monomers ... \n"
          << "          * monomers per line " << apl << std::endl;
    
    delta = mLengthOfSimBox/(apl+1);
    hdelta = 0.2*(delta/2.0);
    xo = yo = zo = delta/2.0;
    hbox = mLengthOfSimBox*0.5;

    if(!mInitType){
        for(long k=0; k<apl && atom_index < mNumberOfMonomers; ++k){
            for(long j=0; j<apl && atom_index < mNumberOfMonomers; ++j){
                for(long i=0; i<apl && atom_index < mNumberOfMonomers; ++i){

                    x = xo + i*delta + mRandom -> nextDouble(-hdelta, +hdelta) - hbox;
                    y = yo + j*delta + mRandom -> nextDouble(-hdelta, +hdelta) - hbox;
                    z = zo + k*delta + mRandom -> nextDouble(-hdelta, +hdelta) - hbox;
                     
                    Monomer a(x, y, z);
                    applyBoundaryCondition(a);
                    a.mType = WATER;
                    a.mID = atom_index;
                    box_index = whichCell(x, y, z);
                    a.mCellID = box_index;
                    a.mLipidWaterID = a.mID;
                    mMonomers[a.mID] = a;
                    mCells[box_index].addMonomer(a.mID);
                    atom_index++;
                }
            }
        }
        *mOut << "          * created monomers " << atom_index << std::endl;
        
        *mOut << "       -> building monomer indexes ... " << std::endl;
        buildMonomerIndexes();

        *mOut << "       -> creating lipids ... " << std::endl;
        long initNumberOfLipids(mNumberOfLipids);
        mNumberOfLipids = 0;
        if(mMaxNumberOfLipids <= 0) mMaxNumberOfLipids = initNumberOfLipids;
        if(mMinNumberOfLipids <= 0) mMinNumberOfLipids = initNumberOfLipids;
        if(initNumberOfLipids == 0 || initNumberOfLipids > mMaxNumberOfLipids || initNumberOfLipids < mMinNumberOfLipids){
            initNumberOfLipids = mRandom->nextLong( mMinNumberOfLipids, mMaxNumberOfLipids );
        }
        

        *mOut << "          * initial number of lipids: " << initNumberOfLipids << std::endl;
        int l=0;
        while(l < initNumberOfLipids && insertLipid()) l++;
        
        if(mNumberOfLipids != initNumberOfLipids){
            errorMsg("init()", "cannot insert enough lipids!");
            exit(1);
        }

    }else{
        if(initConfFile == ""){
            *mOut << "Error @init(...): the init. conf. file cannot be empty, for init. type 1 !!!\n";
            exit(1);
        }
        std::ifstream fin;
        getInputStream(fin, initConfFile);
        readModel(fin);
        *mOut << "          * file name: "  << initConfFile << "\n";
        *mOut << "          * number of monomers: "     << mNumberOfMonomers << "\n";
        *mOut << "          * number of lipids: "       << mNumberOfLipids << std::endl;
        buildMonomerIndexes();
    }

    *mOut << "       -> check model integrity ... " << std::endl;
    checkModelIntegrity();


    mEnergy = potential();
    doMeasure();

    *mOut << "       -> do observables measurement ... \n"
          << "          * initial observables: \t" ;
    for(int d=0; d<mNumberOfObservables; ++d) 
        *mOut<< mObservables[d] << "\t";

    *mOut << "\n       -> initialize adaptive move arrays ... ";
    if(mAdaptiveMoveFlag)
        *mOut << "\n           * adaptive move flag ON!\n";
    else
        *mOut << "\n           * adaptive move flag OFF!\n";

    initAdapTrialMove();
    *mOut << "     => initialize LipidModel sucessfully !" 
          << std::endl;
}

void LipidModel::buildMonomerIndexes(){
    if (mMonomerIndexes != NULL) delete [] mMonomerIndexes;
    mMonomerIndexes = new MonomerIndex[mNumberOfMonomers];
    for (long i = 0; i < mNumberOfCells; ++i) {
        MonomerIndex iter = mCells[i].mMonomers.begin();
        while (iter != mCells[i].mMonomers.end()) {
            mMonomerIndexes[*iter] = iter;
            iter++;
        }
    }
}

void LipidModel::buildCellHelpers(){
    mCellID2Coord = new long*[mNumberOfCells];
    mCellCornerCoord = new MathVector<double>*[mNumberOfCells];

    if(mNumberOfCellsPerLine == 0 || mNumberOfCellsPerLine == 1){
        mCoord2CellID = new long**[1];
        mCoord2CellID[0] = new long*[1];
        mCoord2CellID[0][0] = new long[1];
        mCoord2CellID[0][0][0] = 0;
    }else if(mNumberOfCellsPerLine > 1){
        mCoord2CellID = new long**[mNumberOfCellsPerLine];
        for(int i=0; i<mNumberOfCellsPerLine; ++i){
            mCoord2CellID[i] = new long*[mNumberOfCellsPerLine];
            for(int j=0; j<mNumberOfCellsPerLine; ++j){
                mCoord2CellID[i][j] = new long[mNumberOfCellsPerLine];
            }
        }
    }

    for(int i=0; i<mNumberOfCells; ++i) {
        mCellID2Coord[i] = new long[mDim];
        mCellCornerCoord[i] = new MathVector<double>(mDim);
        if(mDim == 3){
            mCellID2Coord[i][2] = i / mNumberOfCellsPerLine2;
            mCellID2Coord[i][1] = (i - mCellID2Coord[i][2]*mNumberOfCellsPerLine2)/mNumberOfCellsPerLine;
            mCellID2Coord[i][0] = i - mCellID2Coord[i][1]*mNumberOfCellsPerLine - mCellID2Coord[i][2]*mNumberOfCellsPerLine2; 
            
            mCoord2CellID[ mCellID2Coord[i][0] ][ mCellID2Coord[i][1] ][ mCellID2Coord[i][2] ] = i;
            
            for(int j=0; j<mDim; ++j){
                mCellCornerCoord[i] -> set( j, mCellID2Coord[i][j]*mLengthOfCell );
            }
        }else {
            errorMsg("buildCells", "simulation currently only support 3D system!");
            exit(1);
        }
    }
}

void LipidModel::buildCells(){
    if(mCells != NULL) delete [] mCells;
    mLengthOfSimBox = pow( mNumberOfMonomers / mDensity, 1.0 / 3.0 );    
    mLowerBoundOfSimBox = -mLengthOfSimBox/2;
    mHigherBoundOfSimBox = mLengthOfSimBox/2;
    mNumberOfCellsPerLine = (long)floor(mLengthOfSimBox / mMinLengthOfCell);
    mNumberOfCellsPerLine2 = mNumberOfCellsPerLine*mNumberOfCellsPerLine;
    mNumberOfCells = mNumberOfCellsPerLine2*mNumberOfCellsPerLine;
    
    bool* neighborFlags = new bool [mNumberOfCells];

    if(mNumberOfCellsPerLine == 0 || mNumberOfCellsPerLine == 1){
        mLengthOfCell = mLengthOfSimBox;
        mNumberOfCells = 1;
        mNumberOfCellsPerLine2 = 1;
        mNumberOfCells = 1;
        mCells = new Cell[mNumberOfCells];
    } else {
        mLengthOfCell = mLengthOfSimBox / mNumberOfCellsPerLine;
        mCells = new Cell[mNumberOfCells];

        for (long i = 0; i < mNumberOfCells; ++i) {

            for(int ii=0; ii< mNumberOfCells; ++ii) neighborFlags[ii] = false;

            // assign the neighbor for each box
            long z = i / mNumberOfCellsPerLine2;
            long y = floor(double(i % mNumberOfCellsPerLine2) / mNumberOfCellsPerLine);
            long x = i % mNumberOfCellsPerLine;
            bool flag = true;
            for (long j = 0; flag && j < mDim; ++j) {
                // initialize the neighbor of each box
                if (mDim == 2) flag = false;
                long self   = cellIndex(    x, 
                                            y, 
                                            (z + j - 1 + mNumberOfCellsPerLine) % mNumberOfCellsPerLine);
                long up     = cellIndex(    x, 
                                            (y + 1 + mNumberOfCellsPerLine) % mNumberOfCellsPerLine, 
                                            (z + j - 1 + mNumberOfCellsPerLine) % mNumberOfCellsPerLine);
                long down   = cellIndex(    x,  
                                            (y - 1+mNumberOfCellsPerLine) % mNumberOfCellsPerLine,  
                                            (z + j - 1+ mNumberOfCellsPerLine) % mNumberOfCellsPerLine);
                long left   = cellIndex(    (x - 1 + mNumberOfCellsPerLine) % mNumberOfCellsPerLine,  
                                            y,                                                      
                                            (z + j - 1 + mNumberOfCellsPerLine)% mNumberOfCellsPerLine );
                long right  = cellIndex(    (x + 1 + mNumberOfCellsPerLine) % mNumberOfCellsPerLine,  
                                            y,                                                      
                                            (z + j - 1 + mNumberOfCellsPerLine)% mNumberOfCellsPerLine);
                long upl    = cellIndex(    (x - 1 + mNumberOfCellsPerLine) % mNumberOfCellsPerLine,  
                                            (y + 1 + mNumberOfCellsPerLine) % mNumberOfCellsPerLine,  
                                            (z + j - 1 + mNumberOfCellsPerLine)% mNumberOfCellsPerLine);
                long upr    = cellIndex(    (x + 1 + mNumberOfCellsPerLine) % mNumberOfCellsPerLine,  
                                            (y + 1 + mNumberOfCellsPerLine) % mNumberOfCellsPerLine,  
                                            (z + j - 1 + mNumberOfCellsPerLine)% mNumberOfCellsPerLine);
                long downl  = cellIndex(    (x - 1 + mNumberOfCellsPerLine) % mNumberOfCellsPerLine,  
                                            (y - 1 + mNumberOfCellsPerLine) % mNumberOfCellsPerLine,  
                                            (z + j - 1 + mNumberOfCellsPerLine)% mNumberOfCellsPerLine);
                long downr  = cellIndex(    (x + 1 + mNumberOfCellsPerLine) % mNumberOfCellsPerLine,  
                                            (y - 1 + mNumberOfCellsPerLine) % mNumberOfCellsPerLine,  
                                            (z + j - 1 + mNumberOfCellsPerLine)% mNumberOfCellsPerLine);
                
                if (self != i && self < mNumberOfCells) {
                    neighborFlags[self] = true;
                }
                if (up != i && up < mNumberOfCells) {
                    neighborFlags[up] = true;
                }
                if (down != i && down < mNumberOfCells) {
                    neighborFlags[down] = true;
                }
                if (left != i && left < mNumberOfCells) {
                    neighborFlags[left] = true;
                }
                if (right != i && right < mNumberOfCells) {
                    neighborFlags[right] = true;
                }
                if (upl != i && upl < mNumberOfCells) {
                    neighborFlags[upl] = true;
                }
                if (downl != i && downl < mNumberOfCells) {
                    neighborFlags[downl] = true;
                }
                if (upr != i && upr < mNumberOfCells) {
                    neighborFlags[upr] = true;
                }
                if (downr != i && downr < mNumberOfCells) {
                    neighborFlags[downr] = true;
                }
            }
            for(int ij = 0; ij < mNumberOfCells; ++ij){
                if(neighborFlags[ij])
                    mCells[i].addNeighbor(ij);
            }
        }
    }
    delete [] neighborFlags;


    buildCellHelpers();

    *mOut << "          * length of simulation box: " << mHigherBoundOfSimBox << " - " << mLowerBoundOfSimBox << " = "<< mLengthOfSimBox << "\n"
          << "          * length of cell: " << mLengthOfCell << "\n"
          << "          * number of cells: " << mNumberOfCells << "\n"
          << "          * number of cells per line: " << mNumberOfCellsPerLine << "\n";
}


void LipidModel::readModel(std::ifstream& fin){
    using std::string;
    using std::istringstream;
    using std::list;
    long rNoM(-1),  rNoB(-1), atom_index, box_index;
    double tE, x, y, z;
    int bi, ba, bb, bn;
    char type, waste;
    
    list<long> lipidArray;

    string line, word, readStatus("VARIABLES");

    while(getline(fin, line)){
        istringstream sin(line);
        sin >>  word;
        if(readStatus.compare("VARIABLES") == 0 && word.compare("OLMODEL") == 0){
            getline(fin, line);
            istringstream tsin(line);
            tsin >> rNoM >> rNoB;
            if(rNoM != mNumberOfMonomers){
                errorMsg("readModel", "Number of monomers is incorrect!"); 
                exit(1);
            }
        } else if (readStatus.compare("VARIABLES") == 0 && word.compare("Observables:") == 0){
            sin >> tE >> mNumberOfLipids;
            if(rNoB != mNumberOfLipids*2){
                 errorMsg("readModel", "Number of lipids and bonds are inconsistent!"); 
                 *mOut << "Number of Lipids:  " << mNumberOfLipids << "\n"
                       << "Number of Bonds:   " << rNoB << std::endl;
                 exit(1);
            }
        } else if (word.compare("@<TRIPOS>ATOM") == 0){
            readStatus = "ATOM";
        } else if (word.compare("@<TRIPOS>BOND") == 0){
            readStatus = "BOND";
        } else if (readStatus.compare("ATOM") == 0){
            sin.seekg(0);
            sin >> atom_index >> type >> x >> y >> z >> waste;
            Monomer a(x,y,z); 
            applyBoundaryCondition(a);
            if(type == 'W'){
                a.mType = WATER;
            }
            else if(type == 'H'){
                a.mType = HYDROPHOBIC;
            }
            else if(type == 'P'){
                a.mType = POLAR;
            }
            else{
                errorMsg("readModel", "ilegal monomer type");
                *mOut << "Type : "  << type << std::endl;
                exit(1);
            }
            a.mID = atom_index-1;
            box_index = whichCell(x, y, z);
            a.mCellID = box_index;
            a.mLipidWaterID = a.mID; 
            mMonomers[a.mID] = a;
            mCells[box_index].addMonomer(a.mID);

            if(a.mType != WATER){
                lipidArray.push_back(a.mID);
            }
            if(atom_index == mNumberOfMonomers){
                readStatus = "VARIABLES";
            }

        } else if (readStatus.compare("BOND") == 0){
            sin.seekg(0);
            sin >> bi >> ba >> bb >> bn; 
            addBond(ba-1,bb-1);
        }
    }

    // find all the lipids in the system, and rearrange mLipidWater array
    list<long>::iterator iter;
    list<long>::iterator jter;
    list<long>::iterator kter;
    long head,middle,tail;
    long foundLipids(0), a(0);
    iter = lipidArray.begin();
    while(iter != lipidArray.end()){
        head = -1;
        middle = -1;
        tail = -1;
        jter = iter;
        jter++;
        while(jter != lipidArray.end() && !isBonded(*iter,*jter)) { jter++;}
        
        kter = iter;
        kter++;
        while(kter == jter || (kter != lipidArray.end() && !(isBonded(*iter,*kter) || isBonded(*jter,*kter)) ) ) { kter++;}

        if (jter == lipidArray.end() || kter == lipidArray.end()){
            errorMsg("readModel","lipidArray has problem!");
            exit(1);
        }

        if(mMonomers[*iter].mType == POLAR){
            head = *iter;
            if(isBonded(*iter, *jter)){ middle = *jter; tail = *kter;}
            else { middle = *kter; tail = *jter;}
        } else if(mMonomers[*jter].mType == POLAR){
            head = *jter;
            if(isBonded(*iter, *jter)){ middle = *iter; tail = *kter;}
            else { middle = *kter; tail = *iter;}
        } else if(mMonomers[*kter].mType == POLAR){
            head = *kter;
            if(isBonded(*iter, *kter)){ middle = *iter; tail = *jter;}
            else { middle = *jter; tail = *iter;}
        } else{
            errorMsg("readModel", "error in lipid reconstruction from file!");
            exit(1);
        }
        
        a = foundLipids*3;
        
        swapLipidWaterArray(a  , mMonomers[head].mLipidWaterID  );
        swapLipidWaterArray(a+1, mMonomers[middle].mLipidWaterID);
        swapLipidWaterArray(a+2, mMonomers[tail].mLipidWaterID  );

        foundLipids ++;

        lipidArray.erase(jter); 
        lipidArray.erase(kter); 
        iter = lipidArray.erase(iter); 
    }
    
    if(foundLipids != mNumberOfLipids){
        errorMsg("readModel", "foundLipids not equal to mNumberOfLipids!");
        exit(1);
    }

    mEnergy = potential();
    
    *mOut << "       -> read from file ends\n"
          << "          * Energy From File:  " << tE << "\n"
          << "          * Energy Calculated: " << mEnergy << std::endl;

}

void LipidModel::print(std::ostream& out, PrintType type){
    switch(type){

        case PRINT_MODEL_OBSERVABLES_ONLY:
            out << std::setw(25) << std::setprecision(20) << mObservables[0] << " " 
                << std::setw(8)  << std::setprecision(5)  << mObservables[1] << " " << std::endl;
            break;
        case PRINT_MODEL_OBSERVABLES:
            out << std::setw(25) << std::setprecision(20) << mObservables[0] << " " 
                << std::setw(8)  << std::setprecision(5)  << mObservables[1] << " "
                << std::setw(10) << std::setprecision(8) << (stat[0] == 0 ? 0 : (double) (stat[0] - stat[1])/stat[0]) << " " 
                << std::setw(10) << std::setprecision(8) << (stat[2] == 0 ? 0 : (double) (stat[2] - stat[3])/stat[2]) << " " 
                << std::setw(10) << std::setprecision(8) << (stat[4] == 0 ? 0 : (double) (stat[4] - stat[5])/stat[4]) << " " 
                << std::setw(10) << std::setprecision(8) << (stat[4] == 0 ? 0 : (double) (stat[6] - stat[7])/stat[6]) << " " 
                << std::setw(10) << std::setprecision(8) << (stat[4] == 0 ? 0 : (double) (stat[8] - stat[9])/stat[8]) << " " << std::endl; 

            break;
        case PRINT_MODEL_INFORMATION :
            {

                double mean(0), s(0), mm(0), ma(0);

                mm = mMaxMoveDistance[0];
                ma = mMaxMoveDistance[0];
                for(int i=0; i < mHistogram -> size(); ++i){
                    mean += mMaxMoveDistance[i];
                    if(mMaxMoveDistance[i] < mm){
                        mm = mMaxMoveDistance[i];
                    }
                    if(mMaxMoveDistance[i] > ma){
                        ma = mMaxMoveDistance[i];
                    }
                }
                mean /= mHistogram -> size();
                for(int i=0; i < mHistogram -> size(); ++i){
                    s += pow(mMaxMoveDistance[i] - mean,2.0);
                }
                s /= mHistogram -> size();
                s = sqrt(s);

                out << std::left
                    << "\n"
                    << "***********************************************************" << "\n" 
                    << "*                                                         *" << "\n" 
                    << "*                    MODEL INFORMATION                    *" << "\n" 
                    << "*                                                         *" << "\n" 
                    << "***********************************************************" << "\n" 
                    << std::setw(30) << "Model Dimension "           << std::setw(20) << mDim  << "\n"
                    << std::setw(30) << "Density "                   << std::setw(20) << mDensity << "\n"             
                    << std::setw(30) << "Initialization Type "       << std::setw(20) << mInitType << "\n"             
                    << std::setw(30) << "No. Monomers, Lipids "      << std::setw(15) << mNumberOfMonomers  << " " 
                                                                     << std::setw(15) << mNumberOfLipids << "\n"  
                    << std::setw(30) << "No. Cells, Cells Per Line " << std::setw(15) << mNumberOfCells     
                                                              << " " << std::setw(15) << mNumberOfCellsPerLine << "\n"
                    << std::setw(30) << "Length Box, Cell (Min.) "   << std::setw(15) << mHigherBoundOfSimBox << " - " << mLowerBoundOfSimBox << " = "<< mLengthOfSimBox    
                                                              << " " << std::setw(15) << mLengthOfCell << " ( " << std::setw(8) << mMinLengthOfCell << " )"<< "\n"
                    << std::setw(30) << "Init Max. Move Distance "   << std::setw(20) << mInitMaxMoveDistance   << "\n"
                    << std::setw(30) << "Current Max. Move Distance "<< std::setw(10) << mean << "(" << std::setw(8) << s << "), range [" << std::setw(5) << mm << ", " << std::setw(5)<< ma << "]\n"
                    << std::setw(30) << "Observables "               << std::setw(6) << mNumberOfObservables << ": "
                                                                     << std::setw(15) << (mObservables == NULL? -1 : mObservables[0]) << " " 
                                                                     << std::setw(15) << (mObservables == NULL? -1 : mObservables[1])  << "\n" 
                    << std::setw(30) << "Max. Min. Energy "          << std::setw(15) << mMaxEnergy << " " << std::setw(15) << mMinEnergy << "\n"
                    << std::setw(30) << "Max. Min. No. Lipids "      << std::setw(15) << mMaxNumberOfLipids << " " << std::setw(15) << mMinNumberOfLipids << "\n"
                    << std::setw(30) << "Move Fraction "             << std::setw(15) << mMoveFraction[0] << "\t" 
                                                                     << std::setw(15) << mMoveFraction[1] << "\t" 
                                                                     << std::setw(15) << mMoveFraction[2] << "\n"
                    << std::setw(30) << "Largest Cutoff "            << std::setw(15) << mLargestCutoff   << "\n"
                    << "\n"
                    << "****************************\n" 
                    << "*    L-J Potential Info    *\n" 
                    << "****************************\n" 
                    << "LJ_cutoff " << LJ_cutoff[0][0] << std::setw(5) << " " << LJ_cutoff[0][1] << std::setw(5) << " " << LJ_cutoff[0][2] << " \n"
                    << "LJ_cutoff " << LJ_cutoff[1][0] << std::setw(5) << " " << LJ_cutoff[1][1] << std::setw(5) << " " << LJ_cutoff[1][2] << " \n"
                    << "LJ_cutoff " << LJ_cutoff[2][0] << std::setw(5) << " " << LJ_cutoff[2][1] << std::setw(5) << " " << LJ_cutoff[2][2] << " \n\n"
                    << "LJ_sigma "  << LJ_sigma[0][0]  << std::setw(5) << " " << LJ_sigma[0][1]  << std::setw(5) << " " << LJ_sigma[0][2]  << " \n"
                    << "LJ_sigma "  << LJ_sigma[1][0]  << std::setw(5) << " " << LJ_sigma[1][1]  << std::setw(5) << " " << LJ_sigma[1][2]  << " \n"
                    << "LJ_sigma "  << LJ_sigma[2][0]  << std::setw(5) << " " << LJ_sigma[2][1]  << std::setw(5) << " " << LJ_sigma[2][2]  << " \n\n"
                    << "LJ_epsilon "<< LJ_epsilon[0][0]<< std::setw(5) << " "<< LJ_epsilon[0][1] << std::setw(5) << " "<< LJ_epsilon[0][2] << " \n"
                    << "LJ_epsilon "<< LJ_epsilon[1][0]<< std::setw(5) << " "<< LJ_epsilon[1][1] << std::setw(5) << " "<< LJ_epsilon[1][2] << " \n"
                    << "LJ_epsilon "<< LJ_epsilon[2][0]<< std::setw(5) << " "<< LJ_epsilon[2][1] << std::setw(5) << " "<< LJ_epsilon[2][2] << " \n"
                    << "\n"
                    << "****************************\n" 
                    << "*   FENE Potential Info    *\n" 
                    << "****************************\n" 
                    << "FENE_K " << std::setw(20) << FENE_K << "\n"
                    << "FENE_R " << std::setw(20) << FENE_R << "\n"
                    << "FENE_ro "<< std::setw(20) <<FENE_ro << "\n"
                    << "FENE_WCA_cutoff " << FENE_WCA_cutoff[0][0] << std::setw(5) << " " << FENE_WCA_cutoff[0][1] << std::setw(5) << " " << FENE_WCA_cutoff[0][2] << " \n"
                    << "FENE_WCA_cutoff " << FENE_WCA_cutoff[1][0] << std::setw(5) << " " << FENE_WCA_cutoff[1][1] << std::setw(5) << " " << FENE_WCA_cutoff[1][2] << " \n"
                    << "FENE_WCA_cutoff " << FENE_WCA_cutoff[2][0] << std::setw(5) << " " << FENE_WCA_cutoff[2][1] << std::setw(5) << " " << FENE_WCA_cutoff[2][2] << " \n"
                    << "\n"
                    << "****************************\n" 
                    << "*    RSC Potential Info    *\n" 
                    << "****************************\n" 
                    << "RSC_Epsilon " << std::setw(20)   << SCR_epsilon<< "\n" 
                    << "RSC_Sigma "   << std::setw(20)     << SCR_sigma    << "\n" 
                    << "RSC_Cutoff "  << std::setw(20)    << SCR_cutoff  << "\n"
                    << std::endl;
            }
            print(out, PRINT_MODEL_STATISTICS);
            break;

        case PRINT_MAX_MOVE_DISTANCE:
            for(int i=0; i<mHistogram->size(); ++i){
                out << i << "\t" << mMaxMoveDistance[i] << "\n";
            }
            break;

        case PRINT_MODEL_CELL :
            {
                out << "\n"
                    << "Cell ID to Cartesian Coordinates Mapping: \n";
                for (long i = 0; i<mNumberOfCells; ++i){
                    out << "Cell ID = " << std::setw(4) << i << " (";
                    for(long j=0; j<mDim; ++j){
                        out << std::setw(3) << mCellID2Coord[i][j] << ", ";
                    }
                    out << ")" << std::endl;
                }
                out << "\n"
                    << "Cartesian Coordinates to Cell ID Mapping: \n";
                for(long i=0; i<mNumberOfCellsPerLine; ++i){
                    for(long j=0; j<mNumberOfCellsPerLine; ++j){
                        for(long k=0; k<mNumberOfCellsPerLine; ++k){
                            out << std::setw(3) << i << " " 
                                << std::setw(3) << j << " " 
                                << std::setw(3) << k << " --> " 
                                << std::setw(3) << mCoord2CellID[i][j][k] << endl;
                        }
                    }
                }
                out << "\n"
                    << "Cartesian Coordinates of Lowest Corner: \n";
                for(long i=0; i<mNumberOfCells; ++i){
                    out << "Cell ID = " << std::setw(4) << i 
                        << " --> " 
                        << std::setw(5) << mCellCornerCoord[i]->get(0) << " "
                        << std::setw(5) << mCellCornerCoord[i]->get(1) << " "
                        << std::setw(5) << mCellCornerCoord[i]->get(2) << "\n";
                }

                out << "\n"
                    << "Cell Info.: \n";
                for (long i = 0; i < mNumberOfCells; ++i) {
                    out << "\n"
                        << "################\n"
                        << "# Cell ID " << std::setw(5) << i <<"\n"
                        << "# No. Monomers " << std::setw(5) << mCells[i].mNumberOfMonomers << "\n"
                        << "# Neighbor " << mCells[i].mNumberOfNeighbors << "\n";
                    CellIndex biter = mCells[i].mNeighbors.begin();
                    while(biter != mCells[i].mNeighbors.end()){
                        out << std::setw(4) << *biter << ' ';
                        biter++;
                    }
                    out << std::endl;   
                    MonomerIndex iter = mCells[i].mMonomers.begin();
                    while (iter != mCells[i].mMonomers.end()) {
                        out << "Monomer " 
                                << std::setw(8) << mMonomers[*iter].mID
                                << std::setw(5) << mMonomers[*iter].mDim 
                                << std::setw(5) << mMonomers[*iter].mCellID
                                << std::setw(5) << mMonomers[*iter].mLipidWaterID
                                << std::setw(5) << MTToString(mMonomers[*iter].mType) << "\t"
                                << std::setw(12) << std::scientific << mMonomers[*iter].mCoord->get(0) << " "
                                << std::setw(12) << std::scientific << mMonomers[*iter].mCoord->get(1) << " "
                                << std::setw(12) << std::scientific << mMonomers[*iter].mCoord->get(2) << " "
                                << std::endl;
                        iter++;
                    }
                }
            }
            break;
        case PRINT_MODEL_XYZ :
            {
                out << mNumberOfMonomers << "\n\n";
                for(long i=0; i<mNumberOfMonomers; ++i){
                    out << MTToString(mMonomers[i].mType) << " ";
                    for(long d=0; d < mDim; ++d){
                        out << std::scientific << mMonomers[i].mCoord->get(d) << " " ;
                    }
                    out << std::endl;
                }
            }
            break;
        case PRINT_MODEL_MOL2 :
            {
                long count(0);
                for (long i = 0; i < mNumberOfMonomers; ++i) {
                    for (long j = 0; j < MAXBOND; ++j)
                        if (mBonds[i][j] != NOBOND) count++;
                }
                out << "\n"
                     << "@<TRIPOS>MOLECULE\n"
                     << "OLMODEL\n"
                     << std::setw(5) << mNumberOfMonomers << std::setw(5) << count << "\n"
                     << "LIPID\n"
                     << "Observables: " 
                     << std::setw(15) << std::setprecision(10) << mObservables[0] << " " 
                     << std::setw(8)  << std::setprecision(5)  << mObservables[1] << " " 
                     << "\n"
                     << "@<TRIPOS>ATOM\n";
                for (long i = 0; i < mNumberOfMonomers; ++i) {
                    out << " "
                          << i + 1 << " "
                          << MTToString(mMonomers[i].mType) << " "
                          << mMonomers[i].mCoord->get(0) << " "
                          << mMonomers[i].mCoord->get(1) << " "
                          << mMonomers[i].mCoord->get(2) << " "
                          << MTToString(mMonomers[i].mType) << " \n";
                }
                out << "@<TRIPOS>BOND\n";
                count = 0;
                for (long i = 0; i < mNumberOfMonomers; ++i) {
                    for (long j = 0; j < MAXBOND; ++j) {
                        if (mBonds[i][j] != NOBOND)
                            out << ++count << " " << i + 1 << " " << mBonds[i][j] + 1 << " " << 1 << "\n";
                    }
                }
            }
            break;
        case PRINT_MODEL_STATISTICS:
            {
                out << "\n"
                    << "############################\n"
                    << "Statistics: \n"
                    << std::setw(20) << "Move Type"
                    << std::setw(20) << "Total" << std::setw(20) << "Accepted" << std::setw(20) << "Rejected\n"
                    << std::setw(20) << "RandomMove"
                    << std::setw(20) << stat[0] << std::setw(20) << stat[0] - stat[1] << std::setw(20) << stat[1] << "\n"
                    << std::setw(20) << " (ratio)"
                    << std::setw(20) << (double) stat[0] / stat[0] << std::setw(20) << (double) (stat[0] - stat[1]) / stat[0] << std::setw(20) << (double) stat[1] / stat[0] << "\n"
                    << std::setw(20) << "RandomShiftMove"
                    << std::setw(20) << stat[6] << std::setw(20) << stat[6] - stat[7] << std::setw(20) << stat[7] << "\n"
                    << std::setw(20) << " (ratio)"
                    << std::setw(20) << (double) stat[6] / stat[6] << std::setw(20) << (double) (stat[6] - stat[7]) / stat[6] << std::setw(20) << (double) stat[7] / stat[6] << "\n"
                    << std::setw(20) << "ReptationMove"
                    << std::setw(20) << stat[8] << std::setw(20) << stat[8] - stat[9] << std::setw(20) << stat[9] << "\n"
                    << std::setw(20) << " (ratio)"
                    << std::setw(20) << (double) stat[8] / stat[8] << std::setw(20) << (double) (stat[8] - stat[9]) / stat[8] << std::setw(20) << (double) stat[9] / stat[8] << "\n"
                    << std::setw(20) << "InsertionMove"
                    << std::setw(20) << stat[2] << std::setw(20) << stat[2] - stat[3] << std::setw(20) << stat[3] << "\n"
                    << std::setw(20) << " (ratio)"
                    << std::setw(20) << (double) stat[2] / stat[2] << std::setw(20) << (double) (stat[2] - stat[3]) / stat[2] << std::setw(20) << (double) stat[3] / stat[2] << "\n"
                    << std::setw(20) << "DeletionMove"
                    << std::setw(20) << stat[4] << std::setw(20) << stat[4] - stat[5] << std::setw(20) << stat[5] << "\n"
                    << std::setw(20) << " (ratio)"
                    << std::setw(20) << (double) stat[4] / stat[4] << std::setw(20) << (double) (stat[4] - stat[5]) / stat[4] << std::setw(20) << (double) stat[5] / stat[4] << "\n"
                    << std::endl;
            }
            break;
        case PRINT_MODEL_LIPIDWATER_ARRAY:
            out << "### LipidWater Array: [LipidWaterArray ID, mLipidWaterID, mID, mType] \n\n"; 
            for(int i=0; i<mNumberOfMonomers; ++i){
                out << "[" 
                    << std::setw(mOrder) << i << ", " 
                    << std::setw(mOrder) << mMonomers[mLipidWater[i]].mLipidWaterID << ", " 
                    << std::setw(mOrder) << mMonomers[mLipidWater[i]].mID  << ", "
                    << std::setw(1) << MTToString(mMonomers[mLipidWater[i]].mType)  
                    << "]" << "\t";
                if((i+1) % 3 ==0) out << "\n";
            }
            out << "\n" << std::endl;
            break;
        default:
            errorMsg("Print(...)","Illegal Print Type: " + type );
            exit(1);
    }
}


void LipidModel::initBonds(){
    if(mBonds == NULL){
        mBonds = new long*[mNumberOfMonomers];
        for(int i=0; i < mNumberOfMonomers; ++i){
            mBonds[i] = new long[MAXBOND];
            for(int j=0; j < MAXBOND; ++j) mBonds[i][j] = NOBOND;
        }
    }else{
        for(int i=0; i < mNumberOfMonomers; ++i)
            for(int j=0; j < MAXBOND; ++j) mBonds[i][j] = NOBOND;
    }
}

bool LipidModel::isBonded(const long& a, const long& b) {
    if (a > b) {
        for (long i = 0; i < MAXBOND; ++i) {
            if (mBonds[b][i] == a) {
                return true;
            }
        }
    } else if(a < b) {
        for (long i = 0; i < MAXBOND; ++i) {
            if (mBonds[a][i] == b) {
                return true;
            }
        }
    } else{
        errorMsg("isBonded(...)", "cannot bond itself !" );
        exit(1);
    }
    return false;
}

void LipidModel::addBond(const long& a, const long& b) {
    if (a > b) {
        for (long i = 0; i < MAXBOND; ++i) {
            if (mBonds[b][i] == NOBOND) {
                mBonds[b][i] = a;
                return;
            }
        }
        errorMsg("addBond(...)", "bond already existed!");
        exit(1);
    } else if (a < b) {
        for (long i = 0; i < MAXBOND; ++i) {
            if (mBonds[a][i] == NOBOND) {
                mBonds[a][i] = b;
                return;
            }
        }
        errorMsg("addBond(...)", "bond already existed!");
        exit(1);
    } else{
        errorMsg("addBond(...)", "cannot bond itself !" );
        exit(1);
    }

}

void LipidModel::delBond(const long& a, const long& b) {
    if (a > b) {
        for (long i = 0; i < MAXBOND; ++i) {
            if (mBonds[b][i] == a) {
                mBonds[b][i] = NOBOND;
                return;
            }
        }
        errorMsg("delBond(...)", "bond does not exist! " );
        exit(1);
    } else if (a < b) {
        for (long i = 0; i < MAXBOND; ++i) {
            if (mBonds[a][i] == b) {
                mBonds[a][i] = NOBOND;
                return;
            }
        }
        errorMsg("delBond(...)", "bond does not exist! " );
        exit(1);
    } else{
        errorMsg("delBond(...)", "cannot delbond itself !");
        exit(1);
    }
}

void LipidModel::errorMsg(std::string func, std::string msg, bool isWarning){
    if(!isWarning){
        *mOut << "Error @" << func << ", " << msg << std::endl;
    }else{
        *mOut << "Warning @" << func << ", " << msg << std::endl;
    }
}

void LipidModel::updateTrialMoveDistance(long i){
    double p = 1 - (double)mMoveRejectCount[i]/(double)mMoveCount[i];
    mMaxMoveDistance[i] = mMaxMoveDistance[i]*( log(ADAWL_CONSTANT_A * OPTIMAL_ACCEPT_PROB + ADAWL_CONSTANT_B)/log(ADAWL_CONSTANT_A*p + ADAWL_CONSTANT_B) );
    mMoveRejectCount[i] = 0;
    mMoveCount[i] = 0;


   // *mOut << " update max. move distance for bin " << i << " to " << mMaxMoveDistance[i] << std::endl;
}

void LipidModel::initAdapTrialMove(){
    if( mMaxMoveDistance    != NULL ) delete [] mMaxMoveDistance    ; 
    if( mMoveCount          != NULL ) delete [] mMoveCount          ;
    if( mMoveRejectCount    != NULL ) delete [] mMoveRejectCount    ;
    mMaxMoveDistance = new double[mHistogram->size()];
    mMoveCount = new long[mHistogram->size()];
    mMoveRejectCount = new long[mHistogram->size()];
    for(int i=0; i<mHistogram->size(); ++i){
        mMaxMoveDistance[i] = mInitMaxMoveDistance;
        mMoveCount[i] = 0;
        mMoveRejectCount[i] = 0;
    }
}

void LipidModel::setConstraint(const HistogramND* histogram){

    mHistogram = histogram;
    if( histogram -> dim() == 2 ){
        mMaxEnergy = histogram -> max(0);
        mMinEnergy = histogram -> min(0);
        mMaxNumberOfLipids = (int)(histogram -> max(1) + 0.5);
        mMinNumberOfLipids = (int)(histogram -> min(1) + 0.5);
    } else if( histogram -> dim() ==1 ){
        mMaxEnergy = histogram -> max(0);
        mMinEnergy = histogram -> min(0);
        mMaxNumberOfLipids = -1;
        mMinNumberOfLipids = -1;
    }
}

void LipidModel::setRandom(Random& ran){
    mRandom = &ran;
}

double LipidModel::potentialLJ(const Monomer& m1, const Monomer& m2){
    double r = m1.mCoord -> minImageDist(*m2.mCoord, mLengthOfSimBox);
    if (r >= LJ_cutoff[m1.mType][m2.mType]) return 0; // out of cutoff distance, no longeraction
    double r_6 = pow(r, 6.0);
    double r_12 = r_6*r_6;
    return 4 * LJ_epsilon[m1.mType][m2.mType] * LJ_sigma_6[m1.mType][m2.mType] * ( LJ_sigma_6[m1.mType][m2.mType] / r_12 - 1 / r_6 ) - LJ_shift[m1.mType][m2.mType];
}

double LipidModel::potentialSCR(const Monomer& m1, const Monomer& m2){
    double r = m1.mCoord -> minImageDist(*m2.mCoord, mLengthOfSimBox);
    if (r >= SCR_cutoff) return 0;
    return 4 * SCR_epsilon * SCR_sigma_9 / pow(r, 9) - SCR_shift;
}

double LipidModel::potentialFENE(const Monomer& m1, const Monomer& m2){
    double r = m1.mCoord -> minImageDist(*m2.mCoord, mLengthOfSimBox);
    //double potential_wca(0);
    //if ( r < FENE_WCA_cutoff[m1.mType][m2.mType]){
    //    double r_6 = pow(r, 6.0);
    //    double r_12 = r_6*r_6;
    //    potential_wca = 4 * LJ_epsilon[m1.mType][m2.mType] * LJ_sigma_6[m1.mType][m2.mType] * ( LJ_sigma_6[m1.mType][m2.mType] / r_12 - 1 / r_6 ) + LJ_epsilon[m1.mType][m2.mType];
    //}
    return -FENE_K * FENE_R2 * 0.5 * log( 1 - pow(( (r - FENE_ro) / FENE_R), 2.0) );
}

double LipidModel::potential(const Monomer& m, const Monomer& n){
    double e(0);
    switch(m.mType){
        case WATER:
            switch(n.mType){
                case WATER: e = potentialLJ(m,n);       break;
                case POLAR: e = potentialLJ(m,n);       break;
                case HYDROPHOBIC: e = potentialSCR(m,n);      break;
            }
            break;
        case HYDROPHOBIC:
            switch(n.mType){
                case WATER: e = potentialSCR(m,n);       break;
                case HYDROPHOBIC: e = (isBonded(m.mID,n.mID) ? potentialFENE(m,n) : potentialLJ(m,n));  break;
                case POLAR: e = (isBonded(m.mID,n.mID) ? potentialFENE(m,n) : potentialSCR(m,n)); break;
            }
            break;
        case POLAR: 
            switch(n.mType){
                case WATER: e = potentialLJ(m,n); break;
                case POLAR: e = potentialLJ(m,n); break;
                case HYDROPHOBIC: e = (isBonded(m.mID,n.mID) ? potentialFENE(m,n) : potentialSCR(m,n)); break;
            }
            break;
    }
    if(std::isnan(e)) mIsMoveValid = false;
    return e;
}

double LipidModel::potential(const long& i){
   long mi, mj, mk;
   long ci = mMonomers[i].mCellID, ck;
   double tEnergy(0);
   MonomerIndex iptr = mMonomerIndexes[i];
   MonomerIndex jptr = mCells[ ci ].mMonomers.begin();
   MonomerIndex kptr;
   mi = *iptr;
   while( jptr != mCells[ci].mMonomers.end() ){
       if(jptr != iptr){
            mj = *jptr;
            tEnergy += potential(mMonomers[mi], mMonomers[mj]);
       }
       jptr++;

   }

   CellIndex c_iter = mCells[ci].mNeighbors.begin();
   while( c_iter != mCells[ci].mNeighbors.end() ){
       ck = *c_iter;
       kptr = mCells[ck].mMonomers.begin();
       while( kptr != mCells[ck].mMonomers.end() ){
           mk = *kptr;
           tEnergy += potential(mMonomers[mi], mMonomers[mk]); 

           kptr++;
       }
       c_iter++;
   }
   return tEnergy;
}

double LipidModel::potential(){
    double tEnergy(0);
    long mi, mj, mk, i, k;
    MonomerIndex iptr, jptr, kptr;
    CellIndex c_iter;
    for(i = 0; i < mNumberOfCells; ++i){
        iptr = mCells[i].mMonomers.begin();
        while( iptr != mCells[i].mMonomers.end() ){
            mi = *iptr;
            jptr = iptr;
            jptr++;
            // go through everybody inside the cell
            while( jptr != mCells[i].mMonomers.end() ){
                mj = *jptr;
                tEnergy += potential(mMonomers[mi], mMonomers[mj]);
                jptr++;
            }


            // go through everybody inside the neighboring cells
            c_iter = mCells[i].mNeighbors.begin();
            while( c_iter != mCells[i].mNeighbors.end() ){
                k = *c_iter;
                if(i < k){
                    kptr = mCells[k].mMonomers.begin();
                    while( kptr != mCells[k].mMonomers.end() ){
                        mk = *kptr;
                        tEnergy += potential(mMonomers[mi], mMonomers[mk]); 
                        kptr++;
                    }
                }
                c_iter++;
            }

            iptr++;
        }
    }
    return tEnergy;
}

void LipidModel::doMeasure(){
    mObservables[0] = mEnergy;
    mObservables[1] = (double)mNumberOfLipids;
}

void LipidModel::backupMeasure(){
    for (long i = 0; i < mNumberOfObservables; ++i)
        mBackupObservables[i] = mObservables[i];
}

bool LipidModel::insertLipid(){
    long index[3],  numberOfTrials(0);

    RandomAccessNeighborList nlist;
    do{
        index[0] = mLipidWater[ mRandom->nextLong(mNumberOfLipids*3, mNumberOfMonomers) ];
        if( mMonomers[index[0]].mType != WATER ){
            errorMsg("insertLipid()", "illegal type found !");
            exit(1);
        }
        nlist = findWaterNeighborList(index[0]);
        numberOfTrials++; // incase of dead loop
    } while(nlist.size() < 2 && numberOfTrials < mNumberOfMonomers);

    if(numberOfTrials >= mNumberOfMonomers) {
        return false;
    }

    mNeighborWaterAtInsertion = nlist.size();

    std::vector<long> picks = nlist.randomPick(*mRandom, 2);
    
    index[1] = picks[0];
    index[2] = picks[1];


    insertLipid(index[1], index[0], index[2]);

    return true;
}

void LipidModel::insertLipid(const long& head, const long& middle, const long& tail){
    if( mMonomers[head].mType   != WATER   ||
        mMonomers[middle].mType != WATER   ||
        mMonomers[tail].mType   != WATER      ) {
        errorMsg("insertLipid(...)", "cannot created lipid, not all the monomers are W type!");
        exit(1);
    }
    
    // change the middle
    mEnergy -= potential(middle);
    mMonomers[middle].mType = HYDROPHOBIC;
    mEnergy += potential(middle);

    // change the tail
    mEnergy -= potential(tail);
    mMonomers[tail].mType = HYDROPHOBIC;
    addBond(middle, tail);
    mEnergy += potential(tail);
    
    // change the head
    mEnergy -= potential(head);
    mMonomers[head].mType = POLAR;
    addBond(head, middle);
    mEnergy += potential(head);

    long a = mNumberOfLipids*3;

    mLastSwapIDsInLipidWater[0] = mMonomers[head].mLipidWaterID   ;
    mLastSwapIDsInLipidWater[1] = mMonomers[middle].mLipidWaterID ;
    mLastSwapIDsInLipidWater[2] = mMonomers[tail].mLipidWaterID   ;

    swapLipidWaterArray(a  , mMonomers[head].mLipidWaterID  );
    swapLipidWaterArray(a+1, mMonomers[middle].mLipidWaterID);
    swapLipidWaterArray(a+2, mMonomers[tail].mLipidWaterID  );

    mLastCreatedLipidHeadIDInLipidWater = a;
    mNumberOfLipids++;
}

void LipidModel::undoInsertLipid(){

    long i = mLastCreatedLipidHeadIDInLipidWater;
    mLastDeletedLipid[0] = mLipidWater[i];
    mLastDeletedLipid[1] = mLipidWater[i+1];
    mLastDeletedLipid[2] = mLipidWater[i+2];
    
    // change head monomer
    delBond(mLastDeletedLipid[0], mLastDeletedLipid[1]);
    mMonomers[ mLastDeletedLipid[0] ].mType = WATER;

    // change tail monomer
    mMonomers[ mLastDeletedLipid[2] ].mType = WATER;
    delBond( mLastDeletedLipid[1], mLastDeletedLipid[2]);

    // change middle monomer
    mMonomers[ mLastDeletedLipid[1] ].mType = WATER;

    swapLipidWaterArray(i  , mLastSwapIDsInLipidWater[0]);
    swapLipidWaterArray(i+1, mLastSwapIDsInLipidWater[1]);
    swapLipidWaterArray(i+2, mLastSwapIDsInLipidWater[2]);
    
    mEnergy = mBackupObservables[0];
    mNumberOfLipids--;
    mLastDeletedLipid[0] =-1; 
    mLastDeletedLipid[1] =-1;
    mLastDeletedLipid[2] =-1;
}

bool LipidModel::canBeBonded(const Monomer& m, const Monomer& n){
    double dis = m.mCoord -> minImageDist(*n.mCoord, mLengthOfSimBox);
    return (mBondLength[1] < dis && dis < mBondLength[2]);
}

void LipidModel::deleteLipid(long i){
    if (i <= 0) {
        i = mRandom -> nextLong(1, mNumberOfLipids+1);
    }
    
    i = 3*(i-1);
    mLastDeletedLipid[0] = mLipidWater[i];
    mLastDeletedLipid[1] = mLipidWater[i+1];
    mLastDeletedLipid[2] = mLipidWater[i+2];

    // change head monomer
    mEnergy -= potential( mLastDeletedLipid[0] );
    delBond(mLastDeletedLipid[0], mLastDeletedLipid[1]);
    mMonomers[ mLastDeletedLipid[0] ].mType = WATER;
    mEnergy += potential( mLastDeletedLipid[0] );

    // change tail monomer
    mEnergy -= potential( mLastDeletedLipid[2] );
    mMonomers[ mLastDeletedLipid[2] ].mType = WATER;
    delBond( mLastDeletedLipid[1], mLastDeletedLipid[2]);
    mEnergy += potential( mLastDeletedLipid[2] );

    // change middle monomer
    mEnergy -= potential( mLastDeletedLipid[1] );
    mMonomers[ mLastDeletedLipid[1] ].mType = WATER;
    mEnergy += potential( mLastDeletedLipid[1]);

    RandomAccessNeighborList nlist = findWaterNeighborList(mLastDeletedLipid[1]); 
    mNeighborWaterAtDeletion = nlist.size();

    int a = 3*mNumberOfLipids-3;

    mLastSwapIDsInLipidWater[0] = mMonomers[mLastDeletedLipid[0]].mLipidWaterID;
    mLastSwapIDsInLipidWater[1] = mMonomers[mLastDeletedLipid[1]].mLipidWaterID;
    mLastSwapIDsInLipidWater[2] = mMonomers[mLastDeletedLipid[2]].mLipidWaterID;

    swapLipidWaterArray(a,   mMonomers[mLastDeletedLipid[0]].mLipidWaterID);
    swapLipidWaterArray(a+1, mMonomers[mLastDeletedLipid[1]].mLipidWaterID);
    swapLipidWaterArray(a+2, mMonomers[mLastDeletedLipid[2]].mLipidWaterID);
    
    mNumberOfLipids--;
}

void LipidModel::undoDeleteLipid(){

    long head = mLastDeletedLipid[0];
    long middle = mLastDeletedLipid[1];
    long tail = mLastDeletedLipid[2];

    if( mMonomers[head].mType != WATER      ||
        mMonomers[middle].mType != WATER    ||
        mMonomers[tail].mType != WATER      ) {
        errorMsg("undoDeleteLipid(...)", "cannot created lipid, not all the monomers are W type!");
        exit(1);
    }
    
    // change the middle
    mMonomers[middle].mType = HYDROPHOBIC;
    
    // change the tail
    mMonomers[tail].mType = HYDROPHOBIC;
    addBond(middle, tail);
    
    // change the head
    mMonomers[head].mType = POLAR;
    addBond(head, middle);

    long a = mNumberOfLipids*3;

    swapLipidWaterArray(mMonomers[head].mLipidWaterID,   mLastSwapIDsInLipidWater[0]);
    swapLipidWaterArray(mMonomers[middle].mLipidWaterID, mLastSwapIDsInLipidWater[1]);
    swapLipidWaterArray(mMonomers[tail].mLipidWaterID,   mLastSwapIDsInLipidWater[2]);

    mEnergy = mBackupObservables[0];
    mNumberOfLipids++;
}


void LipidModel::doMCMove(){
    backupMeasure();
    
    mIsMoveValid = true;

    if(mAdaptiveMoveFlag)
        mMovedBinIndex = mHistogram -> bIndex(mObservables);

    if(mAdaptiveMoveFlag && mMovedBinIndex >= 0 && mMoveCount[mMovedBinIndex] == ADAWL_RECALCULATE_MOVE_DISTANCE_THRESHOLD){
        updateTrialMoveDistance(mMovedBinIndex);
    } 
    
    if(mRandom->nextDouble()<0.001){
    mMoveProposal = CHANGEVOLUMEMOVE;
    *mOut << "    --Start ChangeVolumeMOve \n" << std::endl;

    doChangeVolumeMove();
    *mOut << "    ----Length of Simulation box" <<mLengthOfSimBox<<"  \n" << std::endl;
    //*mOut << "    ----Density " <<mDensity<<"  \n" << std::endl;
    //*mOut << "    ---- Product of two " <<mDensity *mLengthOfSimBox <<"  \n" << std::endl;

    } 
    



    double p = mRandom->nextDouble() ;
    if (p < mMoveFraction[0]){
        
        // random local displacement trial move
        mMoveProposal = RANDOMMOVE;

        double r = doRandomMove();
        stat[0]++;

        if(mAdaptiveMoveFlag && mMovedBinIndex >= 0) {
            mMoveCount[mMovedBinIndex] ++;
            if(r > mMaxMoveDistance[mHistogram -> bIndex(mObservables)]){
                mMoveProposal = NONEMOVE;
                mMoveRejectCount[mMovedBinIndex]++;
                undoRandomMove();
                stat[1] ++;
            }
        }

    } else if(p < mMoveFraction[1]){

        // random local shift of lipid
        mMoveProposal = RANDOMSHIFTMOVE;
        double r = doRandomShiftMove();
        stat[6]++;
        if(mAdaptiveMoveFlag && mMovedBinIndex >= 0) {
            mMoveCount[mMovedBinIndex] ++;
            if(r > mMaxMoveDistance[mHistogram -> bIndex(mObservables)]){
                mMoveProposal = NONEMOVE;
                mMoveRejectCount[mMovedBinIndex]++;
                undoRandomShiftMove();
                stat[7] ++;
            }
        }
    
    }else if(p < mMoveFraction[2]){
        if(doReptationMove()){
            mMoveProposal = REPTATIONMOVE;
        }else{
            mMoveProposal = NONEMOVE;
            stat[9]++;
        }
        stat[8]++;
    } else {
        // insertion or deletion trial move
        if ( mRandom->nextDouble() < 0.5 ) {
            // insertion trial move
            stat[2]++;
            if (mNumberOfLipids < mMaxNumberOfLipids) {
                if(insertLipid()){
                    mMoveProposal = INSERTIONMOVE;
                } else{ 
                    mMoveProposal = NONEMOVE;
                    stat[3]++;
                }
            } else {
                mMoveProposal = NONEMOVE;
                stat[3]++;
            }
        } else {
            // deletion trial move
            stat[4]++;
            if (mNumberOfLipids > mMinNumberOfLipids) {
                deleteLipid();
                mMoveProposal = DELETIONMOVE;
            } else {
                mMoveProposal = NONEMOVE;
                stat[5]++;
            }
        }
    }
    doMeasure();
}

void LipidModel::undoMCMove(){
    switch (mMoveProposal) {
        case RANDOMMOVE: 
            undoRandomMove(); 
            if(mAdaptiveMoveFlag && mMovedBinIndex >= 0) mMoveRejectCount[mMovedBinIndex]++;
            stat[1]++;
            break;
        case INSERTIONMOVE: 
            undoInsertLipid();
            stat[3]++;
            break;
        case DELETIONMOVE: 
            undoDeleteLipid();
            stat[5]++;
            break;
        case RANDOMSHIFTMOVE:
            undoRandomShiftMove();
            if(mAdaptiveMoveFlag && mMovedBinIndex >= 0) mMoveRejectCount[mMovedBinIndex]++;
            stat[7]++;
            break;
        case REPTATIONMOVE:
            undoReptationMove();
            stat[9]++;
            break;
        case CHANGEVOLUMEMOVE:
           undoChangeVolumeMove();
           break;
        case NONEMOVE:
            break;
        default:
            errorMsg("undoMCMove", "invalid move type!");
            exit(1);
    }
    doMeasure();
}

bool LipidModel::doReptationMove(long index){
    if(index == -1){
        // randomly picked a lipid
        index = mRandom->nextLong(1, mNumberOfLipids+1);
    }
    
    long l[3] = {-1,-1,-1};
    index = 3*(index-1);
    l[0] = mLipidWater[index];
    l[1] = mLipidWater[index+1];
    l[2] = mLipidWater[index+2];
    
    
    if( mMonomers[l[0]].mType != POLAR || 
        mMonomers[l[1]].mType != HYDROPHOBIC ||
        mMonomers[l[2]].mType != HYDROPHOBIC ){
        errorMsg("doReptation(...)", "selected lipid is not correct!");
    }


    // backup energy and monomers of the lipid 
    mPrevEnergy = mEnergy;
    mLastDeletedLipid[0] = l[0];
    mLastDeletedLipid[1] = l[1];
    mLastDeletedLipid[2] = l[2];
    
    RandomAccessNeighborList nlist, mlist;

    if(mRandom->nextDouble() < 0.5){
        // pick head
        mReptationFlag = true;
        nlist = findWaterNeighborList(l[0]);
        if(nlist.size() < 1){
            return false;         
        }
        mlist = findWaterNeighborList(l[1]);
        
        std::vector<long> picks = nlist.randomPick(*mRandom, 1);
        
        // change tail from H to W
        mEnergy -= potential(l[2]);
        mMonomers[l[2]].mType = WATER;
        delBond(l[1],l[2]);
        mEnergy += potential(l[2]);

        // chagen head from P to H
        mEnergy -= potential(l[0]);
        mMonomers[l[0]].mType = HYDROPHOBIC;
        mEnergy += potential(l[0]);

        // change pick from W to P
        mEnergy -= potential(picks[0]);
        mMonomers[picks[0]].mType = POLAR;
        addBond(picks[0], l[0]);
        mEnergy += potential(picks[0]);

        // reorganize lipid-water array
        swapLipidWaterArray(index+2, index+1);
        swapLipidWaterArray(index+1, index);
        swapLipidWaterArray(index, mMonomers[picks[0]].mLipidWaterID);

    }else{
        // pick tail
        mReptationFlag = false;
        nlist = findWaterNeighborList(l[2]);
        if(nlist.size() < 1){
            return false;         
        }
        mlist = findWaterNeighborList(l[1]);
        
        std::vector<long> picks = nlist.randomPick(*mRandom, 1);

        // change head from P to W
        mEnergy -= potential(l[0]);
        mMonomers[l[0]].mType = WATER;
        delBond(l[1],l[0]);
        mEnergy += potential(l[0]);

        // change middle from H to P
        mEnergy -= potential(l[1]);
        mMonomers[l[1]].mType = POLAR;
        mEnergy += potential(l[1]);

        // change pick from W to H
        mEnergy -= potential(picks[0]);
        mMonomers[picks[0]].mType = HYDROPHOBIC;
        addBond(picks[0], l[2]);
        mEnergy += potential(picks[0]);

        // reorganize lipid-water array
        swapLipidWaterArray(index, index+1);
        swapLipidWaterArray(index+1, index+2);
        swapLipidWaterArray(index+2, mMonomers[picks[0]].mLipidWaterID);
    }
    mReptationCFactor = (double)nlist.size() / (double)(mlist.size() + 1);
    return true;
}

void LipidModel::undoReptationMove(){
    // change energy back
    mEnergy = mBackupObservables[0];

    if(mReptationFlag){

        long curr_lipid_head_lipidwater_id = mMonomers[mLastDeletedLipid[0]].mLipidWaterID - 1; 
        long curr_lipid_head = mLipidWater[curr_lipid_head_lipidwater_id];

        // change back bond
        delBond(curr_lipid_head, mLastDeletedLipid[0]) ;
        addBond(mLastDeletedLipid[1], mLastDeletedLipid[2]);

        // change back type
        mMonomers[curr_lipid_head].mType = WATER;
        mMonomers[mLastDeletedLipid[0]].mType = POLAR;
        mMonomers[mLastDeletedLipid[1]].mType = HYDROPHOBIC;
        mMonomers[mLastDeletedLipid[2]].mType = HYDROPHOBIC;
        
        // change back lipidwater array
        swapLipidWaterArray(curr_lipid_head_lipidwater_id, curr_lipid_head_lipidwater_id + 1);
        swapLipidWaterArray(curr_lipid_head_lipidwater_id + 1, curr_lipid_head_lipidwater_id + 2);
        swapLipidWaterArray(curr_lipid_head_lipidwater_id + 2, mMonomers[mLastDeletedLipid[2]].mLipidWaterID);
    }else{

        long curr_lipid_tail_lipidwater_id = mMonomers[mLastDeletedLipid[2]].mLipidWaterID + 1; 
        long curr_lipid_tail = mLipidWater[curr_lipid_tail_lipidwater_id];

        // change back bond
        delBond(curr_lipid_tail, mLastDeletedLipid[2]);
        addBond(mLastDeletedLipid[0], mLastDeletedLipid[1]);

        // change back type
        mMonomers[curr_lipid_tail].mType = WATER;
        mMonomers[mLastDeletedLipid[0]].mType = POLAR;
        mMonomers[mLastDeletedLipid[1]].mType = HYDROPHOBIC;
        mMonomers[mLastDeletedLipid[2]].mType = HYDROPHOBIC;

        // change back lipidwater array
        swapLipidWaterArray(curr_lipid_tail_lipidwater_id, curr_lipid_tail_lipidwater_id - 1);
        swapLipidWaterArray(curr_lipid_tail_lipidwater_id - 1, curr_lipid_tail_lipidwater_id - 2);
        swapLipidWaterArray(curr_lipid_tail_lipidwater_id - 2, mMonomers[mLastDeletedLipid[0]].mLipidWaterID);
    }
}

void LipidModel::undoRandomMove(){
    long fromCell = mMonomers[mBackupMonomer.mID].mCellID;
    long toCell = mBackupMonomer.mCellID;
    mMonomers[mBackupMonomer.mID] = mBackupMonomer;
    mEnergy = mPrevEnergy;
    migrate(mBackupMonomer.mID, fromCell, toCell);
}

double LipidModel::doRandomMove(long index){
    if(index == -1){
        index = mRandom->nextLong(0,mNumberOfMonomers);
    }
    double r(0);
    for(long i=0; i< mDim; ++i){
        if(!mAdaptiveMoveFlag || mMovedBinIndex < 0)
            mChangeOfCoord[i] = mRandom->nextDouble(-mInitMaxMoveDistance, mInitMaxMoveDistance);
        else
            mChangeOfCoord[i] = mRandom->nextDouble(-mMaxMoveDistance[mMovedBinIndex], mMaxMoveDistance[mMovedBinIndex]);
        r += pow(mChangeOfCoord[i], 2.0);
    }
    
    mPrevEnergy = mEnergy;
    mBackupMonomer = mMonomers[index];
    
    moveMonomer(index, mChangeOfCoord);

   return sqrt(r);
}



void LipidModel::doChangeVolumeMove(){
    long index = 1;
    //1. Generate a random number to decide how much we need to chage the volume(-1,1)
    // If larger than 0 then the system expand, if less than 0 , the system shrinks
    //*mOut << "    ----In doChangeVolumeMOve \n" << std::endl;

    double changeRate = mRandom->nextDouble(-1,1)*mVolumeChangeRate;
        //*mOut << "    ----change Rate is  " << changeRate << "\n" << std::endl;
        *mOut << "    ----changeRate "  << changeRate<< "\n" << std::endl;

    for(long i = 0; i < mNumberOfLipids; ++i){

        long l[3] = {-1,-1,-1};
        index = 3*i;
        l[0] = mLipidWater[index];
        l[1] = mLipidWater[index+1];
        l[2] = mLipidWater[index+2];

        //*mOut << "    ----the l0,l1,l2 is  " << l[0]<<l[1]<<l[2] << "\n"<< std::endl;
        //*mOut << "    ----the type of l0,l1,l2 is  "
        // <<mMonomers[l[0]].mType <<mMonomers[l[1]].mType<<mMonomers[l[2]].mType << "\n"<< std::endl;

        if( mMonomers[l[0]].mType != POLAR ||
            mMonomers[l[1]].mType != HYDROPHOBIC ||
            mMonomers[l[2]].mType != HYDROPHOBIC){
            *mOut << "    ---- selected lipid is not correct!\n" << std::endl;

            errorMsg("doVolumeChangeMove(...)", "selected lipid is not correct!");
        }
    //For every lipid calculate the new position based on the change rate. 

        for(long j=0; j<mDim; j++){

        mChangeOfCoord[j] = mMonomers[l[1]].mCoord->get(j) * changeRate;
        if( i ==5){
            *mOut << "    *** Coord change   "  <<      mChangeOfCoord[j] << "\n" << std::endl;

        }
        }

        mPrevEnergy = mEnergy;
        mBackupList[index] = mMonomers[l[0]];
        mBackupList[index+1] = mMonomers[l[1]];
        mBackupList[index+2] = mMonomers[l[2]];

        if(i ==5){
        *mOut << "    ---- Coord Before "  <<      mMonomers[l[0]].mCoord->get(0) << "\n" << std::endl;

                moveMonomerGroup(l,3,mChangeOfCoord);

        *mOut << "    ----Coord After "  <<  mMonomers[l[0]].mCoord->get(0)  << "\n" << std::endl;

        }
        else{
        moveMonomerGroup(l,3,mChangeOfCoord);}

    }

    mLengthOfSimBox = mLengthOfSimBox * (1+changeRate);
    mDensity = mDensity/(1+changeRate);

}

void LipidModel::undoChangeVolumeMove(){
    for(long j=0; j<mNumberOfLipids; ++j){
        for(int i=0; i<3;i++){
            long fromCell = mMonomers[mBackupList[i+3*j].mID].mCellID;
            long toCell = mBackupList[i+3*j].mCellID;
            mMonomers[mBackupList[i+3*j].mID] = mBackupList[i+3*j];
            migrate(mBackupList[i+3*j].mID, fromCell, toCell);
        }
        mEnergy = mPrevEnergy;
    }
}




double LipidModel::doRandomShiftMove(long index){
    if(index == -1){
        index = mRandom->nextLong(1, mNumberOfLipids+1);
    }
    
    long l[3] = {-1,-1,-1};
    index = 3*(index-1);
    l[0] = mLipidWater[index];
    l[1] = mLipidWater[index+1];
    l[2] = mLipidWater[index+2];
    
    if( mMonomers[l[0]].mType != POLAR || 
        mMonomers[l[1]].mType != HYDROPHOBIC ||
        mMonomers[l[2]].mType != HYDROPHOBIC ){
        errorMsg("doRandomShiftMove(...)", "selected lipid is not correct!");
    }

    double r(0);
    for(long i=0; i< mDim; ++i){
        if(!mAdaptiveMoveFlag || mMovedBinIndex < 0)
            mChangeOfCoord[i] = mRandom->nextDouble(-mInitMaxMoveDistance, mInitMaxMoveDistance);
        else
            mChangeOfCoord[i] = mRandom->nextDouble(-mMaxMoveDistance[mMovedBinIndex], mMaxMoveDistance[mMovedBinIndex]);
        r += pow(mChangeOfCoord[i], 2.0);
    }
    
    mPrevEnergy = mEnergy;
    mBackupLipid[0] = mMonomers[l[0]];
    mBackupLipid[1] = mMonomers[l[1]];
    mBackupLipid[2] = mMonomers[l[2]];
    moveMonomerGroup(l, 3, mChangeOfCoord);

    return sqrt(r);
}

void LipidModel::undoRandomShiftMove(){
    for(int i=0; i<3; ++i){
        long fromCell = mMonomers[mBackupLipid[i].mID].mCellID;
        long toCell = mBackupLipid[i].mCellID;
        mMonomers[mBackupLipid[i].mID] = mBackupLipid[i];
        migrate(mBackupLipid[i].mID, fromCell, toCell);
    }
    mEnergy = mPrevEnergy;
}

void LipidModel::moveMonomerGroup(const long* index, const long& size, double* coord){
    for(int i=0; i<size; ++i)
        moveMonomer(index[i], coord);
}

void LipidModel::moveMonomer(const long& index, double* coord){
    
    double oldEnergy = potential(index);
    long oldCellID = mMonomers[index].mCellID;
    for (long i = 0; i < mDim; ++i) {
        mMonomers[index].mCoord->add(i, coord[i]);
    }
    applyBoundaryCondition(mMonomers[index]);
    long newCellID = whichCell(mMonomers[index]);
    migrate(index, oldCellID, newCellID);
    double newEnergy = potential(index);
    mEnergy += (newEnergy - oldEnergy);

}

void LipidModel::migrate(const long& index, const long& oldB, const long& newB) {
    if (oldB != newB) {
        mCells[newB].mMonomers.push_front(index);
        mCells[newB].mNumberOfMonomers++;
        mCells[oldB].delMonomer(mMonomerIndexes[index]);
        mMonomers[index].mCellID = newB;
        mMonomerIndexes[index] = mCells[newB].mMonomers.begin();
    }
}


double LipidModel::cFactor(const double* obs_up, const double* obs_down){
    long up_i = mHistogram -> bIndex(obs_up);
    long down_i = mHistogram -> bIndex(obs_down);
    double rE  = mMaxMoveDistance[up_i];
    double rEp = mMaxMoveDistance[down_i];
    return pow(rE/rEp, 3.0);
}

double LipidModel::cFactor() {
    double ratio(1);
    switch (mMoveProposal) {

        case NONEMOVE:
            ratio = 0;
            break;

        case RANDOMMOVE: // random move 
            if(mAdaptiveMoveFlag && mMovedBinIndex >= 0){
                double rE  = mMaxMoveDistance[mMovedBinIndex];
                double rEp = mMaxMoveDistance[mHistogram -> bIndex(mObservables)];
                ratio = pow(rE/rEp, 3.0);
            } else ratio =1;
            break;

        case RANDOMSHIFTMOVE:
            if(mAdaptiveMoveFlag && mMovedBinIndex >= 0){
                double rE  = mMaxMoveDistance[mMovedBinIndex];
                double rEp = mMaxMoveDistance[mHistogram -> bIndex(mObservables)];
                ratio = pow(rE/rEp, 3.0);
            } else ratio =1;
            break;

        case INSERTIONMOVE: // creation
            if(mNeighborWaterAtInsertion <= 0){
                errorMsg("cFactor()", " INSERTION -- zero neighbor waters !");
                exit(1);
            }
            ratio = (mNumberOfMonomers - 3 * (mNumberOfLipids - 1)) * mNeighborWaterAtInsertion * (mNeighborWaterAtInsertion - 1);
            ratio /= (double) mNumberOfLipids;

            mNeighborWaterAtInsertion = -1;
            break;

        case DELETIONMOVE: // deletion
            if(mNeighborWaterAtDeletion <= 0){
                errorMsg("cFactor()", " DELETIONMOVE -- zero neighbor waters !");
                exit(1);
            }
            ratio = mNumberOfLipids + 1;
            ratio /= (double) (mNumberOfMonomers - 3 * (mNumberOfLipids)) * mNeighborWaterAtDeletion * (mNeighborWaterAtDeletion - 1);

            mNeighborWaterAtDeletion = -1;
            break;

        case REPTATIONMOVE:
            ratio = mReptationCFactor;
            break;

        case CHANGEVOLUMEMOVE:
            ratio = pow(1-mVolumeChangeRate, mNumberOfMonomers);
            break;
        default:
            errorMsg("cFacotr()", "invalid move type!");
            exit(1);
    }
    return ratio;
}

void LipidModel::checkModelIntegrity(){
    int i(0),ii(0), ij(0), ik(0);

    // check mLipidWater Array related quantities
    for(i=0; i<mNumberOfLipids; ++i){
        ii = 3*i;
        ij = ii+1;
        ik = ii+2;
        if(mMonomers[mLipidWater[ii]].mType != POLAR ||
           mMonomers[mLipidWater[ij]].mType != HYDROPHOBIC ||
           mMonomers[mLipidWater[ik]].mType != HYDROPHOBIC){
            
            errorMsg("checkModelIntegrity", "Lipid portion is illegal!");
            exit(1);
        }
        if( mMonomers[mLipidWater[ii]].mLipidWaterID != ii ||
            mMonomers[mLipidWater[ij]].mLipidWaterID != ij ||
            mMonomers[mLipidWater[ik]].mLipidWaterID != ik){

            errorMsg("checkModelIntegrity", "mLipidWaterID is incorrect!");
            exit(1);
        }

        if( mMonomers[mLipidWater[ii]].mID != mLipidWater[ii] ||
            mMonomers[mLipidWater[ij]].mID != mLipidWater[ij] ||
            mMonomers[mLipidWater[ik]].mID != mLipidWater[ik]){
            
            errorMsg("checkModelIntegrity", "mID is incorrect!");
            exit(1);
        }
    }
    for(i=3*mNumberOfLipids; i< mNumberOfMonomers; ++i){
        if(mMonomers[mLipidWater[i]].mType != WATER){
            errorMsg("checkModelIntegrity", "Water portion is illegal!");
            exit(1);
        }

        if( mMonomers[mLipidWater[i]].mLipidWaterID != i){
            errorMsg("checkModelIntegrity", "mLipidWaterID is incorrect!");
            exit(1);
        }

        if(mMonomers[mLipidWater[i]].mID != mLipidWater[i]){
            errorMsg("checkModelIntegrity", "mID is incorrect!");
            exit(1);
        }
    }

    // check monomer indexes
    for(int i=0; i<mNumberOfMonomers; ++i){
        if(mMonomers[*mMonomerIndexes[i]].mID != i){
            errorMsg("checkModelIntegrity", "mMonomerIndexes are incorrect!");
            exit(1);
        }
    }


    // check monomer cells
    int ci(0);
    for(int i=0; i<mNumberOfMonomers; ++i){
        ci = whichCell(mMonomers[i]);
        if(ci != mMonomers[i].mCellID){
            errorMsg("checkModelIntegrity", "mMonomers.mCellID is incorrect!");
            exit(1);
        }
    }

    // check model constraints
    if(mMaxNumberOfLipids >=0 && mMinNumberOfLipids >=0){
       if( !( mMinNumberOfLipids <= mNumberOfLipids && mNumberOfLipids <= mMaxNumberOfLipids ) ){
            errorMsg("checkModelIntegrity", "mNumberOfLipids is out of limits!");
            exit(1);
       }
    }

}

void LipidModel::backup(std::ofstream& fout){
    long sod = sizeof(double);
    long sol = sizeof(long);
    
    // create temporary array for bonds
    long* tBonds = new long[mNumberOfMonomers*MAXBOND];
    for(int i=0; i<mNumberOfMonomers; ++i){
        for(int j=0; j<MAXBOND; ++j){
            tBonds[i*MAXBOND+j] = mBonds[i][j];
        }
    }

    // create temporary arrays for L-J potential parameters
    long numberOfLJs = NUMBEROFMONOMERTYPES * NUMBEROFMONOMERTYPES;
    double* tLJ_cutoff = new double[numberOfLJs];
    double* tLJ_sigma = new double[numberOfLJs];
    double* tLJ_epsilon = new double[numberOfLJs]; 
    double* tFENE_WCA_cutoff = new double[numberOfLJs];
    for(int i=0; i<NUMBEROFMONOMERTYPES; ++i){
        for(int j=0; j<NUMBEROFMONOMERTYPES; ++j){
            tLJ_cutoff[i*NUMBEROFMONOMERTYPES+j]        = LJ_cutoff[i][j];
            tLJ_sigma[i*NUMBEROFMONOMERTYPES+j]         = LJ_sigma[i][j];
            tLJ_epsilon[i*NUMBEROFMONOMERTYPES+j]       = LJ_epsilon[i][j];
            tFENE_WCA_cutoff[i*NUMBEROFMONOMERTYPES+j]  = FENE_WCA_cutoff[i][j];
        }
    }

    if(mHistogram == NULL){
        errorMsg("backup(...)","histogram cannot be NULL!");
        exit(1);
    }
    
    long numOfBins = 0;
    numOfBins = mHistogram -> size();

    fout.write(reinterpret_cast<char*>(&numOfBins                             ), sol                  );
    fout.write(reinterpret_cast<char*>(mMaxMoveDistance                       ), numOfBins*sod        );
    fout.write(reinterpret_cast<char*>(mMoveCount                             ), numOfBins*sol        );
    fout.write(reinterpret_cast<char*>(mMoveRejectCount                       ), numOfBins*sol        );

    fout.write(reinterpret_cast<char*>(stat                                   ), NUMBEROFSTATS*sod              );
    fout.write(reinterpret_cast<char*>(&mInitType                             ), sizeof(bool)                   );
    fout.write(reinterpret_cast<char*>(&mAdaptiveMoveFlag                     ), sizeof(bool)                   );
    fout.write(reinterpret_cast<char*>(&mNumberOfMonomers                     ), sol                            );
    fout.write(reinterpret_cast<char*>(&mOrder                                ), sol                            );
    fout.write(reinterpret_cast<char*>(&mNumberOfObservables                  ), sol                            );
    fout.write(reinterpret_cast<char*>(&mDim                                  ), sol                            );
    fout.write(reinterpret_cast<char*>(mLastDeletedLipid                      ), 3*sol                          );
    fout.write(reinterpret_cast<char*>(&mMaxNumberOfLipids                    ), sol                            );
    fout.write(reinterpret_cast<char*>(&mMinNumberOfLipids                    ), sol                            );
    fout.write(reinterpret_cast<char*>(&mNumberOfLipids                       ), sol                            );
    fout.write(reinterpret_cast<char*>(&mLastCreatedLipidHeadIDInLipidWater   ), sol                            );
    fout.write(reinterpret_cast<char*>(mLipidWater                            ), mNumberOfMonomers*sol          );
    fout.write(reinterpret_cast<char*>(tBonds                                 ), mNumberOfMonomers*MAXBOND*sol  );
    fout.write(reinterpret_cast<char*>(mMoveFraction                          ), 3*sod                          );
    fout.write(reinterpret_cast<char*>(&mMaxEnergy                            ), sod                            );
    fout.write(reinterpret_cast<char*>(&mMinEnergy                            ), sod                            );
    fout.write(reinterpret_cast<char*>(mBondLength                            ), 3*sod                          );
    fout.write(reinterpret_cast<char*>(&mEnergy                               ), sod                            );
    fout.write(reinterpret_cast<char*>(&mPrevEnergy                           ), sod                            );
    fout.write(reinterpret_cast<char*>(&mDensity                              ), sod                            );
    fout.write(reinterpret_cast<char*>(&mMinLengthOfCell                      ), sod                            );
    fout.write(reinterpret_cast<char*>(&mInitMaxMoveDistance                      ), sod                            );
    fout.write(reinterpret_cast<char*>(mChangeOfCoord                         ), mDim*sod                       );
    fout.write(reinterpret_cast<char*>(mObservables                           ), mNumberOfObservables*sod       );
    fout.write(reinterpret_cast<char*>(mBackupObservables                     ), mNumberOfObservables*sod       );
    fout.write(reinterpret_cast<char*>(&mLargestCutoff                         ), sod                            );
    fout.write(reinterpret_cast<char*>(tLJ_cutoff                             ), numberOfLJs*sod                );
    fout.write(reinterpret_cast<char*>(tFENE_WCA_cutoff                       ), numberOfLJs*sod                );
    fout.write(reinterpret_cast<char*>(tLJ_sigma                              ), numberOfLJs*sod                );
    fout.write(reinterpret_cast<char*>(tLJ_epsilon                            ), numberOfLJs*sod                );
    fout.write(reinterpret_cast<char*>(&SCR_epsilon                           ), sod                            );
    fout.write(reinterpret_cast<char*>(&SCR_sigma                             ), sod                            );
    fout.write(reinterpret_cast<char*>(&SCR_cutoff                            ), sod                            );
    fout.write(reinterpret_cast<char*>(&FENE_K                                ), sod                            );
    fout.write(reinterpret_cast<char*>(&FENE_R                                ), sod                            );
    fout.write(reinterpret_cast<char*>(&FENE_ro                               ), sod                            );
    fout.write(reinterpret_cast<char*>(&mMoveProposal                         ), sizeof(TrialMoveType)          );
    fout.write(reinterpret_cast<char*>(&mBoundaryCondition                    ), sizeof(BoundaryConditionType)  );
    for(int i=0; i<mNumberOfCells;      ++i)     mCells[i].backup(fout);
    for(int i=0; i<mNumberOfMonomers;   ++i)     mMonomers[i].backup(fout);
  
    //    mBackupMonomer.backup(fout);

    delete [] tBonds;
    delete [] tLJ_cutoff;
    delete [] tFENE_WCA_cutoff;
    delete [] tLJ_sigma;
    delete [] tLJ_epsilon;
}

void LipidModel::recover(std::ifstream& fin){

    *mOut << "     => LipidModel: start recovery phase ... \n"
          << "       -> reading from file ... " << std::endl;
    long sod = sizeof(double);
    long sol = sizeof(long);

    long numOfBins;
    fin.read(reinterpret_cast<char*>(&numOfBins                             ), sol                  );
    
    mMaxMoveDistance = new double[numOfBins];
    mMoveCount = new long[numOfBins];
    mMoveRejectCount = new long[numOfBins];

    fin.read(reinterpret_cast<char*>(mMaxMoveDistance                       ), numOfBins*sod        );
    fin.read(reinterpret_cast<char*>(mMoveCount                             ), numOfBins*sol        );
    fin.read(reinterpret_cast<char*>(mMoveRejectCount                       ), numOfBins*sol        );

    fin.read(reinterpret_cast<char*>(stat                                   ), NUMBEROFSTATS*sod                );
    fin.read(reinterpret_cast<char*>(&mInitType                             ), sizeof(bool)                     );
    fin.read(reinterpret_cast<char*>(&mAdaptiveMoveFlag                     ), sizeof(bool)                     );
    fin.read(reinterpret_cast<char*>(&mNumberOfMonomers                     ), sol                              );
    fin.read(reinterpret_cast<char*>(&mOrder                                ), sol                              );
    fin.read(reinterpret_cast<char*>(&mNumberOfObservables                  ), sol                              );
    fin.read(reinterpret_cast<char*>(&mDim                                  ), sol                              );
    fin.read(reinterpret_cast<char*>(mLastDeletedLipid                      ), 3*sol                            );
    fin.read(reinterpret_cast<char*>(&mMaxNumberOfLipids                    ), sol                              );
    fin.read(reinterpret_cast<char*>(&mMinNumberOfLipids                    ), sol                              );
    fin.read(reinterpret_cast<char*>(&mNumberOfLipids                       ), sol                              );
    fin.read(reinterpret_cast<char*>(&mLastCreatedLipidHeadIDInLipidWater   ), sol                              );
    
    if(mLipidWater != NULL) delete [] mLipidWater;
    mLipidWater = new long[mNumberOfMonomers];
    
    long* tBonds = new long[mNumberOfMonomers*MAXBOND];
    
    if(mChangeOfCoord != NULL) delete [] mChangeOfCoord;
    mChangeOfCoord = new double[mDim];
    
    if(mObservables != NULL) delete [] mObservables;
    mObservables = new double[mNumberOfObservables];
    
    if(mBackupObservables != NULL) delete [] mBackupObservables;
    mBackupObservables = new double[mNumberOfObservables];
    
    long numberOfLJs = NUMBEROFMONOMERTYPES * NUMBEROFMONOMERTYPES;
    double* tLJ_cutoff = new double[numberOfLJs];
    double* tLJ_sigma = new double[numberOfLJs];
    double* tLJ_epsilon = new double[numberOfLJs]; 
    double* tFENE_WCA_cutoff = new double[numberOfLJs];
    
    
    fin.read(reinterpret_cast<char*>(mLipidWater                            ), mNumberOfMonomers*sol            );
    fin.read(reinterpret_cast<char*>(tBonds                                 ), mNumberOfMonomers*MAXBOND*sol    );
    fin.read(reinterpret_cast<char*>(mMoveFraction                          ), 3*sod                            );
    fin.read(reinterpret_cast<char*>(&mMaxEnergy                            ), sod                              );
    fin.read(reinterpret_cast<char*>(&mMinEnergy                            ), sod                              );
    fin.read(reinterpret_cast<char*>(mBondLength                            ), 3*sod                            );
    fin.read(reinterpret_cast<char*>(&mEnergy                               ), sod                              );
    fin.read(reinterpret_cast<char*>(&mPrevEnergy                           ), sod                              );
    fin.read(reinterpret_cast<char*>(&mDensity                              ), sod                              );
    fin.read(reinterpret_cast<char*>(&mMinLengthOfCell                      ), sod                              );
    fin.read(reinterpret_cast<char*>(&mInitMaxMoveDistance                      ), sod                              );
    fin.read(reinterpret_cast<char*>(mChangeOfCoord                         ), mDim*sod                         );
    fin.read(reinterpret_cast<char*>(mObservables                           ), mNumberOfObservables*sod         );
    fin.read(reinterpret_cast<char*>(mBackupObservables                     ), mNumberOfObservables*sod         );
    fin.read(reinterpret_cast<char*>(&mLargestCutoff                         ), sod                              );
    fin.read(reinterpret_cast<char*>(tLJ_cutoff                             ), numberOfLJs*sod                  );
    fin.read(reinterpret_cast<char*>(tFENE_WCA_cutoff                       ), numberOfLJs*sod                  );
    fin.read(reinterpret_cast<char*>(tLJ_sigma                              ), numberOfLJs*sod                  );
    fin.read(reinterpret_cast<char*>(tLJ_epsilon                            ), numberOfLJs*sod                  );
    fin.read(reinterpret_cast<char*>(&SCR_epsilon                           ), sod                              );
    fin.read(reinterpret_cast<char*>(&SCR_sigma                             ), sod                              );
    fin.read(reinterpret_cast<char*>(&SCR_cutoff                            ), sod                              );
    fin.read(reinterpret_cast<char*>(&FENE_K                                ), sod                              );
    fin.read(reinterpret_cast<char*>(&FENE_R                                ), sod                              );
    fin.read(reinterpret_cast<char*>(&FENE_ro                               ), sod                              );
    fin.read(reinterpret_cast<char*>(&mMoveProposal                         ), sizeof(TrialMoveType)            );
    fin.read(reinterpret_cast<char*>(&mBoundaryCondition                    ), sizeof(BoundaryConditionType)    );
   
    mLengthOfSimBox = pow( mNumberOfMonomers / mDensity, 1.0 / 3.0 );    
    mLowerBoundOfSimBox = -mLengthOfSimBox/2;
    mHigherBoundOfSimBox = mLengthOfSimBox/2;
    mNumberOfCellsPerLine = (long)floor(mLengthOfSimBox / mMinLengthOfCell);
    mNumberOfCellsPerLine2 = mNumberOfCellsPerLine*mNumberOfCellsPerLine;
    mNumberOfCells = mNumberOfCellsPerLine2*mNumberOfCellsPerLine;
    
    if(mNumberOfCellsPerLine == 0 || mNumberOfCellsPerLine == 1){
        mLengthOfCell = mLengthOfSimBox;
        mNumberOfCells = 1;
        mNumberOfCellsPerLine2 = 1;
        mNumberOfCells = 1;
    } else {
        mLengthOfCell = mLengthOfSimBox / mNumberOfCellsPerLine;
    }
    
    if(mCells != NULL) delete [] mCells;
    mCells = new Cell[mNumberOfCells];
    if(mMonomers != NULL) delete [] mMonomers;
    mMonomers = new Monomer[mNumberOfMonomers];
    *mOut << "       -> rebuilding cells " << std::endl;
    for(int i=0; i<mNumberOfCells   ;   ++i)  mCells[i].recover(fin);
    *mOut << "       -> rebuilding monomers " << std::endl;
    for(int i=0; i<mNumberOfMonomers;   ++i)  mMonomers[i].recover(fin);

    buildCellHelpers();
    //  mBackupMonomer.recover(fin);
    
    // recover bonds
    initBonds(); 
    for(int i=0; i<mNumberOfMonomers; ++i){
        for(int j=0; j<MAXBOND; ++j){
            mBonds[i][j] = tBonds[i*MAXBOND+j];
        }
    }
   
    // recover potentials
    for(int i=0; i<NUMBEROFMONOMERTYPES; ++i){
        for(int j=0; j<NUMBEROFMONOMERTYPES; ++j){
            LJ_cutoff[i][j]         = tLJ_cutoff[i*NUMBEROFMONOMERTYPES+j];
            FENE_WCA_cutoff[i][j]   = tFENE_WCA_cutoff[i*NUMBEROFMONOMERTYPES+j];
            LJ_sigma[i][j]          = tLJ_sigma[i*NUMBEROFMONOMERTYPES+j];
            LJ_epsilon[i][j]        = tLJ_epsilon[i*NUMBEROFMONOMERTYPES+j];
        }
    }
    delete [] tBonds;
    delete [] tLJ_cutoff;
    delete [] tFENE_WCA_cutoff;
    delete [] tLJ_sigma;
    delete [] tLJ_epsilon;
    
    *mOut << "       -> building monomer indexes ... " << std::endl;
    buildMonomerIndexes();


    *mOut << "       -> pre-calculating potential parameters ... \n";
    // Lennard-Jones potential
    for(int i=0; i<NUMBEROFMONOMERTYPES; ++i){
        for(int j=0; j<NUMBEROFMONOMERTYPES; ++j){
            LJ_sigma_6[i][j] = pow(LJ_sigma[i][j], 6.0);
            LJ_shift[i][j] = 4 * LJ_epsilon[i][j] * LJ_sigma_6[i][j] * ( LJ_sigma_6[i][j] / pow(LJ_cutoff[i][j], 12.0) - 1 / pow(LJ_cutoff[i][j], 6.0) );
        }
    }
    
    // Soft Core Repulsive potential
    SCR_sigma_9 = pow(SCR_sigma, 9.0);
    SCR_shift = 4 * SCR_epsilon * SCR_sigma_9 / pow(SCR_cutoff, 9);

    // FENE potential
    FENE_R2 = FENE_R*FENE_R;

    *mOut << "       -> check model integrity ... \n";
    checkModelIntegrity();

    double tEnergy[2] = {0,0};
    for(int i=0; i<mNumberOfMonomers; ++i) tEnergy[0]+= potential(i);

    tEnergy[1] = potential();
    *mOut << "       *  Observables   (recover)  : " << std::setprecision(20)  << mObservables[0] << "\t\t" << mObservables[1] << "\n"
          << "       *  Backup energy (recover)  : " << std::setprecision(20)  << mEnergy << "\n" 
          << "       *  Recalculated energy 1    : " << std::setprecision(20)  << tEnergy[0]/2 << "\n"
          << "       *  Recalculated energy 2    : " << std::setprecision(20)  << tEnergy[1]   << std::endl;
    if(abs(tEnergy[0]/2 - mEnergy) > DOUBLE_TOLERANCE || abs(tEnergy[1] - mEnergy) > DOUBLE_TOLERANCE ){
        errorMsg("recover", "energy difference excess DOUBLE TOLERANCE!");
        exit(1);
    }
    *mOut << "     => LipidModel: Recover phase finish successfully!" << std::endl;
}

 
void LipidModel::packup(long* longData, double* doubleData){
    int base(0);
    longData[0] = mNumberOfMonomers;
    longData[1] = mNumberOfLipids;
    
    base = 2;
    for(int i=0; i<mNumberOfMonomers; ++i)
        longData[base+i] = (long)mMonomers[i].mType; 
    
    base += mNumberOfMonomers;
    for(int i=0; i<mNumberOfMonomers; ++i)
        longData[base+i] = mLipidWater[i];

    base = 0;
    doubleData[0] = mEnergy;
    base = 1;
    for(int i=0; i<mNumberOfMonomers; ++i){
        for(int d=0; d < mDim; ++d){
            doubleData[base + i*mDim + d] = mMonomers[i].mCoord -> get(d);
        }
    }
}

void LipidModel::unpack(long* longData, double* doubleData){
    
    *mOut << "     => Unpack package and rebuild the model ... " << std::endl;

    // check constraints first.
    if(    (   (           mMaxNumberOfLipids > 0 && mMinNumberOfLipids > 0          )
             && ( longData[1] > mMaxNumberOfLipids || longData[1] < mMinNumberOfLipids )
            )
        || ( doubleData[0] > mMaxEnergy || doubleData[0] < mMinEnergy)) {
        
        errorMsg("LipidModel::unpack", "replica exceeds the model constraints! Give up unpacking !", true);
        return;
    }

    int base(0);
    if( mNumberOfMonomers != longData[0] ) {
        errorMsg("LipidModel::unpack", "mNumberOfMonomers is inconsistent ! ");
        exit(1);
    }
    mNumberOfLipids = longData[1];
    base = 2;
    *mOut << "       ->  reassign type to monomers ... " << std::endl;
    // reassign type information to monomers
    for(int i=0; i<mNumberOfMonomers; ++i)
        mMonomers[i].mType = (MonomerType)longData[base+i];
    base += mNumberOfMonomers;
    for(int i=0; i<mNumberOfMonomers; ++i){
        mLipidWater[i] = longData[base+i];
    
        // reassign the mLipidWaterID
        mMonomers[mLipidWater[i]].mLipidWaterID = i;
    }

    *mOut << "       ->  rebuild bonds information ... " << std::endl;
    // build up bonds information
    initBonds();
    int l;
    for(int i=0; i<mNumberOfLipids; ++i){
        l = 3*i;
        addBond(mLipidWater[l]   , mLipidWater[l+1]);
        addBond(mLipidWater[l+1] , mLipidWater[l+2]);
    }

    // clear monomers information inside cells
    for(int i=0; i<mNumberOfCells; ++i){
        mCells[i].mMonomers.clear();
        mCells[i].mNumberOfMonomers = 0;
    }


    *mOut << "       ->  reassign cooridnates to monomers and then monomers to cells ... " << std::endl;
    mEnergy = doubleData[0];
    base = 1;
    for(int i=0; i<mNumberOfMonomers; ++i){
        // reassign coordinates
        for(int d=0; d<mDim; ++d){
            mMonomers[i].mCoord -> set(d, doubleData[base + i*mDim + d]);
        }

        // reassign monomer to cell
        mMonomers[i].mCellID = whichCell(mMonomers[i]);
        mCells[mMonomers[i].mCellID].addMonomer(i);
    }

    *mOut << "       ->  rebuild monomer indexes information ... " << std::endl;
    // build up monomer indexes information
    buildMonomerIndexes(); 


    *mOut << "       ->  check energy ... " << std::endl;
    double energy = potential();
    if(abs(energy - mEnergy) > DOUBLE_TOLERANCE){
        errorMsg("LipidModel::unpack", "energy are incorrect !");
        exit(1);
    }
    *mOut << "           * recalculated energy: " << std::setprecision(20) << energy << "\n"
          << "           * energy from package: " << std::setprecision(20) << mEnergy << "\n";

    *mOut << "       -> check model integrity ... " << std::endl;
    // check model integrity
    checkModelIntegrity();

    // do measurement of observables
    doMeasure();
    *mOut << "     => Finish unpacking and rebuilding sucessfully !\n" << std::endl;    
}

void LipidModel::resetStatistics(){
    for(int i=0; i<NUMBEROFSTATS; ++i)
        stat[i] = 0;
}

RandomAccessNeighborList LipidModel::findWaterNeighborList(const long& id){
    RandomAccessNeighborList nlist;
    MonomerIndex iptr, jptr, kptr;
    CellIndex c_kptr;
    long ci, ck;
    iptr = mMonomerIndexes[id];
    ci = mMonomers[id].mCellID;
    jptr = mCells[ci].mMonomers.begin();
    while(jptr != mCells[ci].mMonomers.end()){
        if(jptr != iptr && mMonomers[*jptr].mType == WATER && canBeBonded(mMonomers[*iptr], mMonomers[*jptr])){
            nlist.addNeighbor(*jptr);
        }
        jptr++;
    }

    c_kptr = mCells[ci].mNeighbors.begin();
    while( c_kptr != mCells[ci].mNeighbors.end() ){
        ck = *c_kptr;
        kptr = mCells[ck].mMonomers.begin();
        while(kptr != mCells[ck].mMonomers.end()){
            if(mMonomers[*kptr].mType == WATER && canBeBonded(mMonomers[*iptr], mMonomers[*kptr])){
                nlist.addNeighbor(*kptr);
            }
            kptr++;
        }
        c_kptr++;
    }
    return nlist;
}

bool LipidModel::isCellClose(const long& monomer_id, const long& t_cell_id, const double& r){

/*
 
                  -------- 
                /|        /| 
               / |       / |
               ---------   |
              |   ------|- /
              | /       | /
              |/        |/  
           (*) --------- 

*/
    long me_cell_id = mMonomers[monomer_id].mCellID;

    if(mNumberOfCellsPerLine <=2 || me_cell_id == t_cell_id) return true;

    const long* me_coord = mCellID2Coord[me_cell_id];
    const long* t_coord = mCellID2Coord[t_cell_id];

    double dis(0);
    long* k = new long[mDim];
    long* a = new long[mDim];
    long* d = new long[mDim];


    for(int i=0; i<mDim; ++i)   {
        d[i] = me_coord[i] - t_coord[i];
        k[i] = 0;
        a[i] = 0;
    }

    for(int i=0; i<mDim; ++i){
        if(d[i] != 0){
            k[i] = 1;
            if(d[i] == -1 || d[i] > 1){
                a[i] = 1;
            }else{
                a[i] = 0;
            }
        }
    }

    for(int i=0; i<mDim; ++i){
        dis += k[i]*( pow( mCellCornerCoord[me_cell_id]->get(i) + a[i]*mLengthOfCell - mMonomers[monomer_id].mCoord->get(i), 2.0 ) );
    }

    delete [] k;
    delete [] a;
    delete [] d;

    return (sqrt(dis) < r);
}
