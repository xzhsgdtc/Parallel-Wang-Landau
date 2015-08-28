
#include "HistogramND.h"
#include <cstdlib>
#include <cmath>
#include <iomanip>


/*****************************************
 * :::::::::::::::::::::::::::::::::::::
 * 
 * Implementation of Bin class
 * 
 */

HistogramND::Bin::Bin():count(0),dos(0),dim(0){  
    min = NULL;
    max = NULL;
}

HistogramND::Bin::~Bin(){
    if(min != NULL) delete [] min;
    if(max != NULL) delete [] max;
}

HistogramND::Bin::Bin(const Bin& one){
    dim = one.dim;
    count = one.count;
    dos = one.dos;
    min = new double[dim];
    max = new double[dim];
    for(long i=0; i<dim; ++i){
        min[i] = one.min[i];
        max[i] = one.max[i];
    }
}

void HistogramND::Bin::init(long d){
    dim = d;
    min = new double[dim];
    max = new double[dim];
    count = 0;
    dos = 0;
}

void HistogramND::Bin::set(long d, double s, double e){ 
    min[d] = s; 
    max[d] = e; 
}


HistogramND::Bin& HistogramND::Bin::operator=(Bin& one){
    if(this == &one) return *this;
    dim = one.dim;
    count = one.count;
    dos = one.dos;
    if(min != NULL) delete [] min;
    if(max != NULL) delete [] max;
    min = new double[dim];
    max = new double[dim];
    for(long i=0; i<dim; ++i){
        min[i] = one.min[i];
        max[i] = one.max[i];
    }
    return *this;
}


/****************************************************
 * :::::::::::::::::::::::::::::::::::::::::::::::::
 * 
 *      Implementation of HistogramND class
 * 
 */

HistogramND::HistogramND(long d):mDim(d), mNumberOfBins(0), mFlatnessType(DEFAULT_FLAT), mUpdateType(DEFAULT_UPDATE), mBins(NULL) {
    mStepEachDim   = new double[mDim];
    mMinEachDim    = new double[mDim];
    mMaxEachDim    = new double[mDim];
    mBinsEachDim   = new long[mDim];
    mWeightEachDim = new long[mDim];
    mDiscreteFlag = new bool[mDim];
    mTempIndex  = new long[mDim];
    setOutput(std::cout);
}

HistogramND::HistogramND(const HistogramND& one) {
    mDim = one.mDim;
    mFlatnessType = one.mFlatnessType;
    mUpdateType = one.mUpdateType;
    mBinsEachDim = new long[mDim];
    mStepEachDim = new double[mDim];
    mMinEachDim = new double[mDim];
    mMaxEachDim = new double[mDim];
    mWeightEachDim = new long[mDim];
    mDiscreteFlag = new bool[mDim];
    mTempIndex  = new long[mDim];
    
    for(long i=0; i<mDim; ++i){
        mWeightEachDim[i] = one.mWeightEachDim[i];
        mDiscreteFlag[i] = one.mDiscreteFlag[i];
        mBinsEachDim[i]   = one.mBinsEachDim[i];
        mStepEachDim[i]   = one.mStepEachDim[i];
        mMinEachDim[i]    = one.mMinEachDim[i];
        mMaxEachDim[i]    = one.mMaxEachDim[i];
    }
    
    mNumberOfBins = one.mNumberOfBins;
    
    mBins = new Bin[mNumberOfBins];
    for(long i=0; i<mNumberOfBins; ++i){
        mBins[i] = one.mBins[i];
    }
}

HistogramND& HistogramND::operator =(const HistogramND& one){
    if(this == &one) return *this;
    
    if( mBinsEachDim   != NULL )    delete [] mBinsEachDim;
    if( mStepEachDim   != NULL )    delete [] mStepEachDim;
    if( mMinEachDim    != NULL )    delete [] mMinEachDim;
    if( mMaxEachDim    != NULL )    delete [] mMaxEachDim;
    if( mWeightEachDim != NULL )    delete [] mWeightEachDim;
    if( mDiscreteFlag != NULL )    delete [] mDiscreteFlag;
    if( mTempIndex != NULL )    delete [] mTempIndex;
    
    mDim = one.mDim;
    mFlatnessType = one.mFlatnessType;
    mUpdateType = one.mUpdateType;
    mBinsEachDim = new long[mDim];
    mStepEachDim = new double[mDim];
    mMinEachDim = new double[mDim];
    mMaxEachDim = new double[mDim];
    mWeightEachDim = new long[mDim];
    mDiscreteFlag = new bool[mDim];
    mTempIndex  = new long[mDim];
    
    for(long i=0; i<mDim; ++i){
        mBinsEachDim[i] = one.mBinsEachDim[i];
        mStepEachDim[i] = one.mStepEachDim[i];
        mMinEachDim[i] = one.mMinEachDim[i];
        mMaxEachDim[i] = one.mMaxEachDim[i];
        mWeightEachDim[i] = one.mWeightEachDim[i];
        mDiscreteFlag[i] = one.mDiscreteFlag[i];
    }
    
    mNumberOfBins = one.mNumberOfBins;
    
    if(mBins != NULL) delete [] mBins;
    mBins = new Bin[mNumberOfBins];
    for(long i=0; i<mNumberOfBins; ++i){
        mBins[i] = one.mBins[i];
    }
    return *this;
}

HistogramND::~HistogramND() {
    if(mBins != NULL) delete [] mBins;
    if(mBinsEachDim  != NULL)    delete [] mBinsEachDim;
    if(mStepEachDim  != NULL)    delete [] mStepEachDim;
    if(mMinEachDim   != NULL)    delete [] mMinEachDim;
    if(mMaxEachDim   != NULL)    delete [] mMaxEachDim;
    if(mWeightEachDim!= NULL)    delete [] mWeightEachDim;
    if(mDiscreteFlag!= NULL)    delete [] mDiscreteFlag;
    if( mTempIndex != NULL )    delete [] mTempIndex;
}


void HistogramND::set(long dim, double start, double end, long num, bool flag){
    if(dim > mDim){
        (*mOut) << "Error: set, out of mDim of histogram !" << std::endl;
        exit(1);
    }
    mBinsEachDim[dim] = num;
    mMinEachDim[dim] = start;
    mMaxEachDim[dim] = end;
    mDiscreteFlag[dim] = flag;
    if(mDiscreteFlag[dim]){
        mStepEachDim[dim] = floor( (end-start)/(double)num+0.5 );
    }else{
        mStepEachDim[dim] = (end-start)/(double)num;
    }
}

void HistogramND::build(){
    *mOut << "     => start building histogram object ! \n";
    mNumberOfBins = 1;
    for(long i=0; i<mDim; ++i){
        mNumberOfBins *= mBinsEachDim[i]; 
        if(i == 0) 
            mWeightEachDim[0] = 1;
        else
            mWeightEachDim[i] = mBinsEachDim[i-1];
    }
    
    if(mBins != NULL) delete [] mBins;
    mBins = new Bin[mNumberOfBins];    
    
    double* start = new double[mDim];
    double* end = new double[mDim];
    long* dim_count = new long[mDim];
    
    long temp_index;
    for(long index = 0; index < mNumberOfBins; ++index){
        mBins[index].init(mDim);
        temp_index = index;
        for(long d=0; d<mDim; ++d){
            dim_count[d] = temp_index%mBinsEachDim[d];
            temp_index -= dim_count[d];
            temp_index /= mBinsEachDim[d];
            start[d] = mMinEachDim[d] + mStepEachDim[d]*dim_count[d];
            end[d] = mMinEachDim[d]+mStepEachDim[d]*(dim_count[d]+1);
            if(end[d] > mMaxEachDim[d] || ((mMaxEachDim[d] - end[d]) < mStepEachDim[d] ) && dim_count[d] == mBinsEachDim[d] ) end[d] = mMaxEachDim[d];
            mBins[index].set(d,start[d],end[d]);
        }
    }
    *mOut << "        -> total number of bins: " << mNumberOfBins << "\n"
          << "        -> histogram dimension: " << mDim << "\n"
          << "        -> bins each dim: " ;
    for(int d=0; d<mDim; ++d) *mOut << mBinsEachDim[d] << "  ";
    *mOut <<"\n        -> weights each dim: " ;
    for(int d=0; d<mDim; ++d) *mOut << mWeightEachDim[d] << "  ";
    *mOut <<"\n        -> min. each dim: " ;
    for(int d=0; d<mDim; ++d) *mOut << mMinEachDim[d] << "  ";
    *mOut <<"\n        -> max. each dim: " ;
    for(int d=0; d<mDim; ++d) *mOut << mMaxEachDim[d] << "  ";
    *mOut <<"\n        -> step each dim: " ;
    for(int d=0; d<mDim; ++d) *mOut << mStepEachDim[d] << "  ";
 
    *mOut << "\n     => build Histogam sucessfully ! " << std::endl;
    
    delete [] start;
    delete [] end;
}

bool HistogramND::flat(const double& fcriterion, const double& modfactor){
    bool result = false;
    switch(mFlatnessType){
        case DEFAULT_FLAT:
        {
            double ave(0);
            long min = mBins[0].count;
            for(long i=0; i<mNumberOfBins; ++i){
                ave += (double)mBins[i].count;
                if(mBins[i].count < min) min = (double)mBins[i].count;
            }
            ave /= (double)mNumberOfBins;
            double standard = ave*fcriterion;
            if(min > standard) result = true;
        }
        break;
        case ZHOU_PRE:
        {
            result = true;
            long min_hit = long(1/sqrt( modfactor ));
            for(long i=0; i<mNumberOfBins; ++i){
                if(mBins[i].count < min_hit){
                    result = false;
                    break;
                }
            }
        }
        break;
    }
    return result;
}

void HistogramND::reset(){
    for(long i=0; i<mNumberOfBins; ++i){
        mBins[i].count =0;
    }
}

long HistogramND::update(const double* observables, double mf){
    long updated_bin_index(-1);
    switch(mUpdateType){
        case DEFAULT_UPDATE:

            long i = bIndex(observables);
            switch(i){
                case -1:
                    (*mOut) << "Error: update, at least one mDim is illegel!" << std::endl;
                    (*mOut) << "value = " << observables[i] << std::endl;
                    exit(1);
                    break;
                case -2:
                    (*mOut) << "Error: HistogramND->update, index out of range for at least one dimension!" << std::endl;
                    exit(1);
                    break;
                case -3:
                    (*mOut) << "Error: HistogramND->update, index out of range as a whole!" << std::endl;
                    exit(1);
                    break;
            }

            mBins[i].count++;
            mBins[i].dos += mf;
            updated_bin_index =  i;
        break;
    }
    return updated_bin_index;
}


long HistogramND::bIndex(const double* observables) const {
    // negative return means:
    // -1: observable out of range
    // -2: index out of range for certrain dimension
    // -3: index out of range as a whole
    for(int i=0; i<mDim; ++i){
        if(observables[i] < mMinEachDim[i] || observables[i] > mMaxEachDim[i])
            return -1;
    }
    long i(0);
    for(long d=0; d<mDim; ++d){
        mTempIndex[d] = findIndex(d, observables[d]);
        if(mTempIndex[d] >= mBinsEachDim[d] || mTempIndex[d] < 0) return -2;
        i += mTempIndex[d]*mWeightEachDim[d];
    }
    if(i<0 || i >= mNumberOfBins ) return -3;
    return i;
}

double HistogramND::dosAt(const double* observables){
    long i = bIndex(observables);
    if (i < 0) return i;
    else return mBins[i].dos;
}


long HistogramND::findIndex(long d, double quantity) const{
    if(mDiscreteFlag[d]){
      return floor( (quantity+0.5-mMinEachDim[d]) / mStepEachDim[d]);  
    } else{
        if (quantity > (mMaxEachDim[d] - mStepEachDim[d])  && quantity <= mMaxEachDim[d]) return mBinsEachDim[d]-1;
        return floor((quantity-mMinEachDim[d])/mStepEachDim[d]);
    }
}


void HistogramND::print(std::ostream& out){

    out << "### HISTOGRAMND INFORMATION\n"
              << "\n# UpdateType: " << mUpdateType << "\n"
              << "\n# FlatnessType: " << mFlatnessType << "\n"
              << "\n# Dimension: " << mDim << "\n"
              << "\n# NumberofBins: " << mNumberOfBins << "\n";
    out << "\n# NumberOfBinsEachDimension: ";
    for(long k=0; k<mDim; ++k){
        out << std::setw(8) << mBinsEachDim[k] << " ";
    }
    out << "\n";

    out << "\n# WeightOfEachDimension: ";
    for(long k=0; k<mDim; ++k){
        out << std::setw(8) << mWeightEachDim[k] << " ";
    }
    out << "\n";

    out << "\n# StepOfEachDimension: ";
    for(long k=0; k<mDim; ++k){
        out << std::setw(8) << mStepEachDim[k] << " ";
    }
    out << "\n";

    out << "\n# MaximumOfEachDimension: ";
    for(long k=0; k<mDim; ++k){
        out << std::setw(8) << mMaxEachDim[k] << " ";
    }
    out << "\n";

    out << "\n# MinimumOfEachDimension: ";
    for(long k=0; k<mDim; ++k){
        out << std::setw(8) << mMinEachDim[k] << " ";
    }
    out << "\n" << std::endl;

    long cut = 1;
    for(long i=0; i<mDim-1; ++i){
        cut *= mBinsEachDim[i];
    }
    for(long i=0; i<mNumberOfBins; ++i){
        if(mDim >1 && i!=0 && i % cut == 0) out << std::endl;
        out << std::setw(10) << i << "  " ;
        for (long d = 0; d < mDim; ++d) {
            if(mDiscreteFlag[d]){
                out << std::setw(10) << ceil(mBins[i].min[d]-0.5) << "  "
                    << std::setw(10) << ceil(mBins[i].min[d]-0.5) << "  ";
            }else{
                out << std::fixed << std::setw(10) << mBins[i].min[d] << "  "
                    << std::setw(10) << mBins[i].max[d] << "  ";
            }
        }
        out  << std::setw(10) << mBins[i].count << "  "
             << std::scientific   << std::setw(20) << mBins[i].dos << std::endl;
    }
}


void HistogramND::read(std::string fname){
    std::ifstream fin(fname.c_str());
    if(!fin.is_open()){
        *mOut << "Error: cannot open file: " << fname << std::endl;
        exit(1);
    }
    read(fin);
}

void HistogramND::read(std::ifstream& fin){
    std::string line;
    long count(0);
    double* entry(NULL);
    long nofe(0);
    while(getline(fin, line)){
        std::istringstream ins;
        ins.str(line);
        std::string word;
        if (line[0] == '#'){
        // read in parameters
            ins >> word;
            ins >> word;
            if(word.compare("UpdateType:") == 0){
                int t;
                ins >> t; 
                mUpdateType = (HistogramUpdateType)t;
            }else if(word.compare("FlatnessType:") == 0){
                int t;
                ins >> t;
                mFlatnessType = (FlatnessType)t;
            }else if(word.compare("Dimension:") == 0){
                ins >> mDim;
                nofe = 3+2*mDim;
                entry = new double[nofe];
                if( mBinsEachDim   != NULL )    delete [] mBinsEachDim;
                if( mStepEachDim   != NULL )    delete [] mStepEachDim;
                if( mMinEachDim    != NULL )    delete [] mMinEachDim;
                if( mMaxEachDim    != NULL )    delete [] mMaxEachDim;
                if( mWeightEachDim != NULL )    delete [] mWeightEachDim;
                if( mDiscreteFlag  != NULL )    delete [] mDiscreteFlag;
                mStepEachDim   = new double[mDim];
                mMinEachDim    = new double[mDim];
                mMaxEachDim    = new double[mDim];
                mBinsEachDim   = new long[mDim];
                mWeightEachDim = new long[mDim];
                mDiscreteFlag  = new bool[mDim];
            }else if(word.compare("NumberofBins:") == 0){
                ins >> mNumberOfBins;
            }else if(word.compare("NumberOfBinsEachDimension:") == 0){
                for(int i=0; i<mDim; ++i) ins >> mBinsEachDim[i];
            }else if(word.compare("WeightOfEachDimension:") == 0){
                for(int i=0; i<mDim; ++i) ins >> mWeightEachDim[i];
            }else if(word.compare("StepOfEachDimension:") == 0){
                for(int i=0; i<mDim; ++i) ins >> mStepEachDim[i];
            }else if(word.compare("MaximumOfEachDimension:") == 0){
                for(int i=0; i<mDim; ++i) ins >> mMaxEachDim[i];
            }else if(word.compare("MinimumOfEachDimension:") == 0){
                for(int i=0; i<mDim; ++i) ins >> mMinEachDim[i];
                build();
            }
        }else{
            if(line.size() != 0){
                for(int i=0; i< nofe; ++i) ins >> entry[i];
                if (count != long(entry[0])){
                   *mOut << "Error: error in reading histogram file!!!" << std::endl;
                   exit(1);
                }
                for(int d=0; d<mDim; ++d){
                    if( abs(mBins[count].min[d] - entry[1+d*2]) > HISTOGRAM_DOUBLE_TOLERANCE 
                            || abs(mBins[count].max[d] - entry[2+d*2]) > HISTOGRAM_DOUBLE_TOLERANCE ){
                        *mOut << "Error: error in reading histogram file!!!" << std::endl;
                        exit(1);
                    }
                }                  
                mBins[count].count = entry[nofe-2];
                mBins[count].dos = entry[nofe-1];
                count++;
            }      
        }
    }
    if (entry != NULL) delete [] entry;
}

void HistogramND::dosNorm(){
    double min = -1;
    for(long i=0; i<mNumberOfBins; ++i){
        if(min == -1){
            min = mBins[i].dos;
        }else if(mBins[i].dos < min){
            min = mBins[i].dos;
        }
    }
    for(long i=0; i<mNumberOfBins; ++i){
        mBins[i].dos -= min;
    }
}

void HistogramND::dosCopyAll(double* logge){
    for(int i=0; i<mNumberOfBins; ++i){
        logge[i] = mBins[i].dos;
    }
}

void HistogramND::dosAssignAll(double* logge){
    for(int i=0; i<mNumberOfBins; ++i){
        mBins[i].dos = logge[i];
    }
}

void HistogramND::backup(std::ofstream& fout){
    
    int sod = sizeof(double);
    int sol = sizeof(long);
    fout.write(reinterpret_cast<char*> (&mFlatnessType), sizeof(FlatnessType));
    fout.write(reinterpret_cast<char*> (&mUpdateType), sizeof(HistogramUpdateType));
    fout.write(reinterpret_cast<char*> (&mDim)      ,           sol);
    fout.write(reinterpret_cast<char*> (mBinsEachDim)      ,    mDim * sol);
    fout.write(reinterpret_cast<char*> (mStepEachDim)      ,    mDim * sod);
    fout.write(reinterpret_cast<char*> (mMinEachDim)       ,    mDim * sod);
    fout.write(reinterpret_cast<char*> (mMaxEachDim)       ,    mDim * sod);
    fout.write(reinterpret_cast<char*> (mDiscreteFlag)     ,  mDim * sizeof (bool));
    
    //long s =mDim*sod;
    for(long i=0; i<mNumberOfBins; ++i){
        fout.write(reinterpret_cast<char*>(&mBins[i].dos),sod);
        fout.write(reinterpret_cast<char*>(&mBins[i].count),sol);
    }
    

}

void HistogramND::recover(std::ifstream& fin) {
    // file exists, recovery from it
    int sod = sizeof (double);
    int sol = sizeof (long);
    fin.read(reinterpret_cast<char*> (&mFlatnessType), sizeof(FlatnessType));
    fin.read(reinterpret_cast<char*> (&mUpdateType), sizeof(HistogramUpdateType));
    fin.read(reinterpret_cast<char*> (&mDim), sol);
    fin.read(reinterpret_cast<char*> (mBinsEachDim), mDim * sol);
    fin.read(reinterpret_cast<char*> (mStepEachDim), mDim * sod);
    fin.read(reinterpret_cast<char*> (mMinEachDim), mDim * sod);
    fin.read(reinterpret_cast<char*> (mMaxEachDim), mDim * sod);
    fin.read(reinterpret_cast<char*> (mDiscreteFlag), mDim * sizeof (bool));

    build();
    //long s = mDim*sod;
    for (long i = 0; i < mNumberOfBins; ++i) {
        fin.read(reinterpret_cast<char*>(&mBins[i].dos), sod);
        fin.read(reinterpret_cast<char*>(&mBins[i].count), sol);
    }

}

void HistogramND::setOutput(std::ostream& out){
    mOut = &out;
}



