#include "Monomer.h"

Monomer::Monomer() {
    mCoord = NULL;
    mCellID = 0;
    mID = -1;
    mDim = -1;
    mLipidWaterID = 0;
}

Monomer::Monomer(double _x, double _y, double _z) {
    mDim = 3;
    mCoord = new MathVector<double>(_x, _y, _z);
    mCellID = -1;
    mID = -1;
    mLipidWaterID = 0;
}

Monomer::Monomer(double _x, double _y) {
    mDim = 2;
    mCoord = new MathVector<double>(_x, _y);
    mCellID = -1;
    mID = -1;
    mLipidWaterID = 0;
}

Monomer::Monomer(const Monomer& a) {
    mCoord = new MathVector<double>(a.mDim);
    *mCoord = *a.mCoord;
    mDim = a.mDim;
    mType = a.mType;
    mID = a.mID;
    mCellID = a.mCellID;
    
    mLipidWaterID = a.mLipidWaterID;
}

Monomer& Monomer::operator=(const Monomer& a) {
    if (&a == this) return *this;
    
    if (mCoord != NULL) delete mCoord;
    mCoord = new MathVector<double>(a.mDim);
    *mCoord = *a.mCoord;
    mDim = a.mDim;
    mType = a.mType;
    mID = a.mID;
    mCellID = a.mCellID;
   
    mLipidWaterID = a.mLipidWaterID;
    return *this;
}

Monomer::~Monomer() {
    if (mCoord != NULL) delete mCoord;
}

void Monomer::backup(std::ofstream& fout) {

    long sol = sizeof (long);

    fout.write(reinterpret_cast<char*> (&mType), sizeof (MonomerType));
    fout.write(reinterpret_cast<char*> (&mDim), sol);
    fout.write(reinterpret_cast<char*> (&mID), sol);
    fout.write(reinterpret_cast<char*> (&mCellID), sol);
    fout.write(reinterpret_cast<char*> (&mLipidWaterID), sol);
  
    mCoord->backup(fout);
}

void Monomer::recover(std::ifstream& fin) {
    long sol = sizeof (long);

    fin.read(reinterpret_cast<char*> (&mType), sizeof (MonomerType));
    fin.read(reinterpret_cast<char*> (&mDim), sol);
    fin.read(reinterpret_cast<char*> (&mID), sol);
    fin.read(reinterpret_cast<char*> (&mCellID), sol);
    fin.read(reinterpret_cast<char*> (&mLipidWaterID), sol);
 
    if(mCoord != NULL) delete mCoord;
    mCoord = new MathVector<double>(mDim);
    mCoord->recover(fin);
}

void Monomer::print(std::ostream& out) {
    out << "Type  " << mType << "\n"
            << "Dim  " << mDim << "\n"
            << "ID  " << mID << "\n"
            << "CellID  " << mCellID << "\n"
            << "LipidWaterID  " << mLipidWaterID << "\n"

            << std::endl;
    if(mCoord != NULL) mCoord->print(out);
    out << std::endl;
}

//==================================================//
//       implementation of Cell class below         //
//==================================================// 

Cell::Cell(const Cell& orig) {
    mMonomers = orig.mMonomers;
    mNeighbors = orig.mNeighbors;
    mNumberOfMonomers = orig.mNumberOfMonomers;
    mNumberOfNeighbors = orig.mNumberOfNeighbors;
}

Cell& Cell::operator=(const Cell& orig) {
    if (this == &orig) return *this;
    mMonomers = orig.mMonomers;
    mNeighbors = orig.mNeighbors;
    mNumberOfMonomers = orig.mNumberOfMonomers;
    mNumberOfNeighbors = orig.mNumberOfNeighbors;
    return *this;
}

Cell::~Cell() {
    mMonomers.clear();
    mNeighbors.clear();
}

void Cell::backup(std::ofstream& fout) {
    long sol = sizeof (long);
    long count = 0;
    fout.write(reinterpret_cast<char*> (&mNumberOfMonomers), sol);
    fout.write(reinterpret_cast<char*> (&mNumberOfNeighbors), sol);

    // write mMonomers
    long* tmMonomers = new long[mNumberOfMonomers];
    MonomerIndex m_iter = mMonomers.begin();
    while (m_iter != mMonomers.end()) {
        tmMonomers[count] = *m_iter;
        count++;
        m_iter++;
    }
    if(count != mNumberOfMonomers){
        std::cout << "Error, backup(...) @ Cell, mNumberOfMonomers != count" << std::endl;
        delete [] tmMonomers;
        exit(1);
    }
    fout.write(reinterpret_cast<char*> (tmMonomers), mNumberOfMonomers * sol);
    delete [] tmMonomers;
    
    // write mNeighbors
    count = 0;
    long* tmNeighbors = new long[mNumberOfNeighbors];
    CellIndex c_iter = mNeighbors.begin();
    while (c_iter != mNeighbors.end()) {
        tmNeighbors[count] = *c_iter;
        count++;
        c_iter++;
    }
    if(count != mNumberOfNeighbors){
        std::cout << "Error, backup(...) @ Cell, mNumberOfNeighbors != count" << std::endl;
        delete [] tmNeighbors;
        exit(1);
    }
    fout.write(reinterpret_cast<char*> (tmNeighbors), mNumberOfNeighbors * sol);
    delete [] tmNeighbors;
}

void Cell::recover(std::ifstream& fin) {
    long sol = sizeof (long);
    fin.read(reinterpret_cast<char*> (&mNumberOfMonomers), sol);
    fin.read(reinterpret_cast<char*> (&mNumberOfNeighbors), sol);
    mMonomers.clear();
    mNeighbors.clear();

    // read mMonomers
    long* tmMonomers = new long[mNumberOfMonomers];
    fin.read(reinterpret_cast<char*> (tmMonomers), mNumberOfMonomers * sol);
    for (long i = 0; i < mNumberOfMonomers; ++i) {
        mMonomers.push_back(tmMonomers[i]);
    }
    if(mMonomers.size() != mNumberOfMonomers){
        std::cout << "Error, recover(...) @ Cell, mNumberOfMonomers incorrect!" << std::endl;
        delete [] tmMonomers;
        exit(1);
    }
    delete [] tmMonomers;

    // read mNeighbors
    long* tmNeighbors = new long[mNumberOfNeighbors];
    fin.read(reinterpret_cast<char*> (tmNeighbors), mNumberOfNeighbors * sol);
    for (long i = 0; i < mNumberOfNeighbors; ++i) {
        mNeighbors.push_back(tmNeighbors[i]);
    }
    if(mNeighbors.size() != mNumberOfNeighbors){
        std::cout << "Error, recover(...) @ Cell, mNumberOfNeighbors incorrect!" << std::endl;
        delete [] tmNeighbors;
        exit(1);
    }
    delete [] tmNeighbors;
}

void Cell::print(std::ostream& out){
    out << "Monomers ( " << mNumberOfMonomers << " ):  ";  
    MonomerIndex m_iter = mMonomers.begin();
    while(m_iter != mMonomers.end()){
        out << *m_iter << "\t";
        m_iter++;
    }
    out << std::endl;
    out << "Neighbor Cells ( " << mNumberOfNeighbors << " ):  ";
    CellIndex c_iter = mNeighbors.begin();
    while(c_iter != mNeighbors.end()){
        out << *c_iter << "\t";
        c_iter++;
    }
    out << std::endl;
}

