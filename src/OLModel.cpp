/* 
 * File:   OLModel.cpp
 * Author: jerry
 * 
 * Created on August 16, 2012, 5:09 PM
 */

#include "OLModel.h"
#include "ToString.h"
#include "InputOutput.h"

OLModel::OLModel() {
    mCells = NULL;
    mMonomerIndexes = NULL;
    observables = NULL;
    backup_observables = NULL;
    mMonomers = NULL;
    _dim = 3;
    mLengthOfSimBox = 0;
    numberOfmMonomers = 0;
    initType = true;
    delta_coord = NULL;
    boundaryType = 1;
    moveProposal = -1;
    setOutput(cout);
}

OLModel::OLModel(const OLModel& one) {
    _energy = one._energy;
    _prevEnergy = one._prevEnergy;
    moveProposal = one.moveProposal;
    mLengthOfSimBox = one.mLengthOfSimBox;
    mincellsize = one.mincellsize;
    mLengthOfCell = one.mLengthOfCell;
    mMaxMoveDistance = one.mMaxMoveDistance;
    mNumberOfCells = one.mNumberOfCells;
    mNumberOfCellsPerLine = one.mNumberOfCellsPerLine;
    mNumberOfCellsPerLine2 = one.mNumberOfCellsPerLine2;
    density = one.density;
    initType = one.initType;
    boundaryType = one.boundaryType;
    _dim = one._dim;
    numberOfmMonomers = one.numberOfmMonomers;
    mBackupMonomer = one.mBackupMonomer;
    delta_coord = new double[_dim];
    if(one.observables != NULL){
        observables = new double[numberOfObservables];
        for(long i=0; i< numberOfObservables; ++i)
            observables[i] = one.observables[i];
    }
    if(one.backup_observables != NULL){
        backup_observables = new double[numberOfObservables];
        for(long i=0; i< numberOfObservables; ++i)
            backup_observables[i] = one.backup_observables[i];
    }
    if (one.mCells != NULL && one.mMonomers != NULL) {
            mCells = new Box[mNumberOfCells];
            mMonomers = new Monomer[numberOfmMonomers];
            if (_dim == 2 || _dim == 3) {
                for (long i = 0; i < mNumberOfCells; ++i)
                    mCells[i] = (one.mCells[i]);
                for (long j = 0; j < numberOfmMonomers; ++j)
                    mMonomers[j] = one.mMonomers[j];
            } else {
                (*myout) << "OPERATOR= : Error: _dim is incorrect !" << std::endl;
                exit(1);
            }
    }else {
        return;
    }
    buildmMonomerIndexes();
}

OLModel::~OLModel() {
    if (mCells != NULL)
        delete [] mCells;
    if(mMonomers != NULL)
        delete [] mMonomers;
    if (mMonomerIndexes != NULL)
        delete [] mMonomerIndexes;
    if(observables != NULL)
        delete [] observables;
    if(backup_observables != NULL)
        delete [] backup_observables;
    if(delta_coord != NULL)
        delete [] delta_coord;
}

OLModel& OLModel::operator =(OLModel& one) {
    if (this == &one) return *this;
    _energy = one._energy;
    _prevEnergy = one._prevEnergy;
    mLengthOfSimBox = one.mLengthOfSimBox;
    mincellsize = one.mincellsize;
    mLengthOfCell = one.mLengthOfCell;
    moveProposal = one.moveProposal;
    mMaxMoveDistance = one.mMaxMoveDistance;
    mNumberOfCells = one.mNumberOfCells;
    mNumberOfCellsPerLine = one.mNumberOfCellsPerLine;
    mNumberOfCellsPerLine2 = one.mNumberOfCellsPerLine2;
    density = one.density;
    initType = one.initType;
    boundaryType = one.boundaryType;
    _dim = one._dim;
    numberOfmMonomers = one.numberOfmMonomers;
    mBackupMonomer = one.mBackupMonomer;
    if(delta_coord != NULL) delete [] delta_coord;
    delta_coord = new double[_dim];
    if(observables != NULL) delete [] observables;
    if(one.observables != NULL){
        observables = new double[numberOfObservables];
        for(long i=0; i< numberOfObservables; ++i)
            observables[i] = one.observables[i];
    }
    
    if(backup_observables != NULL) delete [] backup_observables;
    if(one.backup_observables != NULL){
        backup_observables = new double[numberOfObservables];
        for(long i=0; i< numberOfObservables; ++i)
            backup_observables[i] = one.backup_observables[i];
    }
    
    if (one.mCells == NULL || one.mMonomers == NULL) return *this;
    if (mCells != NULL) delete [] mCells;
    if (mMonomers != NULL) delete [] mMonomers;
    mCells = new Box[mNumberOfCells];
    mMonomers = new Monomer[numberOfmMonomers];
    if (_dim == 2 || _dim == 3) {
        for (long i = 0; i < mNumberOfCells; ++i)
            mCells[i] = (one.mCells[i]);
        for(long j=0; j<numberOfmMonomers; ++j)
            mMonomers[j] = one.mMonomers[j];
    } else {
        (*myout) << "OPERATOR= : Error: _dim is incorrect !" << std::endl;
        exit(1);
    }
    buildmMonomerIndexes();
    return *this;
}



void OLModel::printIndex(ostream& out){
    out << setw(10)<< "Index" << '\t' << setw(10)<< "mMonomersIndex" << '\t'
        << setw(10) << "IndexCheck" << '\t'
        << setw(10)<< "Coord-x" << '\t' << setw(10)<< "Coord-y" << '\t' << setw(10)<< "Coord-z" 
        << '\t' << setw(10)<< "Type" << '\t' << setw(10)<< "Box"  << '\t'<< setw(10) << "Box-Check"<< endl;
    for (long i = 0; i < numberOfmMonomers; ++i) {
        long index = *mMonomerIndexes[i];
        long box = whichBox(mMonomers[i]);
        out << setw(10) << i << '\t'
                << setw(10) << index << '\t'
                << setw(10) << (index == i) << '\t'
                << setw(10) << mMonomers[index].x->get(0) << '\t'
                << setw(10) << mMonomers[index].x->get(1) << '\t'
                << setw(10) << mMonomers[index].x->get(2) << '\t'
                << setw(10) << mMonomers[index].type << '\t'
                << setw(10) << mMonomers[index].curr_box << '\t'
                << setw(10) << box  << '\t' 
                << setw(10 )<< (box == mMonomers[index].curr_box)<< endl;
    }
}


void OLModel::printCoord(ostream& out){
}

void OLModel::printStruct(ostream& out) {
}


void OLModel::backup(ofstream& fout){
    fout.write(reinterpret_cast<char*>(&initType), sizeof(bool));
    fout.write(reinterpret_cast<char*>(&moveProposal), sol);
    fout.write(reinterpret_cast<char*>(&numberOfmMonomers),sol );
    fout.write(reinterpret_cast<char*>(&_dim), sol);
    fout.write(reinterpret_cast<char*>(&numberOfObservables), sol);
    fout.write(reinterpret_cast<char*>(&boundaryType), sol);
    fout.write(reinterpret_cast<char*>(observables), numberOfObservables*sod);
    fout.write(reinterpret_cast<char*>(backup_observables), numberOfObservables*sod);
    fout.write(reinterpret_cast<char*>(&_energy), sod);
    fout.write(reinterpret_cast<char*>(&_prevEnergy),sod );
    fout.write(reinterpret_cast<char*>(&mLengthOfSimBox),sod );
    fout.write(reinterpret_cast<char*>(&density), sod);
    fout.write(reinterpret_cast<char*>(&mincellsize), sod);
    fout.write(reinterpret_cast<char*>(&mMaxMoveDistance), sod);
    fout.write(reinterpret_cast<char*>(delta_coord), _dim*sod);
    
    mBackupMonomer.backup(fout);
    
    for(long i=0; i<mNumberOfCells; ++i){
        mCells[i].backup(fout);
    }
    
    for(long i=0; i<numberOfmMonomers; ++i){
        mMonomers[i].backup(fout);
    }
}

void OLModel::recover(ifstream& fin){
    long sod = sizeof(double);
    long sol = sizeof(long);
    
    fin.read(reinterpret_cast<char*>(&initType), sizeof(bool));
    fin.read(reinterpret_cast<char*>(&moveProposal), sol);
    fin.read(reinterpret_cast<char*>(&numberOfmMonomers),sol );
    fin.read(reinterpret_cast<char*>(&_dim), sol);

    delta_coord = new double[_dim];
    
    fin.read(reinterpret_cast<char*>(&numberOfObservables), sol);
    fin.read(reinterpret_cast<char*>(&boundaryType), sol);

    observables = new double[numberOfObservables];
    backup_observables = new double[numberOfObservables];
    
    fin.read(reinterpret_cast<char*>(observables), numberOfObservables*sod);
    fin.read(reinterpret_cast<char*>(backup_observables), numberOfObservables*sod);
    fin.read(reinterpret_cast<char*>(&_energy), sod);
    fin.read(reinterpret_cast<char*>(&_prevEnergy),sod );
    fin.read(reinterpret_cast<char*>(&mLengthOfSimBox),sod );
    fin.read(reinterpret_cast<char*>(&density), sod);
    fin.read(reinterpret_cast<char*>(&mincellsize), sod);
    fin.read(reinterpret_cast<char*>(&mMaxMoveDistance), sod);
    fin.read(reinterpret_cast<char*>(delta_coord), _dim*sod);

    
    mBackupMonomer.recover(fin);

    buildmCells();
    
    if(mMonomers != NULL) delete [] mMonomers;
    mMonomers = new Monomer[numberOfmMonomers];
    
    for(long i=0; i<mNumberOfCells; ++i){
        mCells[i].recover(fin);
    }
    
    for(long i=0; i<numberOfmMonomers; ++i){
        mMonomers[i].recover(fin);
    }
    
    buildmMonomerIndexes();
}

void OLModel::setOutput(ostream& out){
    myout = &out;
}

void OLModel::error(ErrorType et, string func_name, string msg, ostream& error_out){
    switch(et){
        case ILLEGAL:
            error_out << "Error occurs inside FUNCTION: " << func_name << ", " << msg << " is illegal!" << endl;
            break;
        case OUTOFRANGE:
            error_out << "Error occurs inside FUNCTION: " << func_name << ", " << msg << " is out of range!" << endl;
            break;
        case DEFAULT:
            error_out << "Error occurs inside FUNCTION: " << func_name << ", " << msg << endl;
            break;
    }
}
