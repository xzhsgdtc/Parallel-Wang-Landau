/* *********************************************************
 * :::::::::::::::::::::::::::::::::::::::::::::::::::::::::
 *
 *   This file contains the definition of classes Monomer
 *   and Cell.
 *
 *   @File:   Monomer.h/.cpp
 *   @Author: Guangjie Shi (Jerry)
 *
 *   @Version 1.0: Oct. 30, 2012
 *      Definition of classes Monomer and Cell
 *
 * :::::::::::::::::::::::::::::::::::::::::::::::::::::::::
 */


#ifndef MONOMER_H
#define MONOMER_H

#include "MathVector.h"
#include <iostream>
#include <fstream>
#include <list>


const long NUMBEROFMONOMERTYPES = 3;  // number of monomer types

/**
 * enum define the type of monomers involved in simulation
 */
enum MonomerType {
    WATER = 0,         // solution particle like water
    HYDROPHOBIC = 1,         // hydrophobic
    POLAR = 2          // polar 
};


/**
 * @class Monomer @file Monomer.h/.cpp
 *
 * Monomer class contain the coordinate, type, dimension, id, cellid etc. 
 * information of a monomer. 
 * 
 */

class Monomer {
public:

    MathVector<double>* mCoord; // coordinates
    MonomerType mType;  // type
    long mDim;          // dimension
    long mID;           // id
    long mCellID;       // cell id
    long mLipidWaterID;  // index in mLipidWater array

    /**
     * Constructor
     */
    Monomer();

    /**
     * constructor for 3D
     * @param x
     * @param y
     * @param z
     */
    Monomer(double x, double y, double z);

    /**
     * Constructor for 2D
     * @param x
     * @param y
     */
    Monomer(double x, double y);


    /**
     * Copy Constructor
     * @param a
     */
    Monomer(const Monomer& a);

    /**
     * Operator =
     * @param a
     * @return 
     */
    Monomer& operator=(const Monomer& a);


    /**
     * Destructor
     */
    ~Monomer();

    void print(std::ostream& out=std::cout);
    void backup(std::ofstream& fout);
    void recover(std::ifstream& fin);
};


/**
 * @class Cell @file Monomer.h/.cpp
 *
 * Cell class are used for constructing the structre of "bucket-list"
 * 
 */

typedef std::list<long>::iterator MonomerIndex;
typedef std::list<long>::iterator CellIndex;


class Cell {

public:
    
    long mNumberOfMonomers;
    long mNumberOfNeighbors;
    std::list< long > mMonomers; // atoms inside box
    std::list< long > mNeighbors; // neighboring boxes

    /**
     * constructor
     */
    Cell() : mNumberOfMonomers(0),mNumberOfNeighbors(0) {
    }

    /**
     * Copy constructor
     * @param orig
     */
    Cell(const Cell& orig);

    /**
     * Override operator =
     * @param orig
     * @return 
     */
    Cell& operator=(const Cell& orig);

    /**
     * Destructor
     */
    ~Cell();
    
    inline void addMonomer(const long& m){
        mMonomers.push_back(m);
        mNumberOfMonomers++;
    }

    inline void delMonomer(const MonomerIndex& mi){
        mMonomers.erase(mi);
        mNumberOfMonomers--;
    }

    inline void addNeighbor(const long& n){
        mNeighbors.push_back(n);
        mNumberOfNeighbors++;
    }
    
    inline void delNeighbor(const CellIndex& ci){
        mNeighbors.erase(ci);
        mNumberOfNeighbors--;
    }
    
    void print(std::ostream& out = std::cout);
    void backup(std::ofstream& fout);
    void recover(std::ifstream& fin);
};

#endif // MONOMER_H
