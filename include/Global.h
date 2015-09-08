/* *********************************************************
 * :::::::::::::::::::::::::::::::::::::::::::::::::::::::::
 * 
 *  This file contain some global functions 
 *
 *   @File:   Global.h/.cpp
 *   @Author: Guangjie Shi (Jerry)
 *
 *   @Version 1.0: Oct. 12, 2012
 *      definition of two classes Monomer and Cell
 *
 *   @Version 2.0: Nov.  1. 2013
 *      clear out all the definition of classes, and only
 *      include some globally used functions. This file 
 *      does not need to consider about backup and recover
 * 
 * :::::::::::::::::::::::::::::::::::::::::::::::::::::::::
 */

#ifndef GLOBAL_H
#define	GLOBAL_H

#include <cstdlib>
#include <iomanip>
#include <string>
#include <iostream>
#include <fstream>
#include "Monomer.h"

const std::string OUTPUT_DIR="./record/";   // directory for outputing result files
const std::string INPUT_DIR="./input/"  ;   // directory for input files

enum PrintType{
    PRINT_MODEL_INFORMATION = 0,
    PRINT_MODEL_CELL,
    PRINT_MODEL_XYZ,
    PRINT_MODEL_MOL2,
    PRINT_MODEL_STATISTICS, 
    PRINT_MODEL_LIPIDWATER_ARRAY,
    PRINT_MODEL_OBSERVABLES,
    PRINT_MODEL_OBSERVABLES_ONLY,
    PRINT_WL_HISTOGRAM,
    PRINT_WL_INFORMATION,
    PRINT_WL_TRACKING_INFORMATION,
    PRINT_MAX_MOVE_DISTANCE
};

/**
 * convert int type to string type
 * @param number
 * @param digits, number of digits (eg. 00001 digits = 5)
 * @return according string value
 */
std::string intToString(const int& number, const int& digits = 0);

/**
 * convert double type to string type
 * @param number
 * @param after number of digits after decimal point "%.(after)f"
 * @return according string value
 */
std::string doubleToString(const double& number, const int& after );

/**
 * convert MonomerType to string type
 * @param a MonomerType
 * @return according string value
 */
std::string MTToString(MonomerType a);

/**
 * get an output file stream object
 * @param fout output file stream
 * @param filename 
 * @param app flag to decide whether append to the end of file
 */
void getOutputStream(std::ofstream& fout, std::string filename, bool app=false);
    
/**
 * get an input file stream object
 * @param fin input file stream
 * @param filename
 */
void getInputStream(std::ifstream& fin, std::string filename);


#endif	/* GLOBAL_H */

