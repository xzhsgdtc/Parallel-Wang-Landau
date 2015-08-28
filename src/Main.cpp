
/* *********************************************************
 * :::::::::::::::::::::::::::::::::::::::::::::::::::::::::
 * 
 *   This file contains the Frame for running 
 *   Multi-Dimensional Wang-Landau sampling
 *
 *   @File:   PWLSimulator.cpp
 *   @Author: Guangjie Shi (Jerry)
 *
 *   @version 1.0: Nov. 18, 2013
 *      This version only support 1D and 2D of DOS
 *      
 * 
 * :::::::::::::::::::::::::::::::::::::::::::::::::::::::::
 */

#include <cstdlib>
#include <mpi.h>
#include "WLFrame.h"
#include "LipidModel.h"

typedef WLFrame<LipidModel> SimulationType;


enum INTER_WIN_COMM_TYPE{
    CROSS_COMM_TYPE = 0,
    SQUARE_COMM_TYPE
};

INTER_WIN_COMM_TYPE interWinCommType = CROSS_COMM_TYPE;   // cross type by default

long replicaExchangeCount[2] = {0,0};
long replicaExchangeAcceptCount[2] = {0,0};


// constants for helping retrieve process ids
// currently at most consider 2 dimension of density of states
const int PIDWORLDCOMM                  = 0;
const int PIDINTRACOMM                  = 1;
//const int PIDINTERCOMM_LEFT             = 2;
//const int PIDINTERCOMM_RIGHT            = 3;
//const int PIDINTERCOMM_UP               = 4;
//const int PIDINTERCOMM_DOWN             = 5;
//const int PIDINTERCOMM_LEFT_UP          = 6;
//const int PIDINTERCOMM_RIGHT_UP         = 7;
//const int PIDINTERCOMM_LEFT_DOWN        = 8;
//const int PIDINTERCOMM_RIGHT_DOWN       = 9;

const int INTERCOMM_LEFT                = 0;
const int INTERCOMM_RIGHT               = 1;
const int INTERCOMM_UP                  = 2;
const int INTERCOMM_DOWN                = 3;
const int INTERCOMM_LEFT_UP             = 4;
const int INTERCOMM_RIGHT_UP            = 5;
const int INTERCOMM_LEFT_DOWN           = 6;
const int INTERCOMM_RIGHT_DOWN          = 7;

int        status                       = 0;
int        dim                          = 0;
int        totalNumberOfProcs           = 0;
int        totalNumberOfWins            = 0;
int        numberOfProcsPerWin          = 0;
int        numberOfTotalInterWinComms   = 0;
int        numberOfInterCrossWinComms   = 0;
int        numberOfInterDiagWinComms    = 0;

int*       interCrossWinCommsEachDim    = NULL; // 0 (horizontal), 1 (vertical)
int*       interDiagWinCommsEachDim     = NULL; // 0 (left diagonal), 1 (right diagonal)
int*       winsEachDim                  = NULL;
int*       myIDs                        = NULL;
int*       myInterCrossWinCommIDs       = NULL;
int*       myInterDiagWinCommIDs        = NULL;
int*       myWinCoord                   = NULL;
double*    overlapEachDim               = NULL;

MPI_Comm   intraWinComm                 = NULL;
MPI_Group  worldGroup                   = NULL;
MPI_Comm*  interCrossWinComms           = NULL;
MPI_Comm*  interDiagWinComms            = NULL;
MPI_Group* interCrossWinGroups          = NULL;
MPI_Group* interDiagWinGroups           = NULL;


std::ofstream out;

// declare routines
void buildInterWinComms(); 
void buildInterCrossWinComms(); 
void buildInterDiagWinComms(); 
void buildIntraWinComm();
void cleanup();

void errorMsgAndQuit(std::string mes, int s){
    out << "Proc " << myIDs[PIDWORLDCOMM] << " encoutner Error @" << mes << std::endl;
    MPI_Abort(MPI_COMM_WORLD, s);
}

int main(int argc, char* argv[]){
    /**
     * argv[0]              : program name
     * argv[1]              : input file for Wang-Landau simulator
     * argv[2]              : n dimension of density of states
     * argv[3   ... 3+n-1 ] : No. of wins for each dimension
     * argv[3+n ... 3+2n-1] : overlap for each dimension
     * argv[3+2n]           : inter windows communication type
     * argv[3+2n+1]         : random seed for direction selector
     */

    if (argc < 7) {
        std::cout   << "Usage: " 
                    << std::string(argv[0]) 
                    << "  [input file]  [dimension of DOS]  [ ... wins for each dim ... ] [ ... overlap for dim ... ] [inter win. comm. type: 0 (cross), 1(square)] [random seed for direction selector]" 
                    << std::endl;
        exit(1);
    }

    int tid, tsize;
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &tsize);
    MPI_Comm_rank(MPI_COMM_WORLD, &tid);
    
    // open log file, append to the end
    std::string fname ="./log/P" + intToString(tid,7) + ".log"; 
    out.open(fname.c_str(), std::ios::app);
    if(!out.is_open()){
        errorMsgAndQuit("cannot open log file: " + fname, MPI_ERR_OTHER);
    }

    // processing input parameters
    out  << "\n# Parallel Frame for Wang-Landau sampling ... \n"
               << "  >> processing input parameters ... " << std::endl;
    
    std::string inputFilename = argv[1];
    dim = atoi(argv[2]);    // dimension of density of states
    winsEachDim     = new int[dim];
    overlapEachDim  = new double[dim];
    totalNumberOfWins = 1;
    for(int i=0; i<dim; ++i){
        winsEachDim[i] = atoi(argv[3+i]);
        overlapEachDim[i] = atof(argv[3+dim+i]);
        totalNumberOfWins *= winsEachDim[i];
    }

    myWinCoord = new int[dim];
    totalNumberOfProcs = tsize;
    numberOfProcsPerWin = totalNumberOfProcs / totalNumberOfWins; 

    long selector_seed = atoi(argv[3+2*dim+1]);
    Random random_selector(selector_seed);

    switch(dim){
        case 1:
            
            myIDs = new int[1+1+2*dim];   // world id, intra comm id, inter comm id
            myIDs[PIDWORLDCOMM] = tid;
            myInterCrossWinCommIDs = new int[2*dim];
            myWinCoord[0] = (int) myIDs[PIDWORLDCOMM]/numberOfProcsPerWin; 
            
            interWinCommType = CROSS_COMM_TYPE;
            interCrossWinCommsEachDim = new int[dim];
            interCrossWinCommsEachDim[0] = winsEachDim[0]-1;
            numberOfInterCrossWinComms = interCrossWinCommsEachDim[0];
            numberOfTotalInterWinComms = numberOfInterCrossWinComms;
            break;
        case 2:
            switch(atoi(argv[3+2*dim])){
                case 0: //CROSS TYPE INTER WIN. COMM.
                    
                    myIDs = new int[1+1+2*dim];   // world id, intra comm id, inter comm id
                    myIDs[PIDWORLDCOMM] = tid;
                    myInterCrossWinCommIDs = new int[2*dim];
                    myWinCoord[0] = (int) (myIDs[PIDWORLDCOMM]/numberOfProcsPerWin) % winsEachDim[0];
                    myWinCoord[1] = (int)  myIDs[PIDWORLDCOMM]/(numberOfProcsPerWin * winsEachDim[0]);
                    
                    interWinCommType = CROSS_COMM_TYPE;
                    interCrossWinCommsEachDim = new int[dim];
                    interCrossWinCommsEachDim[0] = (winsEachDim[0]-1)*winsEachDim[1];
                    interCrossWinCommsEachDim[1] = (winsEachDim[1]-1)*winsEachDim[0];
                    numberOfInterCrossWinComms = interCrossWinCommsEachDim[0] + interCrossWinCommsEachDim[1];
                    numberOfTotalInterWinComms = numberOfInterCrossWinComms;
                    break;

                case 1: //SQUARE TYPE INTER WIN. COMM.
                   
                    myIDs = new int[1+1+4*dim];   // world id, intra comm id, inter comm id
                    myIDs[PIDWORLDCOMM] = tid;
                    myInterCrossWinCommIDs = new int[2*dim];
                    myInterDiagWinCommIDs = new int[2*dim];
                    myWinCoord[0] = (int) (myIDs[PIDWORLDCOMM]/numberOfProcsPerWin) % winsEachDim[0];
                    myWinCoord[1] = (int)  myIDs[PIDWORLDCOMM]/(numberOfProcsPerWin * winsEachDim[0]);
                    
                    interWinCommType = SQUARE_COMM_TYPE;
                    interCrossWinCommsEachDim = new int[dim];
                    interDiagWinCommsEachDim = new int[dim];
                    interCrossWinCommsEachDim[0] = (winsEachDim[0]-1)*winsEachDim[1];
                    interCrossWinCommsEachDim[1] = (winsEachDim[1]-1)*winsEachDim[0];
                    interDiagWinCommsEachDim[0]  = (winsEachDim[0] -1)*(winsEachDim[1]-1);
                    interDiagWinCommsEachDim[1]  = (winsEachDim[0] -1)*(winsEachDim[1]-1);
                    numberOfInterCrossWinComms = interCrossWinCommsEachDim[0] + interCrossWinCommsEachDim[1];
                    numberOfInterDiagWinComms  = interDiagWinCommsEachDim[0] + interDiagWinCommsEachDim[1];
                    numberOfTotalInterWinComms = numberOfInterCrossWinComms + numberOfInterDiagWinComms;
                    break;
                default:
                    errorMsgAndQuit("illegal type of inter win. comm. ! ", MPI_ERR_OTHER);
            }
            break;
        default:
            errorMsgAndQuit("illegal dimensionality of density of states ! ", MPI_ERR_OTHER);
    }

    
    out << "     -> input file: " << inputFilename << "\n"
        << "     -> random seed for direction selector: " << selector_seed << "\n"
        << "     -> dimension of DOS: " << dim << "\n"
        << "     -> windows each dimension: ";
    for(int i=0; i<dim; ++i) out << winsEachDim[i] << "  ";
    out << "\n"
        << "     -> overlap each dimension: ";
    for(int i=0; i<dim; ++i) out << overlapEachDim[i] << "  ";
    out << "\n     -> number of processes: " << totalNumberOfProcs << "\n"
        << "     -> number of wins: " << totalNumberOfWins << "\n"
        << "     -> number of proc. per win.: " << numberOfProcsPerWin << "\n"
        << "     -> current procees id: " << myIDs[PIDWORLDCOMM] << "\n"
        << "     -> current win. coords.: ";
    for(int i=0; i<dim; ++i)
        out << myWinCoord[i] << " ";
    out << "\n     -> number of inter win. comm.: " << numberOfTotalInterWinComms << "  (";
    for(int i=0; i<dim; ++i){
        if (i != dim-1)
            out << interCrossWinCommsEachDim[i] << ", ";
        else 
            out << interCrossWinCommsEachDim[i];
    }
    if(interWinCommType == SQUARE_COMM_TYPE){
        for(int i=0; i<dim; ++i){
            if (i != dim-1)
                out << ", "<< interDiagWinCommsEachDim[i] << ", ";
            else 
                out << interDiagWinCommsEachDim[i];
        }
    }
    out << ")\n"
        << "        * inter win. comm. type: " << interWinCommType << "\n"
        << "  >> finish sucessfully! " << std::endl;

    MPI_Barrier(MPI_COMM_WORLD);
    
    // build inter-wins communicator
    status = MPI_Comm_group(MPI_COMM_WORLD, &worldGroup);
    
    if(status != MPI_SUCCESS){
        errorMsgAndQuit("MPI_Comm_group(MPI_COMM_WORLD, worldGroup)", status);
    }
    
    out << "  >> building communicators ... \n";
    out << "     -> intra wins. comms. ... " << std::endl;
    buildIntraWinComm();
    out << "     -> inter wins. comms. ... " << std::endl;
    buildInterWinComms();
    out << "  >> finish successfully!";
    MPI_Barrier(MPI_COMM_WORLD);

    // start simulation
    SimulationType* run  = new SimulationType(inputFilename.c_str(), myIDs[PIDWORLDCOMM]);
    
    run -> setOutput(out);
    run -> PWL_init(dim, winsEachDim, myWinCoord, overlapEachDim);
    

    long NOF = 3; // NOF = Number Of Files
    
    /** output file names
     *  0: tracking energy
     *  1: information for histogram
     *  2: information for each iteration
     */
    std::string filename[NOF];

    filename[0] = run->addPrefix("_Track_"    + intToString(run -> iteration(),4) + ".dat");
    filename[1] = run->addPrefix("_Histogram_"  + intToString(run -> iteration(),4) + ".dat");
    filename[2] = run->addPrefix("_Iteration_"  + intToString(run -> iteration(),4) + ".dat");

    std::ofstream fout[NOF];
    
    getOutputStream(fout[0], filename[0], true);

    int direction(-1), myDirection(-1), totalDirection(0);
    if(dim == 1){
        totalDirection = 2;
    }else{
        if(interWinCommType == SQUARE_COMM_TYPE){
            totalDirection = 8;
        }else{
            totalDirection = 4;
        }
    }
    
    do{

        
        if(run -> PWL_isReadyReplicaExchange()){
            // decide exchange direction
            switch(dim){
                case 1:
                    // randomly pick a direction to perform replica exchange
                    // 0 for left, 1 for right
                    direction = random_selector.nextLong(0, totalDirection);   
                    switch(direction){
                        case INTERCOMM_LEFT       : myWinCoord[0] % 2 == 0 ? myDirection = INTERCOMM_LEFT        : myDirection = INTERCOMM_RIGHT        ; break;
                        case INTERCOMM_RIGHT      : myWinCoord[0] % 2 == 0 ? myDirection = INTERCOMM_RIGHT       : myDirection = INTERCOMM_LEFT         ; break;
                    }
                    break;

                case 2:
                    // randomly pick a direction to perform replica exchange
                    // 0: left, 1: right, 2: up  3: down
                    // 4: left up, 5: right up, 6: left down, 7: right down
                    direction = random_selector.nextLong(0, totalDirection);  
                    switch(direction){
                        case INTERCOMM_LEFT       : myWinCoord[0] % 2 == 0 ? myDirection = INTERCOMM_LEFT        : myDirection = INTERCOMM_RIGHT        ; break;
                        case INTERCOMM_RIGHT      : myWinCoord[0] % 2 == 0 ? myDirection = INTERCOMM_RIGHT       : myDirection = INTERCOMM_LEFT         ; break;
                        case INTERCOMM_UP         : myWinCoord[1] % 2 == 0 ? myDirection = INTERCOMM_UP          : myDirection = INTERCOMM_DOWN         ; break;
                        case INTERCOMM_DOWN       : myWinCoord[1] % 2 == 0 ? myDirection = INTERCOMM_DOWN        : myDirection = INTERCOMM_UP           ; break;
                        case INTERCOMM_LEFT_UP    : myWinCoord[0] % 2 == 0 ? myDirection = INTERCOMM_LEFT_UP     : myDirection = INTERCOMM_RIGHT_DOWN   ; break;
                        case INTERCOMM_RIGHT_UP   : myWinCoord[0] % 2 == 0 ? myDirection = INTERCOMM_RIGHT_UP    : myDirection = INTERCOMM_LEFT_DOWN    ; break;
                        case INTERCOMM_LEFT_DOWN  : myWinCoord[0] % 2 == 0 ? myDirection = INTERCOMM_LEFT_DOWN   : myDirection = INTERCOMM_RIGHT_UP     ; break;
                        case INTERCOMM_RIGHT_DOWN : myWinCoord[0] % 2 == 0 ? myDirection = INTERCOMM_RIGHT_DOWN  : myDirection = INTERCOMM_LEFT_UP      ; break;
                    }
                    break;
            } 

            // propose replica exchange 
            if(direction <= 3){
                out << "      -> propose replica exchange myID "<< myIDs[2+myDirection] << ", myDirection:  " << myDirection << "   Cross COMM. No." << myInterCrossWinCommIDs[myDirection] << std::endl;
               if( run -> PWL_proposeReplicaExchange(numberOfProcsPerWin, myIDs[2+myDirection], myInterCrossWinCommIDs[myDirection], interCrossWinComms)){
                   replicaExchangeAcceptCount[0]++;
               }
                replicaExchangeCount[0]++;
                out << "         * Done!\n"  
                    << "         -> ( " << replicaExchangeAcceptCount[0] << " / " << replicaExchangeCount[0] << " ), ( " 
                                        << replicaExchangeAcceptCount[1] << " / " << replicaExchangeCount[1] <<" ) " << std::endl;
            }else{
                out << "      -> propose replica exchange myID "<< myIDs[2+myDirection] << ", myDirection:  " << myDirection  << "  Diag. COMM. No." << myInterDiagWinCommIDs[myDirection-4] << std::endl;
                if (run -> PWL_proposeReplicaExchange(numberOfProcsPerWin, myIDs[2+myDirection], myInterDiagWinCommIDs[myDirection-4], interDiagWinComms)){
                   replicaExchangeAcceptCount[1]++;
                }
                replicaExchangeCount[1]++;
                out << "         * Done!\n"  
                    << "         -> ( " << replicaExchangeAcceptCount[0] << " / " << replicaExchangeCount[0] << " ), ( " 
                                        << replicaExchangeAcceptCount[1] << " / " << replicaExchangeCount[1] <<" ) " << std::endl;

            }
        } else{
            run -> doMCMove();
            if(!run -> isMoveAccepted()){
                run -> undoMCMove();
            }
        }
        run -> updateHistogram();
        
        if(run -> isReadyPrint())   run -> print(fout[0], PRINT_WL_TRACKING_INFORMATION);

        if(run -> isReadyCheckHistogram()){
            if(run -> PWL_isHistogramFlat(&intraWinComm)){
                run -> PWL_dosMerge(&intraWinComm, numberOfProcsPerWin);
                
                getOutputStream(fout[1], filename[1]);
                getOutputStream(fout[2], filename[2]);

                run -> print(fout[1], PRINT_WL_HISTOGRAM);
                run -> print(fout[2], PRINT_WL_INFORMATION);
                run -> updateModificationFactor();
                run -> resetHistogram();

                for(int i=0; i<NOF; ++i) fout[i].close();

                filename[0] = run->addPrefix("_Track_"    + intToString(run -> iteration(),4) + ".dat");
                filename[1] = run->addPrefix("_Histogram_"  + intToString(run -> iteration(),4) + ".dat");
                filename[2] = run->addPrefix("_Iteration_"  + intToString(run -> iteration(),4) + ".dat");
                getOutputStream(fout[0], filename[0], true);
            }
        }


    } while(   ! (  run -> PWL_isReadyCheckComplete() 
                    &&  run -> PWL_isSimulationComplete() )  );


    MPI_Barrier(MPI_COMM_WORLD);
    out << "\n  >> cleanup and finiliaze system ... " << std::endl;
    cleanup();
    MPI_Finalize();
    out << "# Wang-Landau Simulation END!" << std::endl;
    return 0;
}


void buildInterWinComms() {
    out << "        <*>  inter cross wins. comms. " << std::endl;
    buildInterCrossWinComms();
    if(interWinCommType == SQUARE_COMM_TYPE){
        out << "        <*>  inter diag. wins. comms. " << std::endl;
        buildInterDiagWinComms();
    }
}

void buildInterDiagWinComms(){
    //::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    // Build Up Inter-Diagonal-Windows Communicators and Groups
    interDiagWinComms  = new MPI_Comm[numberOfInterDiagWinComms];
    interDiagWinGroups = new MPI_Group[numberOfInterDiagWinComms];

    int status;
    int* ranks;
    int numberOfProcsPerInterWinComm = 2*numberOfProcsPerWin;
    int w;
    ranks = new int[numberOfProcsPerInterWinComm];

    for(int d=0; d < dim; ++d){
        for(int i=0; i < interDiagWinCommsEachDim[d]; ++i){
            for (int j = 0; j < numberOfProcsPerInterWinComm; ++j) {

                if (d == 0) {
                    int m = i + 1 + i / (winsEachDim[0] - 1);
                    if(j < numberOfProcsPerWin){
                        ranks[j] = m * numberOfProcsPerWin + j;
                    }
                    else{
                        ranks[j] = ( m + winsEachDim[0]-2 ) * numberOfProcsPerWin + j;
                    }
                } else {
                    int m = i + i / (winsEachDim[0] - 1);
                    if(j < numberOfProcsPerWin){
                        ranks[j] = m * numberOfProcsPerWin + j;
                    }
                    else{
                        ranks[j] = ( m + winsEachDim[0] ) * numberOfProcsPerWin + j;
                    }
                }
            }
                
            w =  i + (d == 0 ?  0: interDiagWinCommsEachDim[d-1]);
            status = MPI_Group_incl(worldGroup, numberOfProcsPerInterWinComm, ranks, &interDiagWinGroups[w]);
            if(status != MPI_SUCCESS){
                errorMsgAndQuit("buildInterDiagWinComms() -> MPI_Group_incl(...)", status);
            }

            status = MPI_Comm_create(MPI_COMM_WORLD, interDiagWinGroups[w], &interDiagWinComms[w]);
            if(status != MPI_SUCCESS){
                errorMsgAndQuit("buildInterDiagWinComms() -> MPI_Comm_create(...)", status);
            }
            
            out << "            * interDiagWinComm " << std::setw(4)<< w << ":  ";
            for(int j=0; j < numberOfProcsPerInterWinComm; ++j){
                out << std::setw(5) << ranks[j] << " ";
            }
            out << std::endl;
        }
    }
    delete [] ranks;

    switch(dim){
        case 1:
            errorMsgAndQuit("dim = 1 cannot have SQUARE_COMM_TYPE!", MPI_ERR_OTHER);
            break;

        case 2:

            myInterDiagWinCommIDs[INTERCOMM_LEFT_UP       -4]  = ((myWinCoord[0] == 0                  || myWinCoord[1] == winsEachDim[1] -1) ? 
                                                                 -1 : (winsEachDim[0] - 1) * myWinCoord[1] + myWinCoord[0] - 1);

            myInterDiagWinCommIDs[INTERCOMM_RIGHT_UP      -4]  = ((myWinCoord[0] == winsEachDim[0] - 1 || myWinCoord[1] == winsEachDim[1] -1) ? 
                                                                 -1 : (winsEachDim[0] - 1) * myWinCoord[1] + myWinCoord[0]);

            myInterDiagWinCommIDs[INTERCOMM_LEFT_DOWN     -4]  = ((myWinCoord[0] == 0                  || myWinCoord[1] == 0) ? 
                                                                 -1 : (winsEachDim[0] - 1) * (myWinCoord[1]-1) + myWinCoord[0] - 1); 

            myInterDiagWinCommIDs[INTERCOMM_RIGHT_DOWN    -4]  = ((myWinCoord[0] == winsEachDim[0] - 1 || myWinCoord[1] == 0) ? 
                                                                 -1 : (winsEachDim[0] - 1) * (myWinCoord[1]-1) + myWinCoord[0]);


            if (myInterDiagWinCommIDs[INTERCOMM_RIGHT_UP  -4]   != -1)  myInterDiagWinCommIDs[INTERCOMM_RIGHT_UP - 4]  += interDiagWinCommsEachDim[0];
            if (myInterDiagWinCommIDs[INTERCOMM_LEFT_DOWN -4]   != -1)  myInterDiagWinCommIDs[INTERCOMM_LEFT_DOWN -4]  += interDiagWinCommsEachDim[0]; 

            //out << myInterDiagWinCommIDs[INTERCOMM_LEFT_UP       -4]  << "\t"
            //    << myInterDiagWinCommIDs[INTERCOMM_RIGHT_UP      -4]  << "\t" 
            //    << myInterDiagWinCommIDs[INTERCOMM_LEFT_DOWN     -4]  << "\t"
            //    << myInterDiagWinCommIDs[INTERCOMM_RIGHT_DOWN    -4]  <<  std::endl;            
    }

    for (int i = 0; i < 4; ++i) {
        if(myInterDiagWinCommIDs[i] != -1){
            status = MPI_Comm_rank(interDiagWinComms[myInterDiagWinCommIDs[i]], &myIDs[6+i]);
            if(status != MPI_SUCCESS){
                errorMsgAndQuit("buildInterDiagWinComms() -> MPI_Comm_rank(...)", status);
            }
        }else{
            myIDs[6+i] = -1;
        }
        out << "            * my id at interDiagWinComm "  << std::setw(4) << myInterDiagWinCommIDs[i] << ":  " << std::setw(5) << myIDs[6+i] << std::endl; 
    }

}

void buildInterCrossWinComms() {
    //::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    // Build Up Inter-Cross-Windows Communicators and Groups
   
    interCrossWinComms  = new MPI_Comm[numberOfInterCrossWinComms];
    interCrossWinGroups = new MPI_Group[numberOfInterCrossWinComms];

    int status;
    int* ranks;
    int numberOfProcsPerInterWinComm =2*numberOfProcsPerWin; 
    int w;
    ranks = new int[numberOfProcsPerInterWinComm];

    for(int d=0; d < dim; ++d){
        for(int i=0; i < interCrossWinCommsEachDim[d]; ++i){
            for (int j = 0; j < numberOfProcsPerInterWinComm; ++j) {

                if (d == 0) {
                    int m = i + i / (winsEachDim[d] - 1);
                    ranks[j] = m * numberOfProcsPerWin + j;
                } else {
                    int m = i;
                    if (j < numberOfProcsPerWin)
                        ranks[j] = m * numberOfProcsPerWin + j;
                    else
                        ranks[j] = (m + winsEachDim[0] - 1) * numberOfProcsPerWin + j;
                }

            }
                
            w =  i + (d == 0 ?  0: interCrossWinCommsEachDim[d-1]);
            status = MPI_Group_incl(worldGroup, numberOfProcsPerInterWinComm, ranks, &interCrossWinGroups[w]);
            if(status != MPI_SUCCESS){
                errorMsgAndQuit("buildInterCrossWinComms() -> MPI_Group_incl(...)", status);
            }

            status = MPI_Comm_create(MPI_COMM_WORLD, interCrossWinGroups[w], &interCrossWinComms[w]);
            if(status != MPI_SUCCESS){
                errorMsgAndQuit("buildInterCrossWinComms() -> MPI_Comm_create(...)", status);
            }
            
            out << "            * interCrossWinComm " << std::setw(4)<< w << ":  ";
            for(int j=0; j < numberOfProcsPerInterWinComm; ++j){
                out << std::setw(5) << ranks[j] << " ";
            }
            out << std::endl;
        }
    }
    delete [] ranks;
    
    switch(dim){
        case 1:
            myInterCrossWinCommIDs[INTERCOMM_LEFT]  = ((myWinCoord[0] == 0)                  ? -1 : myWinCoord[0]-1                                                  );
            myInterCrossWinCommIDs[INTERCOMM_RIGHT] = ((myWinCoord[0] == winsEachDim[0] - 1) ? -1 : myWinCoord[0]                                                    );
            break;

        case 2:
            myInterCrossWinCommIDs[INTERCOMM_LEFT]  = ((myWinCoord[0] == 0)                  ? -1 : (winsEachDim[0] - 1) *  myWinCoord[1]      + myWinCoord[0] - 1   );
            myInterCrossWinCommIDs[INTERCOMM_RIGHT] = ((myWinCoord[0] == winsEachDim[0] - 1) ? -1 : (winsEachDim[0] - 1) *  myWinCoord[1]      + myWinCoord[0]       );
            myInterCrossWinCommIDs[INTERCOMM_UP]    = ((myWinCoord[1] == winsEachDim[1] - 1) ? -1 :  winsEachDim[0]      *  myWinCoord[1]      + myWinCoord[0]       );
            myInterCrossWinCommIDs[INTERCOMM_DOWN]  = ((myWinCoord[1] == 0)                  ? -1 :  winsEachDim[0]      * (myWinCoord[1] - 1) + myWinCoord[0]       );
            
            if (myInterCrossWinCommIDs[INTERCOMM_UP]   != -1)  myInterCrossWinCommIDs[INTERCOMM_UP]    += interCrossWinCommsEachDim[0]; 
            if (myInterCrossWinCommIDs[INTERCOMM_DOWN] != -1)  myInterCrossWinCommIDs[INTERCOMM_DOWN]  += interCrossWinCommsEachDim[0];
    }

    for (int i = 0; i < 2*dim; ++i) {
        if(myInterCrossWinCommIDs[i] != -1){
            status = MPI_Comm_rank(interCrossWinComms[myInterCrossWinCommIDs[i]], &myIDs[2+i]);
            if(status != MPI_SUCCESS){
                errorMsgAndQuit("buildInterCrossWinComms() -> MPI_Comm_rank(...)", status);
            }
        }else{
            myIDs[2+i] = -1;
        }
        out << "            * my id at interCrossWinComm "  << std::setw(4) << myInterCrossWinCommIDs[i] << ":  " << std::setw(5) << myIDs[2+i] << std::endl; 
    }

}

void buildIntraWinComm() {
    //::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    // Build Up Intra-Win Communicator
    int color = (int) myIDs[PIDWORLDCOMM] / numberOfProcsPerWin;
    MPI_Comm_split(MPI_COMM_WORLD, color, 0, &intraWinComm);
    MPI_Comm_rank(intraWinComm, &myIDs[PIDINTRACOMM]);
    out << "            * color of intra win. comm.: " << color << std::endl;
}


void cleanup(){
    delete []  interCrossWinCommsEachDim;   
    delete []  interDiagWinCommsEachDim ;   
    delete []  winsEachDim              ;  
    delete []  myIDs                    ;  
    delete []  myInterCrossWinCommIDs   ;  
    delete []  myInterDiagWinCommIDs    ;  
    delete []  myWinCoord               ;  
    delete []  overlapEachDim           ;  
    for(int i=0; i<numberOfInterCrossWinComms; ++i) {
        if (interCrossWinComms[i] != MPI_COMM_NULL){
            out << i << std::endl;
            MPI_Comm_free(&interCrossWinComms[i]);
            MPI_Group_free(&interCrossWinGroups[i]);
        }
    }
    for(int i=0; i<numberOfInterDiagWinComms; ++i) {
        if (interDiagWinComms[i] != MPI_COMM_NULL){
            out << i << std::endl;
            MPI_Comm_free(&interDiagWinComms[i]);
            MPI_Group_free(&interDiagWinGroups[i]);
        }
    }
    delete [] interCrossWinComms;
    delete [] interCrossWinGroups;
    delete [] interDiagWinComms;
    delete [] interDiagWinGroups;
    MPI_Comm_free(&intraWinComm);
}

