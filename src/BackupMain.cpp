
/* *********************************************************
 * :::::::::::::::::::::::::::::::::::::::::::::::::::::::::
 * 
 *   File:   Main.cpp
 *   Author: Guangjie Shi (Jerry)
 *
 *   Created on March 30, 2013, 1:58 PM
 * 
 *   Description:
 *      This file contains the Frame for running 
 *      Multi-Dimensional Wang-Landau sampling
 * 
 * :::::::::::::::::::::::::::::::::::::::::::::::::::::::::
 */


#include <mpi.h>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <sstream>
#include "WLFrame.h"
#include "LipidModel.h"

typedef WLFrame<LipidModel> SimulationType;

const int DIM = 2;

// Dim: 0
const int Left = 0;
const int Right = 1;

// Dim: 1
const int Up = 2;
const int Down = 3;

int num_procs; // total number of processes in simulation
int num_procs_per_win; // number of processes per window
int pid_line;
int pid_row;
int pid_world_comm; // process id in MPI_Comm_World
int pid_inter_win_comm[4]; // process id in inter-win communicator
int pid_intra_win_comm; // process id in intra-win communicator
int inter_win_comm_id[4]; // window id of inter-win communicator
int intra_win_comm_id; // window id of intra-win communicator

//::::::::::::::::::::::::::::::::::::::::::::::::::::
// 0: horizontal direction (Energy)
// 1: vertical direction (No. of Lipids)
int num_wins_dim[DIM];
int num_inter_win_comm[DIM];
double overlap[DIM];

MPI_Comm **mpi_inter_win_comm;
MPI_Group **mpi_inter_win_group;
MPI_Comm *mpi_intra_win_comm;
MPI_Group world_group;
MPI_Status status;


void build_inter_win_comm();
void build_intra_win_comm();

int opposite(int direc){
    int result = direc;
    switch(direc){
        case Left: result = Right; break;
        case Right: result = Left; break;
        case Up: result = Down; break;
        case Down: result = Up; break;
    }
    return result;
}

string direction_tostring(int direc){
    string result;
    switch(direc){
        case Left:  result = "Left" ; break;
        case Right: result = "Right"; break;
        case Up:    result = "Up"   ; break;
        case Down:  result = "Down" ; break;
    }
    return result;
}


int main(int argc, char* argv[]) {

    /**
     * argv[0] : program name
     * argv[1] : input file for Wang-Landau simulator
     * argv[2] : No. of wins for first dimension (Energy dimension)
     * argv[3] : No. of wins for second dimension (Lipid No. dimension)
     * argv[4] : overlap of first dimension
     * argv[5] : overlap of second dimension
     */
    if (argc != 6) {
        cout << "Input error, format should be [Pro name] [input_wanglandau.txt]" << endl;
        exit(1);
    }

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &num_procs);
    MPI_Comm_rank(MPI_COMM_WORLD, &pid_world_comm);
   

    //::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    // Deal with arguments
    string inputfile = argv[1];
    num_wins_dim[0] = atoi(argv[2]);
    num_wins_dim[1] = atoi(argv[3]);
    overlap[0] = atof(argv[4]);
    overlap[1] = atof(argv[5]);

    num_inter_win_comm[0] = (num_wins_dim[0] - 1) * num_wins_dim[1];
    num_inter_win_comm[1] = (num_wins_dim[1] - 1) * num_wins_dim[0];
    num_procs_per_win = num_procs / (num_wins_dim[0] * num_wins_dim[1]);
    pid_line = (int) (pid_world_comm / (num_wins_dim[0] * num_procs_per_win));
    pid_row = (int) ((pid_world_comm / num_procs_per_win) % num_wins_dim[0]);
 
    MPI_Barrier(MPI_COMM_WORLD);

    //::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    // output information
    if (pid_world_comm == 0) {
        cout << "\n" 
             << "##########  BEGIN: Structure Info.  ##########\n\n"
             << setw(30) << "No. of Procs:  "
             << setw(10) << num_procs << endl
             << setw(30) << "No. of Procs Per Win:  "
             << setw(10) << num_procs_per_win << endl
             << setw(30) << "No. of Wins:  ";
        for (int d = 0; d < DIM; ++d) {
            cout << setw(10) << "Dim-" << d << ": " << num_wins_dim[d] << "\t";
        }
        cout << endl;
        cout << setw(30) << "No. of Inter-Wins Comm:  ";
        for (int d = 0; d < DIM; ++d) {
            cout << setw(10) << "Dim-" << d << ": " << num_inter_win_comm[d] << "\t";
        }
        cout << "\n\n" << endl;
    }

    //::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    // build up communicators
    MPI_Comm_group(MPI_COMM_WORLD, &world_group);
    if(pid_world_comm == 0){
        cout << "    Build up inter-windows communication ... ";
    }

    build_inter_win_comm();
    
    if(pid_world_comm == 0){
        cout << "DONE!" << endl;
    }

    if(pid_world_comm == 0){
        cout << "    Build up intra-window communication ... ";
    }

    build_intra_win_comm();
    
    if(pid_world_comm == 0){
        cout << "DONE!" << endl;
    }
    
    //::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    // output structure infomation
    const int len_pBuffer = 13;
    int pBuffer[len_pBuffer];
    pBuffer[0] = pid_world_comm;
    pBuffer[1] = pid_line;
    pBuffer[2] = pid_row;
    pBuffer[3] = pid_inter_win_comm[0];
    pBuffer[4] = pid_inter_win_comm[1];
    pBuffer[5] = pid_inter_win_comm[2];
    pBuffer[6] = pid_inter_win_comm[3];
    pBuffer[7] = inter_win_comm_id[0];
    pBuffer[8] = inter_win_comm_id[1];
    pBuffer[9] = inter_win_comm_id[2];
    pBuffer[10]= inter_win_comm_id[3];
    pBuffer[11]= pid_intra_win_comm;
    pBuffer[12]= intra_win_comm_id;
    if(pid_world_comm == 0){
       cout << "\n"
            << "    Proc-"  << pBuffer[0]  << " ( L: " << pBuffer[1] << ", R: " << pBuffer[2] << " )\n"
            << "      Intra     " << setw(5)<< pBuffer[11] << setw(5)<< "   (AT) " << setw(5)<< pBuffer[12] << "\n"
            << "      Inter-D0  " << setw(5)<< pBuffer[3]  << setw(5)<< "   (AT) " << setw(5)<< pBuffer[7] 
            << "   (AND) "            << setw(5)<< pBuffer[4]  << setw(5)<< "   (AT) " << setw(5)<< pBuffer[8]  <<"\n"
            << "      Inter-D1  " << setw(5)<< pBuffer[5]  << setw(5)<< "   (AT) " << setw(5)<< pBuffer[9] 
            << "   (AND) "            << setw(5)<< pBuffer[6]  << setw(5)<< "   (AT) " << setw(5)<< pBuffer[10] <<"\n";    
       for(int i=1; i<num_procs; ++i){
            MPI_Recv(pBuffer, len_pBuffer, MPI_INT, i, 0, MPI_COMM_WORLD, &status );
            cout << "\n"
                 << "    Proc-"  << pBuffer[0]  << " ( L: " << pBuffer[1] << ", R: " << pBuffer[2] << " )\n"
                 << "      Intra     " << setw(5)<< pBuffer[11] << setw(5)<< "   (AT) " << setw(5)<< pBuffer[12] << "\n"
                 << "      Inter-D0  " << setw(5)<< pBuffer[3]  << setw(5)<< "   (AT) " << setw(5)<< pBuffer[7] 
                 << "   (AND) "            << setw(5)<< pBuffer[4]  << setw(5)<< "   (AT) " << setw(5)<< pBuffer[8]  <<"\n"
                 << "      Inter-D1  " << setw(5)<< pBuffer[5]  << setw(5)<< "   (AT) " << setw(5)<< pBuffer[9] 
                 << "   (AND) "            << setw(5)<< pBuffer[6]  << setw(5)<< "   (AT) " << setw(5)<< pBuffer[10] <<"\n";    
       }
    }else{
        MPI_Send(pBuffer, len_pBuffer, MPI_INT, 0, 0, MPI_COMM_WORLD); 
    }

  /*  //::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    // output structure infomation
    MPI_File log_file;
    MPI_File_open(MPI_COMM_WORLD, (char*)logfilename.c_str(), MPI_MODE_RDWR | MPI_MODE_CREATE, MPI_INFO_NULL, &log_file);
    string mes = toString();
    int size_mes = 350*sizeof(char)*pid_world_comm;
    MPI_File_set_view(log_file, size_mes, MPI_CHAR,MPI_CHAR, "native", MPI_INFO_NULL);
    MPI_File_write(log_file,(char*)mes.c_str(), mes.length(), MPI_CHAR, &status);
    MPI_File_close(&log_file);
    */

    //::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    // Monte Carlo initialization
    if(pid_world_comm == 0)
        cout << "\n##########  END: Structure Info.  ##########\n\n"
             << "\n### Simulation Starts ... \n"
             << "    Initialize Monte Carlo Simulator ... ";

    SimulationType* mcrun = new SimulationType(inputfile, pid_world_comm);
    int id[DIM] = {pid_row, pid_line};
    mcrun->PWL_init(DIM, num_wins_dim, id,overlap);
   
    MPI_Barrier(MPI_COMM_WORLD);
    
    if(pid_world_comm == 0)
        cout << "DONE!" << endl;

    
    long NOF = 3; // NOF = Number Of Files
    
    /** output file names
     *  0: tracking energy
     *  1: information for histogram
     *  2: information for each iteration
     */
    string filename[NOF];
    string his_curr_file;

    filename[0] = mcrun->add_prefix("energy_track_" + ToString::intToString(mcrun -> getIterCount(),7) + ".dat");
    filename[1] = mcrun->add_prefix("histogram_data_" + ToString::intToString(mcrun -> getIterCount(),7)+".dat");
    filename[2] = mcrun->add_prefix("iteration_data_" + ToString::intToString(mcrun -> getIterCount(),7)+".dat");
    
    ofstream fout[NOF];
    InputOutput::getOutputStream(fout[0],filename[0],true);

    double poss;
    int direction = 0;
    do {

        if (mcrun -> doMCMove() != -1) {
            poss = Random::nextDouble(0, 1);
            if (poss > mcrun -> prob_accept()) {
                mcrun -> undoMCMove();
            }
        }

        mcrun -> his_update();

        mcrun -> track_print(fout[0]);

        if (mcrun -> ready_check_his()) {
            mcrun -> myout << "Check histogram flatness ... ";
            if (mcrun->PWL_his_flat(mpi_intra_win_comm)) {
                mcrun -> myout << "    histogram flat ! \n";
                mcrun->PWL_merge_dos(mpi_intra_win_comm, num_procs_per_win);
                InputOutput::getOutputStream(fout[1],filename[1],true);
                InputOutput::getOutputStream(fout[2],filename[2],true);
               
                mcrun -> myout  << "     -> merge dos \n"
                                << "     -> print histogram to file "           << filename[1] << " \n"
                                << "     -> print current parameters to file "  << filename[2] << "\n"
                                << "    Iteration " << mcrun -> getIterCount() << " finished !" << endl;

                // output histogram file
                mcrun -> his_print(fout[1]);
                // output iteration records to file
                mcrun -> print(fout[2]);

                mcrun -> mf_update();
                mcrun -> his_reset();
                                
                
                for(long i=0; i<NOF; ++i){
                    fout[i].close();
                    fout[i].clear();
                }

                filename[0] = mcrun->add_prefix("energy_track_" + ToString::intToString(mcrun -> getIterCount(),7) + ".dat");
                filename[1] = mcrun->add_prefix("histogram_data_" + ToString::intToString(mcrun -> getIterCount(),7)+".dat");
                filename[2] = mcrun->add_prefix("iteration_data_" + ToString::intToString(mcrun -> getIterCount(),7)+".dat");

                InputOutput::getOutputStream(fout[0],filename[0],true);
            } else{
                mcrun -> myout << "criterion is not satisfied, continue! " << endl;
                ofstream tfout;
                his_curr_file = mcrun->add_prefix("histogram_curr_" + ToString::intToString(mcrun -> getIterCount(),7)+".dat");
                InputOutput::getOutputStream(tfout,his_curr_file);
                mcrun->his_print(tfout);
                tfout.close();
            }
        }

        if (mcrun->PWL_ready_replica_exchange()) {
            int my_direction(-1);
            int dim(-1);
            switch(direction){
                case Left:
                case Right:
                    dim=0;
                    if(pid_row % 2 == 0) my_direction  = opposite(direction); //even
                    else my_direction = direction; //odd
                    break;
                case Up:
                case Down:
                    dim=1;
                    if(pid_line % 2 == 0) my_direction  = direction; //even
                    else my_direction = opposite(direction); //odd
                    break;
                default:
                    cout << "Illegal direction !" << endl;
                    exit(1);
            }
            
          /*  mcrun->myout << "Synchronizing ... ";
            MPI_Barrier(MPI_COMM_WORLD);
            mcrun->myout << "DONE!" << endl;*/
            mcrun -> myout << "Propose replica exchange ... ("
                           << "Dim: "  << setw(4) << dim
                           << ",Direct: " << setw(4) << direction
                           << ",MyDirect: " << setw(4) << my_direction
                           << ",InterCommId:  " << setw(4) << inter_win_comm_id[my_direction]
                           << ",InterPid:  " << setw(4) << pid_inter_win_comm[my_direction]
                           << ") ... ";
            mcrun->PWL_propose_replica_exchange(    num_procs_per_win                    ,
                                                    pid_inter_win_comm[my_direction]     ,
                                                    inter_win_comm_id[my_direction]      ,    
                                                    mpi_inter_win_comm[dim]
                                                );
            
          /*  mcrun->myout << "Synchronizing ... ";
            MPI_Barrier(MPI_COMM_WORLD);
            mcrun->myout << "DONE!" << endl; */
            mcrun -> myout << "DONE!" << endl;
            direction = (direction + 1) % 4;
        }

        if(mcrun->PWL_ready_check_complete()){
            mcrun -> myout << "Check Complete Criterion ... " ;
            if(mcrun->PWL_complete()) {
                mcrun -> myout << "COMPLETE!" << endl;   
                break;
            } else{
                mcrun -> myout << "Not complete, CONTINUE!" << endl;
            }
        }

    } while (1);
    

    //:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    // clean up program
    for(long i=0; i<NOF; ++i){
        fout[i].close();
    }

    MPI_Barrier(MPI_COMM_WORLD);
    delete mcrun;
    if(pid_world_comm == 0)
        cout << "\n### Simulation Ends ! \n";
    MPI_Finalize();
    return 0;
}

void build_inter_win_comm() {
    //::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    // Build Up Inter-Windows Communicators and Groups
    int* ranks;
    ranks = new int[2 * num_procs_per_win];
    mpi_inter_win_comm = new MPI_Comm*[DIM];
    mpi_inter_win_group = new MPI_Group*[DIM];
    for (int i = 0; i < DIM; ++i) {
        mpi_inter_win_comm[i] = new MPI_Comm[num_inter_win_comm[i]];
        mpi_inter_win_group[i] = new MPI_Group[num_inter_win_comm[i]];
    }

    for (int d = 0; d < DIM; ++d) {
        for (int i = 0; i < num_inter_win_comm[d]; ++i) {
            for (int j = 0; j < 2 * num_procs_per_win; ++j) {

                if (d == 0) {
                    int m = i + i / (num_wins_dim[d] - 1);
                    ranks[j] = m * num_procs_per_win + j;
                } else {
                    int m = i;
                    if (j < num_procs_per_win)
                        ranks[j] = m * num_procs_per_win + j;
                    else
                        ranks[j] = (m + num_wins_dim[0] - 1) * num_procs_per_win + j;
                }
            }
            MPI_Group_incl(world_group, 2 * num_procs_per_win, ranks, &mpi_inter_win_group[d][i]);
            MPI_Comm_create(MPI_COMM_WORLD, mpi_inter_win_group[d][i], &mpi_inter_win_comm[d][i]);
        }

    }
    delete [] ranks;

    //::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    // Get Own ID in Inter-Win Communicator
    int line, row;
    line = (int) (pid_world_comm / (num_wins_dim[0] * num_procs_per_win));
    row = (int) ((pid_world_comm / num_procs_per_win) % num_wins_dim[0]);

    inter_win_comm_id[Right]  = ((row == num_wins_dim[0] - 1) ? -1 : (num_wins_dim[0] - 1) * line + row);
    inter_win_comm_id[Left]   = ((row == 0) ? -1 : (num_wins_dim[0] - 1) * line + row - 1);
    inter_win_comm_id[Up]     = ((line == num_wins_dim[1] - 1) ? -1 : num_wins_dim[0] * line + row);
    inter_win_comm_id[Down]   = ((line == 0) ? -1 : num_wins_dim[0] * (line - 1) + row);

    for (int i = 0; i < pow(2.0, DIM); ++i) {
        if (inter_win_comm_id[i] != -1) {
            MPI_Comm_rank(mpi_inter_win_comm[i / DIM][inter_win_comm_id[i]], &pid_inter_win_comm[i]);
        } else {
            pid_inter_win_comm[i] = -1;
        }
    }
    MPI_Barrier(MPI_COMM_WORLD);
    
}

void build_intra_win_comm() {
    //::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    // Build Up Intra-Win Communicator
    intra_win_comm_id = (int)pid_world_comm / num_procs_per_win;
    mpi_intra_win_comm = new MPI_Comm();
    MPI_Comm_split(MPI_COMM_WORLD, intra_win_comm_id, 0, mpi_intra_win_comm);
    MPI_Comm_rank(*mpi_intra_win_comm, &pid_intra_win_comm);
    MPI_Barrier(MPI_COMM_WORLD);
}


