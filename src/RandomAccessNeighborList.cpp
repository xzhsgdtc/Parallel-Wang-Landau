#include "RandomAccessNeighborList.h"


RandomAccessNeighborList::RandomAccessNeighborList(const long guessNumberOfNeighbors){
    if (guessNumberOfNeighbors != -1){
        mNeighbors.reserve(guessNumberOfNeighbors);
    }else{
        mNeighbors.reserve(DEFAULT_NEIGHBOR_SIZE);
    }
}


RandomAccessNeighborList::RandomAccessNeighborList(const RandomAccessNeighborList& orig){
    mNeighbors = orig.mNeighbors;
}


RandomAccessNeighborList::~RandomAccessNeighborList(){
    mNeighbors.clear();
}

void RandomAccessNeighborList::addNeighbor(const long id){
    mNeighbors.push_back(id);
}
    
vector<long> RandomAccessNeighborList::randomPick(Random& ran, const long num){
    // currently use dumpest method to random pick

    vector<long> result;
    if (mNeighbors.empty()) return result;
    long size = mNeighbors.size();
    if (num == size) return mNeighbors;
    
    result.push_back(mNeighbors[ran.nextLong(0,size)]);
    long numberPicked(1), pick;
    bool flag(true);

    while (numberPicked < num){
        pick = mNeighbors[ran.nextLong(0,size)];
        flag = true;
        for(int i=0; i<numberPicked; ++i){
            if(pick == result[i]) {
                flag = false;
                break;
            }
        }    
        if(flag){
            result.push_back(pick);
            numberPicked ++;
        }
    }
    return result;
}

void RandomAccessNeighborList::print(ostream& out){
    long size = mNeighbors.size();
    for(int i=0; i<size; ++i)
        out << mNeighbors[i] << " ";
    out << endl;
}
