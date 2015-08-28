#ifndef RANDOMACCESSNEIGHBORLIST_H
#define RANDOMACCESSNEIGHBORLIST_H

#include <vector>
#include <iostream>
#include "Random.h"


using std::vector;
using std::endl;
using std::cout;
using std::ostream;

const long DEFAULT_NEIGHBOR_SIZE = 10;

class RandomAccessNeighborList{
    public:
        RandomAccessNeighborList(const long guessNumberOfNeighbors = -1);
        RandomAccessNeighborList(const RandomAccessNeighborList& orig);
        ~RandomAccessNeighborList();
        void addNeighbor(const long id);
        vector<long> randomPick(Random& ran, const long num);

        void print(ostream& out);
        long size(){
            return mNeighbors.size();
        }

    private:
        vector<long> mNeighbors;
};


#endif // RANDOMACCESSNEIGHBORLIST_H
