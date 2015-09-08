/* *********************************************************
 * :::::::::::::::::::::::::::::::::::::::::::::::::::::::::
 *
 *  This file contains the definition of template class
 *  MathVecotr, which is generic purpose class can be
 *  used in multiple situations.
 *
 *   @File:   MathVector.h/.cpp
 *   @Author: Guangjie Shi (Jerry)
 *
 *   @Version 1.0: Aug. 24, 2012
 *       Definition of class MathVector
 *
 *   @Version 2.0: Oct. 30, 2013
 *      Changed to template class, and make the periodic 
 *      boundary condition much clearer by separating
 *      it into two cases: 
 *          1. position periodic 2. distance periodic
 *   @Version 2.1: Feb. 5, 2015
 *      1. Corrected the function for calculating periodi
 *         boundary condition
 *      2. Instead of offering the method to generate min.
 *         image vector, this version directly calculate
 *         the min. image distance.
 *
 *
 * :::::::::::::::::::::::::::::::::::::::::::::::::::::::::
 */

#ifndef MATHVECTOR_H
#define	MATHVECTOR_H
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <cstdlib>

const double _default_equal_threshold_ = 1e-6; 

template<class T>
class MathVector {

private:
    
    long mDim; /**< dimension of vector */
    T* mCoord; /**< array to store mCoordinates of vector in cartesion mCoordinate */

public:    
    
    /**
     * Default copy constructor
     */
    MathVector();

    /**
     * Constructor
     * @param d: dimension of vector, default value is 3
     */
    MathVector(long d);
    
    /**
     * define a 3-D math vector
     * @param x
     * @param y
     * @param z
     */
    MathVector(T x, T y, T z);
    
    /**
     * define a 2-D math vector
     * @param orig
     */
    MathVector(T x, T y);
    
    /**
     * Copy constructor
     * @param orig
     */
    MathVector(const MathVector<T>& orig);
    
    /**
     * Destructor 
     */
     virtual ~MathVector();
    
    
     
    /*********************************************************************
     * :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
     * 
     * Methods of reading and writing
     * 
     */    

     /**
      *  @return value given dimension
      */
    T get(long d) { return mCoord[d];  }
    
    /**
     * @return dimensionality of vector
     */
    long dim() {return mDim;}
    
    /**
     * set the value of given dimension
     */
    void set(long d, T value) { mCoord[d] = value;  }
    
    /**
     * modify the value of given dimension
     */
    void add(long d, T value) { mCoord[d] += value; }
    
    
    
     
    /*********************************************************************
     *::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
     * 
     *  Methods of Math operation
     * 
     */
     
     /**
      * @return module of vector 
      */
    double module();
    
    /**
     * @return sum over square of each element 
     */
    T sumSquare();
    
    /**
     * Apply periodic boundary condition (PBC) for given dimension
     * @param d dimension
     * @param l left boundary coordinate
     * @param r right boundary coordinate
     * @return 
     */
    inline void imposePBC(const long& d, const T& l, const T& r);


    /**
     * Apply periodic boundary for all dimension
     * @param l left boundary coordinate
     * @param r right boundary coordinate
     */
    inline void imposePBC(const T& l, const T& r);


    /**
     * Calculate the distance between two vectors, based
     * on minimum image rule (Periodic Boundary Condition)
     */
    inline T minImageDist(const MathVector& m, const T& length); 
    

  //  /**
  //   * Apply periodic boundary conditions to calculate the minimum image vector
  //   * @param length length of the cubic box
  //   */
  //  inline void minImage(const T& length);

  //  /**
  //   * Apply periodic boundary conditions to calculate the minimum image vector for given dimension
  //   * @param d dimension index
  //   * @param length length of the box in given dimension
  //   */
  //  inline void minImage(const long& d, const T& length);

    /**
     * To test whether two vector are compatible
     * @param one
     * @return true if compatible, false if not
     */
    bool compatible(const MathVector& one);
    
    /**
     * Override the operator= for vector
     * @param one
     * @return 
     */
    MathVector& operator=(const MathVector& one);
    
    /**
     * 
     * @param a
     * @param equal threshold, default value is set if it is not given
     * @return true if equal, false if not
     */
    bool equal(const MathVector& a, const double& epsilon = _default_equal_threshold_);
    
    
    /**
     * Override the operator- for vector
     * @param a
     * @return result vector
     */
    MathVector operator-(const MathVector& a);
    
    /**
     * Override the operator+ for vector
     * @param aMathVector<T>.cpp:58:15: error: ‘exit’ was not declared in this scope
     * @param a
     * @return result vector
     */
    MathVector operator+(const MathVector& a);
    
    /**
     * Dot product of two vector 
     * @param a
     * @return result of dot product
     */
    T dot(const MathVector& a);
    
    /**
     * Cross product, only works for 3D vector (for now). 
     * @param a
     * @return result vector
     */
    MathVector cross(const MathVector& a);
    
    /**
     * Override the operator/ for vector
     * @param a
     * @return result vector
     */
    MathVector operator/(const MathVector& a);
    
    /**
     * Plus every element of vector by a constant a
     * @param constant a
     * @return result vector
     */
    MathVector operator+(const T& a);


    /**
     * Minus every element of vector by a constant a
     * @param constant a
     * @return result vector
     */
    MathVector operator-(const T& a);


    /**
     * Multiply every element of vector by a constant a
     * @param constant a
     * @return result vector
     */
    MathVector operator*(const T& a);


    /**
     * Divide every element of vector by a constant a
     * @param constant a
     * @return result vector
     */
    MathVector operator/(const T& a);


    /**
     * Pow every element of vector by a constant a
     * @param constant a
     * @return result vector
     */
    MathVector operator^(const T& a);


    /**
     * Print out vector
     * @param output stream
     * @param precision precision of ouput, default value is 15
     */
    void print(std::ostream& out = std::cout, long precision = 15);

    /**
     * Backup whole object into a binary file.
     * @param fout output file stream
     */
    void backup(std::ofstream& fout);

    /**
     * Read in a binary file, in order to recover.
     * @param fin input file stream
     */
    void recover(std::ifstream& fin);
    
};


template<class T>
MathVector<T>::MathVector():mDim(-1),mCoord(NULL){

}

template<class T>
MathVector<T>::MathVector(long d){
    mDim = d;
    mCoord = new T[mDim];
}

template<class T>
MathVector<T>::MathVector(T x, T y){
    mDim =2;
    mCoord = new T[mDim];
    mCoord[0]= x;
    mCoord[1]= y;
}

template<class T>
MathVector<T>::MathVector(T x, T y, T z){
    mDim =3;
    mCoord = new T[mDim];
    mCoord[0]= x;
    mCoord[1]= y;
    mCoord[2]= z;
}

template<class T>
MathVector<T>::MathVector(const MathVector& orig) {
    mDim = orig.mDim;
    if(mCoord != NULL) delete [] mCoord;
    if(orig.mCoord == NULL){
        mCoord = NULL;
    }else{
        mCoord = new T[mDim];
        for(long i=0; i<mDim; ++i){
            mCoord[i] = orig.mCoord[i];
        }
    }
}

template<class T>
MathVector<T>::~MathVector() {
    if(mCoord != NULL){
        delete [] mCoord;
    }
}

template<class T>
MathVector<T>& MathVector<T>::operator=(const MathVector<T>& one){
    
    if(this == &one) return *this;
    mDim = one.mDim;
    if(mCoord != NULL) delete [] mCoord;

    if(one.mCoord == NULL){
        mCoord = NULL;
    }else{
        mCoord = new T[mDim];
        for(long i=0; i<mDim; ++i){
            mCoord[i] = one.mCoord[i];
        }
    }
    return *this;
}

template<class T>
MathVector<T> MathVector<T>::operator+(const T& a){
    MathVector<T> c(mDim);
    for(long i=0; i<mDim; ++i){
        c.mCoord[i] = mCoord[i]+a;
    }
    return c;
}

template<class T>
MathVector<T> MathVector<T>::operator-(const T& a){
    MathVector<T> c(mDim);
    for(long i=0; i<mDim; ++i){
        c.mCoord[i] = mCoord[i] - a;
    }
    return c;
}

template<class T>
MathVector<T> MathVector<T>::operator*(const T& a){
    MathVector<T> c(mDim);
    for(long i=0; i<mDim; ++i){
        c.mCoord[i] = mCoord[i] * a;
    }
    return c;
}

template<class T>
MathVector<T> MathVector<T>::operator/(const T& a){
    MathVector<T> c(mDim);
    for(long i=0; i<mDim; ++i){
        c.mCoord[i] = mCoord[i] / a;
    }
    return c;
}

template<class T>
MathVector<T> MathVector<T>::operator^(const T& a){
    MathVector<T> c(mDim);
    for(long i=0; i<mDim; ++i){
        c.mCoord[i] = pow(mCoord[i],a);
    }
    return c;
}

template<class T>
MathVector<T> MathVector<T>::operator-(const MathVector& a){
    if(!this->compatible(a)){
        std::cout << "Error in MathVector<T> operation: Minus, two MathVector<T> are not compatible!" << std::endl;
        exit(1);
    }
    MathVector<T> c(a.mDim);
    for(long i=0; i<c.mDim; ++i){
        c.mCoord[i] = mCoord[i] - a.mCoord[i];
    }
    return c;
}

template<class T>
MathVector<T> MathVector<T>::operator+(const MathVector& a){
    if(!this->compatible(a)){
        std::cout << "Error in MathVector<T> operation: +, two MathVector<T> are not compatible!" << std::endl;
        exit(1);
    }
    MathVector<T> c(a.mDim);
    for(long i=0; i<c.mDim; ++i){
        c.mCoord[i] = mCoord[i] + a.mCoord[i];
    }
    return c;
}

template<class T>
T MathVector<T>::dot(const MathVector& a){
    if(!this->compatible(a)){
        std::cout << "Error in MathVector<T> operation: *, two MathVector<T> are not compatible!" << std::endl;
        exit(1);
    }
    T sum(0);
    for(long i=0; i<a.mDim; ++i){
        sum += mCoord[i] * a.mCoord[i];
    }
    return sum;
}

template<class T>
MathVector<T> MathVector<T>::operator/(const MathVector& a){
    if(!this->compatible(a)){
        std::cout << "Error in MathVector<T> operation: /, two MathVector<T> are not compatible!" << std::endl;
        exit(1);
    }
    MathVector<T> c(a.mDim);
    for(long i=0; i<c.mDim; ++i){
        c.mCoord[i] = (T)mCoord[i] / (T)a.mCoord[i];
    }
    return c;
}

template<class T>
MathVector<T> MathVector<T>::cross(const MathVector& a){
    if(!this->compatible(a)){
        std::cout << "Error in MathVector<T> operation: ^, two MathVector<T> are not compatible!" << std::endl;
        exit(1);
    }
    MathVector<T> c(a.mDim);
    if(c.mDim == 3){
        c.mCoord[0] = mCoord[1]*a.mCoord[2] - mCoord[2]*a.mCoord[1];
        c.mCoord[1] = mCoord[2]*a.mCoord[0] - mCoord[0]*a.mCoord[2];
        c.mCoord[2] = mCoord[0]*a.mCoord[1] - mCoord[1]*a.mCoord[0];
    }else{
        std::cout << "OOPs! Sorry we only provide 3 dimentional cross product !" << std::endl;
        exit(1);
    }
    return c;
}

template<class T>
bool MathVector<T>::equal(const MathVector& a, const double& epsilon){
    if(!this->compatible(a)){
        std::cout << "Error in MathVector<T> operation: ==, two MathVector<T> are not compatible!" << std::endl;
        exit(1);
    }
    for(long i=0; i<a.mDim; ++i){
        if( (mCoord[i] - a.mCoord[i]) > epsilon ) return false;
    }
    return true;
}

template<class T>
bool MathVector<T>::compatible(const MathVector& one){
    if(mDim != one.mDim) return false;
    return true;
}

template<class T>
double MathVector<T>::module(){
    T m = sumSquare();
    return sqrt(m);
}

template<class T>
T MathVector<T>::sumSquare(){
    T s(0);
    for(long d=0; d<mDim; ++d){
        s += pow(mCoord[d],2.0);
    }
    return s;
}

template<class T>
void MathVector<T>::imposePBC(const long& d, const T& l, const T& r){
    T L = r - l; // distance between left and right boundary 
    mCoord[d] = mCoord[d] - L * floor( (mCoord[d] - l)/L ); 
}

template<class T>
void MathVector<T>::imposePBC(const T& l, const T& r){
    T L = r - l; // distance between left and right boundary 
    for(int d=0; d<mDim; ++d){
        mCoord[d] = mCoord[d] - L * floor( (mCoord[d] - l)/L ); 
    }
}

template<class T>
T MathVector<T>::minImageDist(const MathVector& m, const T& length){
    if(!this->compatible(m)){
        std::cout << "Error in MathVector<T> function: minImageDist, two MathVector<T> are not compatible!" << std::endl;
        exit(1);
    }
    T dist(0), ad;
    T hl = length*0.5;
    for(int d=0; d<mDim; ++d){
       ad = fabs(mCoord[d] - m.mCoord[d]);  
       if(ad > hl) ad = length - ad;
       dist += pow(ad,2.0);
    }
    return sqrt(dist);
}

//template<class T>
//void MathVector<T>::minImage(const long& d, const T& length){
//    T hl = length*0.5;
//    if(mCoord[d] > hl){
//        mCoord[d] -= length;
//    }else if(mCoord[d] < -hl){
//        mCoord[d] += length;
//    }
//}
//
//template<class T>
//void MathVector<T>::minImage(const T& length){
//    T hl = length*0.5;
//    for(long d=0; d<mDim; ++d){
//        if(mCoord[d] > hl){
//            mCoord[d] -= length;
//        }else if(mCoord[d] < -hl){
//            mCoord[d] += length;
//        }
//    }
//}

template<class T>
void MathVector<T>::print(std::ostream& out, long precision){
    out << "[ ";
    for(long d=0; d<mDim-1; ++d){
        out << std::setprecision(precision) << mCoord[d] << ", ";
    }
    if(mDim > 0){
        out << std::setprecision(precision) << mCoord[mDim-1] ;
    }
    out << " ]" << std::endl;
}

template<class T>
void MathVector<T>::backup(std::ofstream& fout){
    fout.write(reinterpret_cast<char*>(mCoord), mDim*sizeof(T));
}

template<class T>
void MathVector<T>::recover(std::ifstream& fin){
    fin.read(reinterpret_cast<char*>(mCoord), mDim*sizeof(T));
}


#endif	/* MATHVECTOR_H */

