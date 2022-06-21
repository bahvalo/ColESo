#pragma once
#ifndef BASE_CONTAINERS_HEADER
#define BASE_CONTAINERS_HEADER

#include "base_lib.h"

//======================================================================================================================
// Fixed-size data block - unallocatable, contains its data inside, gives explicit access check
// Intended to be used with STL vectors vector<tFixBlock<> > when it is needed to store vector of fixed-size blocks
// - fixed size given by template. no allocations.
// - explicit access check at safe mode
//======================================================================================================================
template <typename T, unsigned int N> class tFixBlock {
private:
    typedef tFixBlock<T,N> BaseType;
    T V[N]; // data array

public:
    inline tFixBlock(){ for(unsigned int i=0; i<N; i++) V[i] = 0; }
    inline tFixBlock(const BaseType &b){for(unsigned int i=0; i<N; i++) V[i]=b.V[i];}
    inline tFixBlock(T  v){ for(unsigned int i=0; i<N; i++) V[i]=v;}
    inline tFixBlock(T *v){ for(unsigned int i=0; i<N; i++) V[i]=v[i];}
    inline~tFixBlock(){}

    inline operator       T*()       {return V;}
    inline operator const T*() const {return V;}

    inline int Size() const {return N;}

    inline       T& operator[](int i){     SAFE_ASSERT((i>=0&&i<(int)N),"tFixBlock: wrong index [%d] %d",i,N); return V[i];}
    inline const T& operator[](int i)const{SAFE_ASSERT((i>=0&&i<(int)N),"tFixBlock: wrong index [%d] %d",i,N); return V[i];}
    inline bool eq(const BaseType &b) const{
        for(unsigned int i=0; i<N; i++) if(V[i]!=b.V[i]) return false;
        return true;
    }
    inline bool gt(const BaseType &b) const{
        for(unsigned int i=0; i<N; i++){
            if(V[i]>b.V[i]) return true;
            if(V[i]<b.V[i]) return false;
        }
        return false;
    }
    inline bool ge(const BaseType &b) const{
        for(unsigned int i=0; i<N; i++){
            if(V[i]>b.V[i]) return true;
            if(V[i]<b.V[i]) return false;
        }
        return true;
    }
    inline bool lt(const BaseType &b) const{
        for(unsigned int i=0; i<N; i++){
            if(V[i]<b.V[i]) return true;
            if(V[i]>b.V[i]) return false;
        }
        return false;
    }
    inline bool le(const BaseType &b) const{
        for(unsigned int i=0; i<N; i++){
            if(V[i]<b.V[i]) return true;
            if(V[i]>b.V[i]) return false;
        }
        return true;
    }
    inline BaseType& operator=(T x[]){for(unsigned int i=0; i<N; i++) V[i] = x[i]; return *this; }
    inline BaseType& operator=( T x){ for(unsigned int i=0; i<N; i++) V[i] = x; return *this; }
    inline BaseType& operator*=(T x){ for(unsigned int i=0; i<N; i++) V[i]*= x; return *this; }
    inline BaseType& operator+=(T x){ for(unsigned int i=0; i<N; i++) V[i]+= x; return *this; }
    inline BaseType& operator-=(T x){ for(unsigned int i=0; i<N; i++) V[i]-= x; return *this; }
    inline BaseType& operator+=(const BaseType& o){ for(unsigned int i=0; i<N; i++) V[i] += o.V[i]; return *this; }
    inline BaseType& operator-=(const BaseType& o){ for(unsigned int i=0; i<N; i++) V[i] -= o.V[i]; return *this; }
    inline BaseType& operator*=(const BaseType& o){ for(unsigned int i=0; i<N; i++) V[i] *= o.V[i]; return *this; }
};

template <typename T, unsigned int N>
inline bool operator==(const tFixBlock<T,N> &a, const tFixBlock<T,N> &b){return a.eq(b);}
template <typename T, unsigned int N>
inline bool operator>(const tFixBlock<T,N> &a, const tFixBlock<T,N> &b){ return a.gt(b);}
template <typename T, unsigned int N>
inline bool operator>=(const tFixBlock<T,N> &a, const tFixBlock<T,N> &b){return a.ge(b);}
template <typename T, unsigned int N>
inline bool operator<(const tFixBlock<T,N> &a, const tFixBlock<T,N> &b){ return a.lt(b);}
template <typename T, unsigned int N>
inline bool operator<=(const tFixBlock<T,N> &a, const tFixBlock<T,N> &b){return a.le(b);}

template<typename T, unsigned int N> 
inline tFixBlock<T,N> operator*(const tFixBlock<T,N> &a, const tFixBlock<T,N> &b){tFixBlock<T,N> R(a); return R *= b;}
template<typename T, unsigned int N> 
inline tFixBlock<T,N> operator+(const tFixBlock<T,N> &a, const tFixBlock<T,N> &b){tFixBlock<T,N> R(a); return R += b;}
template<typename T, unsigned int N> 
inline tFixBlock<T,N> operator-(const tFixBlock<T,N> &a, const tFixBlock<T,N> &b){tFixBlock<T,N> R(a); return R -= b;}
 
template<typename T, unsigned int N> 
inline tFixBlock<T,N> operator*(const tFixBlock<T,N> &a, T b){tFixBlock<T,N> R(a); return R *= b;}
template<typename T, unsigned int N> 
inline tFixBlock<T,N> operator+(const tFixBlock<T,N> &a, T b){tFixBlock<T,N> R(a); return R += b;}
template<typename T, unsigned int N> 
inline tFixBlock<T,N> operator-(const tFixBlock<T,N> &a, T b){tFixBlock<T,N> R(a); return R -= b;}
 
template<typename T, unsigned int N> 
inline tFixBlock<T,N> operator*(T a, const tFixBlock<T,N> &b){tFixBlock<T,N> R(a); return R *= b;}
template<typename T, unsigned int N> 
inline tFixBlock<T,N> operator+(T a, const tFixBlock<T,N> &b){tFixBlock<T,N> R(a); return R += b;}
template<typename T, unsigned int N> 
inline tFixBlock<T,N> operator-(T a, const tFixBlock<T,N> &b){tFixBlock<T,N> R(a); return R -= b;}
 
#endif
