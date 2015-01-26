/******************************************************************************
 *  Copyright (c) 2014. All rights reserved.
 *
 *  Project: Mesh Simplification
 *  Filename: triple.h 
 *  Version: 1.0
 *  Author: Jinming Hu
 *  E-mail: hjm211324@gmail.com
 *  Date: Jun. 24, 2014
 *  Time: 20:47:09
 *  Description: triple class, just like pair
 *****************************************************************************/
#ifndef TRIPLE_H
#define TRIPLE_H

template < class T1, class T2, class T3 >
struct Triple {
    T1 first;
    T2 second;
    T3 third;

    Triple(): first(T1()), second(T2()), third(T3()) { }
    Triple(const T1& t1, const T2& t2, const T3& t3): 
        first(t1), second(t2), third(t3) { }
};

template < class T1, class T2, class T3 >
inline bool operator==(const Triple<T1, T2, T3>& x, const Triple<T1, T2, T3>& y) {
    return x.first == y.first && x.second == y.second && x.third == y.third;
}

template < class T1, class T2, class T3 >
inline bool operator<(const Triple<T1, T2, T3>& x, const Triple<T1, T2, T3>& y) {
    return x.first < y.first ||
           (!(y.first < x.first) && x.second < y.second) ||
           (!(y.first < x.first) && !(y.second < x.second) && x.third < y.third);
}

template < class T1, class T2, class T3 >
inline bool operator!=(const Triple<T1, T2, T3>& x, const Triple<T1, T2, T3>& y) {
    return !(x == y);
}

template < class T1, class T2, class T3 >
inline bool operator>(const Triple<T1, T2, T3>& x, const Triple<T1, T2, T3>& y) {
    return y < x;
}

template < class T1, class T2, class T3 >
inline bool operator<=(const Triple<T1, T2, T3>& x, const Triple<T1, T2, T3>& y) {
    return !(y < x);
}


template < class T1, class T2, class T3 >
inline bool operator>=(const Triple<T1, T2, T3>& x, const Triple<T1, T2, T3>& y) {
    return !(x < y);
}

template < class T1, class T2, class T3 >
inline Triple<T1, T2, T3> make_triple(const T1& t1, const T2& t2, const T3& t3) {
    return Triple<T1, T2, T3>(t1, t2, t3);
}

#endif /* TRIPLE_H */
