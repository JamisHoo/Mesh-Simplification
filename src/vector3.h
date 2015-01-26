/******************************************************************************
 *  Copyright (c) 2014. All rights reserved.
 *
 *  Project: Mesh Simplification
 *  Filename: vector3.h 
 *  Version: 1.0
 *  Author: Jinming Hu
 *  E-mail: hjm211324@gmail.com
 *  Date: Apr. 07, 2014
 *  Time: 14:50:31
 *  Description: Vector3 is a three-dimensional vector, 
 *               providing add, subtract, multiply functions.
 *****************************************************************************/
#ifndef VECTOR3_H
#define VECTOR3_H

#include <cassert>
#include <cstring>
#include <cmath>
#include <iostream>

template <class ELEMTYPE, class COUNTTYPE>
class Vector3 {
    ELEMTYPE elem[3];

public:
    Vector3(const ELEMTYPE a, const ELEMTYPE b, const ELEMTYPE c) {
        elem[0] = a, elem[1] = b, elem[2] = c;
    }

    Vector3() { memset(elem, 0, sizeof(ELEMTYPE) * 3); }

    ELEMTYPE norm() const {
        return sqrt(*this * *this);
    }

    Vector3 normalize() const {
        assert(norm() != 0);
        return *this * (1 / norm());
    }

    Vector3 operator+(const Vector3& v) const {
        return Vector3(elem[0] + v.elem[0],
                       elem[1] + v.elem[1],
                       elem[2] + v.elem[2]);
    }

    Vector3 operator+=(const Vector3& v) {
        return Vector3(elem[0] += v.elem[0],
                       elem[1] += v.elem[1],
                       elem[2] += v.elem[2]);
    }

    Vector3 operator-(const Vector3& v) const {
        return Vector3(elem[0] - v.elem[0],
                       elem[1] - v.elem[1],
                       elem[2] - v.elem[2]);
    }

    Vector3 operator-=(const Vector3& v) {
        return Vector3(elem[0] -= v.elem[0],
                       elem[1] -= v.elem[1],
                       elem[2] -= v.elem[2]);
    }

    ELEMTYPE operator*(const Vector3& v) const {
        return innerProduct(*this, v);
    }

    Vector3 operator*(const ELEMTYPE t) const {
        return Vector3(elem[0] * t, elem[1] * t, elem[2] * t);
    }

    Vector3 operator*=(const ELEMTYPE t) {
        return Vector3(elem[0] *= t,
                       elem[1] *= t,
                       elem[2] *= t);
    }

    ELEMTYPE& operator[](const COUNTTYPE k) {
        assert(k < 3);
        return elem[k];
    }

    ELEMTYPE operator[](const COUNTTYPE k) const {
        assert(k < 3);
        return elem[k];
    }
    
    bool operator==(const Vector3& v) const {
        return elem[0] == v.elem[0] && 
               elem[1] == v.elem[1] &&
               elem[2] == v.elem[2];
    }
    
};

#ifdef DEBUG
template < class ELEMTYPE, class COUNTTYPE >
std::ostream& operator<<(std::ostream& out, const Vector3<ELEMTYPE, COUNTTYPE>& v) {
    out << '(' << v[0] << ", " << v[1] << ", " << v[2] << ')';
    return out;
}
#endif

template < class ELEMTYPE, class COUNTTYPE >
inline Vector3<ELEMTYPE, COUNTTYPE> operator*(const ELEMTYPE t, const Vector3<ELEMTYPE, COUNTTYPE>& v) {
    return Vector3<ELEMTYPE, COUNTTYPE>(v[0] * t, v[1] * t, v[2] * t);
}

template < class ELEMTYPE, class COUNTTYPE >
inline ELEMTYPE innerProduct(const Vector3<ELEMTYPE, COUNTTYPE>& v0, const Vector3<ELEMTYPE, COUNTTYPE>& v1) {
    return v0[0] * v1[0] +
           v0[1] * v1[1] +
           v0[2] * v1[2];
}

template < class ELEMTYPE, class COUNTTYPE >
inline Vector3<ELEMTYPE, COUNTTYPE> crossProduct(const Vector3<ELEMTYPE, COUNTTYPE>& v0, const Vector3<ELEMTYPE, COUNTTYPE>& v1) {
    return Vector3<ELEMTYPE, COUNTTYPE>(v0[1] * v1[2] - v1[1] * v0[2],
                                        v0[2] * v1[0] - v1[2] * v0[0],
                                        v0[0] * v1[1] - v1[0] * v0[1]);
}

template < class ELEMTYPE, class COUNTTYPE >
inline bool operator<(const Vector3<ELEMTYPE, COUNTTYPE>& v1, const Vector3<ELEMTYPE, COUNTTYPE>& v2) {
    return v1[0] < v2[0] || 
           (!(v2[0] < v1[0]) && v1[1] < v2[1]) ||
           (!(v2[0] < v1[0]) && !(v2[1] < v1[1]) && v1[2] < v2[2]);
}


#endif /* VECTOR3_H */
