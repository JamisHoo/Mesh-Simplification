/******************************************************************************
 *  Copyright (c) 2014. All rights reserved.
 *
 *  Project: Mesh Simplification
 *  Filename: matrix.h 
 *  Version: 1.0
 *  Author: Jinming Hu
 *  E-mail: hjm211324@gmail.com
 *  Date: Jun. 24, 2014
 *  Time: 19:32:14
 *  Description: square matrix class, matrix addition
 *****************************************************************************/
#ifndef MATRIX_H
#define MATRIX_H

#include <cstring>

template < class REALTYPE, class INTTYPE, INTTYPE ORDER >
class SquareMatrix {
    REALTYPE elem[ORDER][ORDER];
    INTTYPE _order;
public:
    INTTYPE order() const { return _order; }
    
    SquareMatrix() {
        _order = ORDER;
        memset(elem, 0, sizeof(elem));
    }

    void clear() { memset(elem, 0, sizeof(elem)); }

    bool operator==(const SquareMatrix& m) const {
        for (INTTYPE i = 0; i < ORDER; ++i)
            for (INTTYPE j = 0; j < ORDER; ++j)
                if (elem[i][j] != m(i, j)) return false;
        return true;
    }

    REALTYPE& operator()(const INTTYPE i, const INTTYPE j) {
        return elem[i][j];
    }

    REALTYPE operator()(const INTTYPE i, const INTTYPE j) const {
        return elem[i][j];
    }

    SquareMatrix operator+(const SquareMatrix& m) const {
        SquareMatrix newMatrix;
        for (INTTYPE i = 0; i < ORDER; ++i)
            for (INTTYPE j = 0; j < ORDER; ++j)
                newMatrix(i, j) = elem[i][j] + m(i, j);
        return newMatrix;
    }

    SquareMatrix operator+=(const SquareMatrix& m) {
        for (INTTYPE i = 0; i < ORDER; ++i)
            for (INTTYPE j = 0; j < ORDER; ++j)
                elem[i][j] += m(i, j);
        return *this;
    }

    SquareMatrix operator*(const REALTYPE k) const {
        SquareMatrix newMatrix;
        for (INTTYPE i = 0; i < ORDER; ++i)
            for (INTTYPE j = 0; j < ORDER; ++j)
                newMatrix(i, j) = elem[i][j] * k;
        return newMatrix;
    }

    SquareMatrix operator*=(const REALTYPE k) {
        for (INTTYPE i = 0; i < ORDER; ++i)
            for (INTTYPE j = 0; j < ORDER; ++j)
                elem[i][j] *= k;
        return *this;
    }
};


#endif /* MATRIX_H */
