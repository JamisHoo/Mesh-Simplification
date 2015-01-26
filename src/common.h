/******************************************************************************
 *  Copyright (c) 2014. All rights reserved.
 *
 *  Project: Mesh Simplification
 *  Filename: common.h 
 *  Version: 1.0
 *  Author: Jinming Hu
 *  E-mail: hjm211324@gmail.com
 *  Date: Jun. 24, 2014
 *  Time: 14:14:47
 *  Description: basic declaration, statement, public variables ...
 *****************************************************************************/
#ifndef COMMON_H
#define COMMON_H

#include "vector3.h"
#include "triple.h"
#include "matrix.h"
#include <vector>

#ifdef DEBUG
#include <iostream>
using std::cout;
using std::cerr;
using std::endl;
using std::flush;
#endif

namespace MeshSimplification {
    using INTTYPE = int;
    using REALTYPE = double;
    using Vector = Vector3<REALTYPE, INTTYPE>;
    using Matrix = SquareMatrix<REALTYPE, INTTYPE, 4>;

    // faces, normals
    using Result = std::pair< std::vector<REALTYPE>, std::vector<REALTYPE> >;

    class ObjParser;
    class QuadricMethod;
    class ImprovedQuadricMethod;

}

#endif /* COMMON_H */
