/******************************************************************************
 *  Copyright (c) 2014. All rights reserved.
 *
 *  Project: Mesh Simplification
 *  Filename: parser.h 
 *  Version: 1.0
 *  Author: Jinming Hu
 *  E-mail: hjm211324@gmail.com
 *  Date: Jun. 24, 2014
 *  Time: 14:13:53
 *  Description: parse obj file, only support v, f and #
 *****************************************************************************/
#ifndef PARSER_H
#define PARSER_H

#include "common.h"
#include <sstream>
#include <cstdio>
#include <vector>
#include <map>
#include <limits>

class MeshSimplification::ObjParser {
    std::string removeSpaces(const std::string& str) {
        // delete space characters
        std::string spaces = "\n\r\t\v\f ";
        auto s = str.find_first_not_of(spaces);
        auto e = str.find_last_not_of(spaces);
        if (s == std::string::npos || e == std::string::npos) return "";
        return str.substr(s, e - s + 1);
    }
public:
    template < class VERTEXCALLBACKFUNC, class FACECALLBACKFUNC >
    ObjParser(std::string& file, Vector& minPoint, Vector& maxPoint, 
              VERTEXCALLBACKFUNC vertexCallback, FACECALLBACKFUNC faceCallback) {
        std::stringstream cont(file);
        REALTYPE minx, miny, minz, maxx, maxy, maxz;
        minx = miny = minz = std::numeric_limits<REALTYPE>::max();
        maxx = maxy = maxz = std::numeric_limits<REALTYPE>::min();
        std::string line;
        while (getline(cont, line)) {
            line = removeSpaces(line);
            if (line.length() == 0) continue; // empty line
            
            if (line[0] == '#') continue; // comments
            if (line.length() >= 2 && line.substr(0, 2) == "v ") { // vertex
                REALTYPE x, y, z;
                if (sizeof(REALTYPE) == 4)
                    sscanf(line.c_str(), "v %f %f %f", &x, &y, &z);
                else if (sizeof(REALTYPE) == 8)
                    sscanf(line.c_str(), "v %lf %lf %lf", &x, &y, &z);
                else if (sizeof(REALTYPE) == 16)
                    sscanf(line.c_str(), "v %Lf %Lf %Lf", &x, &y, &z);
                // add other types if necessary
                else assert(0);
                
                vertexCallback(x, y, z);
                minx = std::min(x, minx);
                miny = std::min(y, miny);
                minz = std::min(z, minz);
                maxx = std::max(x, maxx);
                maxy = std::max(y, maxy);
                maxz = std::max(z, maxz);
            }
            else if (line.length() >= 2 && line.substr(0, 2) == "f ") { // face
                INTTYPE v0, v1, v2;
                sscanf(line.c_str(), "f %d %d %d", &v0, &v1, &v2);
                faceCallback(v0 - 1, v1 - 1, v2 - 1);
            }
            // otherwise, not supported for now, ignore
        }
        minPoint = Vector(minx, miny, minz);
        maxPoint = Vector(maxx, maxy, maxz);
    }
    
    void output(const std::vector<REALTYPE>& pts, std::ostream& out) {
        assert(pts.size() % 9 == 0);
        INTTYPE numFaces = pts.size() / 9;
        std::map<Vector, INTTYPE> points;
        for (INTTYPE i = 0; i < pts.size(); i += 3) 
            points.insert(std::make_pair(Vector(pts[i], pts[i + 1], pts[i + 2]), points.size()));
        INTTYPE numPoints = points.size();

        out << "# number of vertexes: " << numPoints
            << ", number of faces: " << numFaces << std::endl;
        
        INTTYPE n = 0;
        for (auto &pt: points) {
            out << "v " << pt.first[0] << ' ' 
                        << pt.first[1] << ' ' 
                        << pt.first[2] << std::endl;
            pt.second = n++;
        }

        for (INTTYPE i = 0; i < pts.size(); i += 9) 
            out << "f "
                << points[Vector(pts[i    ], pts[i + 1], pts[i + 2])] + 1 << ' '
                << points[Vector(pts[i + 3], pts[i + 4], pts[i + 5])] + 1 << ' '
                << points[Vector(pts[i + 6], pts[i + 7], pts[i + 8])] + 1 << std::endl;
    }
};

#endif /* PARSER_H */
