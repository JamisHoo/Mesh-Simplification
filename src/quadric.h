/******************************************************************************
 *  Copyright (c) 2014. All rights reserved.
 *
 *  Project: Mesh Simplification
 *  Filename: quadric.h 
 *  Version: 1.0
 *  Author: Jinming Hu
 *  E-mail: hjm211324@gmail.com
 *  Date: Jun. 24, 2014
 *  Time: 14:52:21
 *  Description: input vertexes and faces, simplify the mesh with quadric method
 *****************************************************************************/
#ifndef QUADRIC_H
#define QUADRIC_H

#include "common.h"
#include <vector>
#include <set>
#include <map>
#include <algorithm>

class MeshSimplification::QuadricMethod {
    class Vertex;

    struct Face {
        // don't change vertexes other than using replace
        Vertex* vertexes[3];
        Vector normal;
        bool dirty;

        Face(Vertex* v0, Vertex* v1, Vertex* v2): dirty(0) {
            vertexes[0] = v0;
            vertexes[1] = v1;
            vertexes[2] = v2;
            calcNormal();
        }

        Vertex* operator[](const INTTYPE k) const {
            assert(k == 0 || k == 1 || k == 2);
            return vertexes[k];
        }
        
        void calcNormal() {
            normal = crossProduct(vertexes[1] -> coord - vertexes[0] -> coord, 
                                  vertexes[2] -> coord - vertexes[0] -> coord).normalize();
        }
        
        INTTYPE replace(Vertex* o0, Vertex* o1, Vertex* n) {
            INTTYPE numReplaces = 0;
            for (INTTYPE i = 0; i < 3; ++i) {
                if (vertexes[i] == o0) 
                    vertexes[i] = n, ++numReplaces;
                if (vertexes[i] == o1) 
                    vertexes[i] = n, ++numReplaces;
            }
            assert(numReplaces == 1 || numReplaces == 2);
            if (numReplaces == 1) calcNormal();
            return numReplaces;
        }
    };

    struct Vertex {
        Vector coord;
        bool dirty;
        std::set<Vertex*> adjacentVertexes;
        std::set<Face*> adjacentFaces;

        Matrix Q;

        std::set< Triple<REALTYPE, Vertex*, Vertex*> >::iterator costIte;
        bool iteValid;


        Vector contractVec;

        Vertex(const REALTYPE x, const REALTYPE y, const REALTYPE z): 
            coord(x, y, z), dirty(0), iteValid(0) { }
        Vertex(const Vector& v): coord(v[0], v[1], v[2]), dirty(0), iteValid(0) { }
        
        REALTYPE& operator[](const int k) { return coord[k]; }
        REALTYPE operator[](const int k) const { return coord[k]; }

        void replace(Vertex* o, Vertex* n) {
            INTTYPE numReplaces = adjacentVertexes.erase(o);
            assert(numReplaces == 1);
            adjacentVertexes.insert(n);
        }

        void calcMatrixQ() {
            Q.clear();
            for (auto f: adjacentFaces) {
                Vector fn = f -> normal;
                REALTYPE d = coord * fn * -1;
                Matrix q;
                for (INTTYPE i = 0; i < q.order(); ++i)
                    for (INTTYPE j = 0; j < q.order(); ++j)
                        q(i, j) = (i == 3? d: fn[i]) * (j == 3? d: fn[j]);
                Q += q;
            }
        }

        Vector minimizeDeltaV(const Matrix& m, const Vector& u, const Vector& v) {
            REALTYPE deno = m(0, 0) * m(1, 1) * m(2, 2) - 
                            m(0, 0) * m(1, 2) * m(2, 1) -
                            m(0, 1) * m(1, 0) * m(2, 2) +
                            m(0, 1) * m(1, 2) * m(2, 0) + 
                            m(0, 2) * m(1, 0) * m(2, 1) -
                            m(0, 2) * m(1, 1) * m(2, 0);
            if (deno <= 1e-3) return (u + v) * REALTYPE(0.5);
            REALTYPE x = m(0, 1) * m(1, 3) * m(2, 2) - 
                         m(0, 1) * m(1, 2) * m(2, 3) +
                         m(0, 2) * m(1, 1) * m(2, 3) -
                         m(0, 2) * m(1, 3) * m(2, 1) -
                         m(0, 3) * m(1, 1) * m(2, 2) +
                         m(0, 3) * m(1, 2) * m(2, 1);
            REALTYPE y = m(0, 0) * m(1, 2) * m(2, 3) -
                         m(0, 0) * m(1, 3) * m(2, 2) -
                         m(0, 2) * m(1, 0) * m(2, 3) +
                         m(0, 2) * m(1, 3) * m(2, 0) +
                         m(0, 3) * m(1, 0) * m(2, 2) -
                         m(0, 3) * m(1, 2) * m(2, 0);
            REALTYPE z = m(0, 0) * m(1, 3) * m(2, 1) -
                         m(0, 0) * m(1, 1) * m(2, 3) +
                         m(0, 1) * m(1, 0) * m(2, 3) -
                         m(0, 1) * m(1, 3) * m(2, 0) -
                         m(0, 3) * m(1, 0) * m(2, 1) +
                         m(0, 3) * m(1, 1) * m(2, 0);
            return Vector(x / deno, y / deno, z / deno);
        }


        Triple<REALTYPE, Vertex*, Vertex*> findLeastCost() {
            REALTYPE leastCost;
            Vertex* bestAdj = nullptr;
            // suppose this is not a isolated vertex
            assert(adjacentVertexes.size());
            Vector v;
            for (auto adj: adjacentVertexes) {
                Matrix q = Q + adj -> Q;
                v = (coord + adj -> coord) * REALTYPE(0.5);
                //v = minimizeDeltaV(q, coord, adj -> coord);
                REALTYPE deltaV = 0;
                for (INTTYPE i = 0; i < q.order(); ++i) 
                    for (INTTYPE j = 0; j < q.order(); ++j)
                        deltaV += (i == 3? REALTYPE(1): v[i]) * q(i, j) * (j == 3? REALTYPE(1): v[j]);
                if (!bestAdj || deltaV < leastCost)
                    leastCost = deltaV, bestAdj = adj, contractVec = v;
            }
            // as long as this is not a isolated vertex
            assert(bestAdj);
            
            return make_triple(leastCost, this, bestAdj);
        }

    };

    void removeEdge(Vertex* u, Vertex* v, Vertex* n) {
        // update their adjacent faces
        for (auto adjf: u -> adjacentFaces) {
            INTTYPE numR = adjf -> replace(u, v, n);
            // face that is reserved
            // adjacent faces of new vertex
            if (numR == 1) n -> adjacentFaces.insert(adjf);
            // face that is removed
            else if (numR == 2) {
                for (INTTYPE i = 0; i < 3; ++i) 
                    if ((*adjf)[i] != n) 
                        assert((*adjf)[i] -> adjacentFaces.erase(adjf) == 1);
                // delete face from faces set
                adjf -> dirty = 1;
                --numFaces;
            }
        }

        for (auto adjf : v -> adjacentFaces) {
            if (adjf -> dirty) continue;
            INTTYPE numR = adjf -> replace(u, v, n);
            if (numR == 1) n -> adjacentFaces.insert(adjf);
            else if (numR == 2) {
                for (INTTYPE i = 0; i < 3; ++i)
                    if ((*adjf)[i] != n)
                        assert((*adjf)[i] -> adjacentFaces.erase(adjf) == 1);
                adjf -> dirty = 1;
                --numFaces;
            }
        }

        // update their adjacent vertexes
        for (auto adjv: u -> adjacentVertexes) {
            if (adjv == v) continue;
            // replace u or v with new vertex
            adjv -> replace(u, n);
            // remove them from costs set
            if (adjv -> iteValid) {
                costs.erase(adjv -> costIte);
                adjv -> iteValid = 0;
            }
            // adjacent vertexes of new vertex
            n -> adjacentVertexes.insert(adjv);
        }
        for (auto adjv: v -> adjacentVertexes) {
            if (adjv == u) continue;
            adjv -> replace(v, n);
            if (adjv -> iteValid) {
                costs.erase(adjv -> costIte);
                adjv -> iteValid = 0;
            }
            n -> adjacentVertexes.insert(adjv);
        }
        for (auto adjv: u -> adjacentVertexes) {
            if (adjv == v) continue;
            // update matrix Q
            adjv -> calcMatrixQ();
        }
        for (auto adjv: v -> adjacentVertexes) {
            if (adjv == u) continue;
            adjv -> calcMatrixQ();
        }

        // update costs set
        for (auto adjv: n -> adjacentVertexes) {
            auto p = costs.insert(adjv -> findLeastCost());
            assert(p.second);
            adjv -> costIte = p.first;
            adjv -> iteValid = 1;
        }

        // calculate cost of new vertex and add it to costs set
        n -> calcMatrixQ();
        auto p = costs.insert(n -> findLeastCost());
        assert(p.second);
        n -> costIte = p.first;
        n -> iteValid = 1;

        vertexes.insert(n);
        // delete old vertexes
        u -> dirty = v -> dirty = 1;
    }

    std::set<Vertex*> vertexes;
    std::set<Face*> faces;
    // only used when inserting
    std::vector<Vertex*> tmpVertexes;
    std::vector<Face*> tmpFaces;

    std::set< Triple<REALTYPE, Vertex*, Vertex*> > costs;
    bool insertSwitch;
    INTTYPE numVertexes;
    INTTYPE numFaces;

public:
    void simplify(const std::vector<REALTYPE>& samplingThresholds,
                  std::vector< Result >& results) {
        assert(insertSwitch == 0);
        if (!samplingThresholds.size()) return;
        costs.clear();

        // calculate matrix Q for each vertex
        for (auto v: vertexes) v -> calcMatrixQ();

        // calculate costs of contraction and make a set
        for (auto v: vertexes) {
            auto p = costs.insert(v -> findLeastCost());
            assert(p.second);
            v -> costIte = p.first;
            v -> iteValid = 1;
        }
        
        INTTYPE oldNumFaces = numFaces;
        auto th = samplingThresholds.begin();
        // continuously delete vertexes and faces
        while (numFaces) {
            if (numFaces <= oldNumFaces * *th) {
                // sample
                sample(results);
                if (++th == samplingThresholds.end()) break;
            }
            assert(costs.size());
            // pick the first element in the set, which has the least cost
            // (cost, vertex u, vertex v)
            auto edge = *costs.begin();
            costs.erase(costs.begin());
            
            if (edge.second -> dirty || edge.third -> dirty) continue;

            // construct a new vertex
            // coordinate of new vertex
            Vertex* newV = new Vertex(edge.second -> contractVec);
            
            // remove edge
            removeEdge(edge.second, edge.third, newV);
        }
    }

    QuadricMethod(): insertSwitch(1) { }

    ~QuadricMethod() {
        for (auto v: vertexes) delete v;
        for (auto f: faces) delete f;
    }

    INTTYPE numberOfVertexes() const { return numVertexes; }
    INTTYPE numOfFaces() const { return numFaces; }

    void switchOff() {
        insertSwitch = 0;
        for (auto v: tmpVertexes) vertexes.insert(v);
        for (auto f: tmpFaces) faces.insert(f);
        numVertexes = vertexes.size();
        numFaces = faces.size();
        tmpVertexes.clear();
        tmpFaces.clear();
    }

    void insertVertex(const REALTYPE x, const REALTYPE y, const REALTYPE z) {
        assert(insertSwitch);
        assert(!isnan(x));
        assert(!isnan(y));
        assert(!isnan(z));
        assert(!isinf(x));
        assert(!isinf(y));
        assert(!isinf(z));
        tmpVertexes.push_back(new Vertex(x, y, z));
    }

    void insertFace(const INTTYPE v0, const INTTYPE v1, const INTTYPE v2) {
        assert(insertSwitch);
        assert(v0 < tmpVertexes.size() && v1 < tmpVertexes.size() && v2 < tmpVertexes.size());
        assert(v0 >= 0 && v1 >= 0 && v2 >= 0);

        tmpFaces.push_back(new Face(tmpVertexes[v0], tmpVertexes[v1], tmpVertexes[v2]));
        
        tmpVertexes[v0] -> adjacentFaces.insert(tmpFaces.back());
        tmpVertexes[v1] -> adjacentFaces.insert(tmpFaces.back());
        tmpVertexes[v2] -> adjacentFaces.insert(tmpFaces.back());

        tmpVertexes[v0] -> adjacentVertexes.insert(tmpVertexes[v1]);
        tmpVertexes[v0] -> adjacentVertexes.insert(tmpVertexes[v2]);
        tmpVertexes[v1] -> adjacentVertexes.insert(tmpVertexes[v0]);
        tmpVertexes[v1] -> adjacentVertexes.insert(tmpVertexes[v2]);
        tmpVertexes[v2] -> adjacentVertexes.insert(tmpVertexes[v0]);
        tmpVertexes[v2] -> adjacentVertexes.insert(tmpVertexes[v1]);
    }

    void sample(std::vector<Result>& results) const {
        std::vector<REALTYPE> fs;
        std::vector<REALTYPE> fns;
        fs.reserve(numFaces * 9);
        fns.reserve(numFaces * 3);

        for (auto f: faces) {
            if (f -> dirty) continue;
            fs.push_back((*(*f)[0])[0]);
            fs.push_back((*(*f)[0])[1]);
            fs.push_back((*(*f)[0])[2]);
            fs.push_back((*(*f)[1])[0]);
            fs.push_back((*(*f)[1])[1]);
            fs.push_back((*(*f)[1])[2]);
            fs.push_back((*(*f)[2])[0]);
            fs.push_back((*(*f)[2])[1]);
            fs.push_back((*(*f)[2])[2]);

            fns.push_back(f -> normal[0]);
            fns.push_back(f -> normal[1]);
            fns.push_back(f -> normal[2]);
            fns.push_back(f -> normal[0]);
            fns.push_back(f -> normal[1]);
            fns.push_back(f -> normal[2]);
            fns.push_back(f -> normal[0]);
            fns.push_back(f -> normal[1]);
            fns.push_back(f -> normal[2]);
        }
        results.push_back(std::make_pair(fs, fns));

    }
};

#endif /* QUADRIC_H */
