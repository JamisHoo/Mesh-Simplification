/******************************************************************************
 *  Copyright (c) 2014. All rights reserved.
 *
 *  Project: 
 *  Filename: improved_quadric.h 
 *  Version: 1.0
 *  Author: Jinming Hu
 *  E-mail: hjm211324@gmail.com
 *  Date: Jun. 28, 2014
 *  Time: 08:21:30
 *  Description: 
 *****************************************************************************/
#ifndef IMPROVED_QUADRIC_H
#define IMPROVED_QUADRIC_H

#include "common.h"
#include <vector>
#include <algorithm>
#include <cmath>

class MeshSimplification::ImprovedQuadricMethod {
    class Vertex;

    struct Face {
        INTTYPE vertexes[3];
        REALTYPE error[4];
        Vector normal;
        bool dirty;
        bool deleted;

        Face(): dirty(0), deleted(0), error{0, 0, 0, 0} {
        }

        Face(const INTTYPE v0, const INTTYPE v1, const INTTYPE v2): 
            vertexes{v0, v1, v2}, dirty(0), deleted(0), error{0, 0, 0, 0} {
        }
    
        
    };

    struct Vertex {
        Vector coord;

        Vertex(const REALTYPE x, const REALTYPE y, const REALTYPE z): coord(x, y, z) {
        }
        Vertex(): coord(0, 0, 0) { }

        REALTYPE operator[](const INTTYPE k) const {
            assert(k >= 0 && k < 3);
            return coord[k];
        }

        INTTYPE tstart;
        INTTYPE tcount;
        INTTYPE border;

        Matrix Q;
    };

    struct Ref {
        INTTYPE tid;
        INTTYPE tvertex;
    };

    std::vector<Face> faces;
    std::vector<Vertex> vertexes;
    std::vector<Ref> refs;

    bool insertSwitch;
    
    REALTYPE vertexError(const Matrix& q, 
                         const REALTYPE x, const REALTYPE y, const REALTYPE z) {
        REALTYPE err = 0;
        REALTYPE vec[4] = { x, y, z, 1 };
        for (INTTYPE i = 0; i < q.order(); ++i)
            for (INTTYPE j = 0; j < q.order(); ++j)
                err += vec[i] * q(i, j) * vec[j];
        return err;
    }

    REALTYPE calcError(const Vertex& v1, const Vertex& v2, Vector& result) {
        Matrix q = v1.Q + v2.Q;
        REALTYPE err = 0;
        REALTYPE det = + q(0, 0) * q(1, 1) * q(2, 2) - q(0, 2) * q(1, 1) * q(2, 0)
                       + q(1, 0) * q(2, 1) * q(0, 2) - q(0, 1) * q(1, 0) * q(2, 2)
                       + q(2, 0) * q(0, 1) * q(1, 2) - q(0, 0) * q(1, 2) * q(2, 1);
        
        if (det) {
            result[0] = REALTYPE(-1) / det * (
                        + q(0, 1) * q(1, 2) * q(2, 3) - q(0, 3) * q(1, 2) * q(2, 1)
                        + q(0, 2) * q(1, 3) * q(2, 1) - q(0, 2) * q(1, 1) * q(2, 3)
                        + q(0, 3) * q(1, 1) * q(2, 2) - q(0, 1) * q(1, 3) * q(2, 2));
            result[1] = REALTYPE(1) / det * (
                        + q(0, 0) * q(1, 2) * q(2, 3) - q(0, 3) * q(1, 2) * q(2, 0)
                        + q(1, 0) * q(2, 2) * q(0, 3) - q(0, 2) * q(1, 0) * q(2, 3)
                        + q(2, 0) * q(0, 2) * q(1, 3) - q(0, 0) * q(1, 3) * q(2, 2));
            result[2] = REALTYPE(-1) / det * (
                        + q(0, 0) * q(1, 1) * q(2, 3) - q(0, 3) * q(1, 1) * q(2, 0)
                        + q(1, 0) * q(2, 1) * q(0, 3) - q(1, 3) * q(2, 1) * q(0, 0)
                        + q(2, 0) * q(0, 1) * q(1, 3) - q(2, 3) * q(0, 1) * q(1, 0));
            err = vertexError(q, result[0], result[1], result[2]);
        }
        else {
            Vector u = v1.coord;
            Vector v = v2.coord;
            Vector v12 = (u + v) * 0.5;
            REALTYPE err1 = vertexError(q, u[0], u[1], u[2]);
            REALTYPE err2 = vertexError(q, v[0], v[1], v[2]);
            REALTYPE err3 = vertexError(q, v12[0], v12[1], v12[2]);
            err = std::min(err1, std::min(err2, err3));
            if (err1 == err) result = u;
            if (err2 == err) result = v;
            if (err3 == err) result = v12;
        }
        return err;
    }

    bool flipped(const Vector& pos, const INTTYPE index0, const INTTYPE index1, 
                 const Vertex& v0, const Vertex& v1, std::vector<INTTYPE>& deleted) {
        INTTYPE borderCount = 0;
        for (INTTYPE i = 0; i < v0.tcount; ++i) {
            const Face& f = faces[refs[v0.tstart + i].tid];
            if (f.deleted) continue;

            INTTYPE s = refs[v0.tstart + i].tvertex;
            INTTYPE id1 = f.vertexes[(s + 1) % 3];
            INTTYPE id2 = f.vertexes[(s + 2) % 3];

            if (id1 == index1 || id2 == index1) {
                ++borderCount;
                deleted[i] = 1;
                continue;
            }
            Vector d1 = (vertexes[id1].coord - pos).normalize();
            Vector d2 = (vertexes[id2].coord - pos).normalize();

            if (std::abs(d1 * d2 > 0.9999999)) return true;

            Vector n = crossProduct(d1, d2).normalize();

            deleted[i] = 0;
            if (n * f.normal < 0.2) return true;
        }
        return false;
    }

    void updateFaces(const INTTYPE index, const Vertex& v, 
                     const std::vector<INTTYPE>& deleted, INTTYPE& numDeleted) {
        Vector p;
        for (INTTYPE i = 0; i < v.tcount; ++i) {
            Ref& r = refs[v.tstart + i];
            Face& f = faces[r.tid];
            if (f.deleted) continue;
            if (deleted[i]) {
                f.deleted = 1;
                ++numDeleted;
                continue;
            }
            f.vertexes[r.tvertex] = index;
            f.dirty = 1;
            f.error[0] = calcError(vertexes[f.vertexes[0]], vertexes[f.vertexes[1]], p);
            f.error[1] = calcError(vertexes[f.vertexes[1]], vertexes[f.vertexes[2]], p);
            f.error[2] = calcError(vertexes[f.vertexes[2]], vertexes[f.vertexes[0]], p);
            f.error[3] = std::min(f.error[0], std::min(f.error[1], f.error[2]));
            refs.push_back(r);
        }
    }

    void updateMesh(const INTTYPE iteration) {
        if (iteration > 0) {
            INTTYPE dst = 0;
            for (INTTYPE i = 0; i < faces.size(); ++i) 
                if (!faces[i].deleted)
                    faces[dst++] = faces[i];
            faces.resize(dst);
        }

        if (iteration == 0) {
            for (auto &v: vertexes) v.Q.clear();
            for (auto &f: faces) {
                Vector n, p[3];
                for (INTTYPE j = 0; j < 3; ++j) 
                    p[j] = vertexes[f.vertexes[j]].coord;
                n = crossProduct(p[1] - p[0], p[2] - p[0]).normalize();
                f.normal = n;
                
                Matrix k;
                REALTYPE pp[4] = { n[0], n[1], n[2], -1.0 * n * p[0]};
                for (INTTYPE r = 0; r < k.order(); ++r)
                    for (INTTYPE c = 0; c < k.order(); ++c)
                        k(r, c) = pp[r] * pp[c];
                for (INTTYPE j = 0; j < 3; ++j)
                    vertexes[f.vertexes[j]].Q += k;
            }
            for (auto &f: faces) {
                Vector p;
                for (INTTYPE j = 0; j < 3; ++j)
                    f.error[j] = calcError(vertexes[f.vertexes[j]], 
                                           vertexes[f.vertexes[(j + 1) % 3]], p);
        
                f.error[3] = *std::min_element(f.error, f.error + 3);
            }
        }
        

        for (auto &v: vertexes) 
            v.tcount = v.tstart = 0;

        for (auto const &f: faces) 
            for (INTTYPE j = 0; j < 3; ++j) 
                vertexes[f.vertexes[j]].tcount++;
        INTTYPE tstart = 0;
        for (auto &v: vertexes) {
            v.tstart = tstart;
            tstart += v.tcount;
            v.tcount = 0;
        }
        refs.resize(faces.size() * 3);
        for (INTTYPE i = 0; i < faces.size(); ++i) 
            for (INTTYPE j = 0; j < 3; ++j) {
                Vertex& v = vertexes[faces[i].vertexes[j]];
                refs[v.tstart + v.tcount].tid = i;
                refs[v.tstart + v.tcount].tvertex = j;
                ++v.tcount;
            }

        if (iteration == 0) { 
            std::vector<INTTYPE> vcount, vids;
            for (auto &v: vertexes) v.border = 0;

            for (INTTYPE i = 0; i < vertexes.size(); ++i) {
                Vertex& v = vertexes[i];
                vcount.clear();
                vids.clear();
                for (INTTYPE j = 0; j < v.tcount; ++j) {
                    INTTYPE k = refs[v.tstart + j].tid;
                    Face& f = faces[k];
                    for (INTTYPE l = 0; l < 3; ++l) {
                        INTTYPE ofs = 0, id = f.vertexes[k];
                        while (ofs < vcount.size()) {
                            if (vids[ofs] == id) break;
                            ++ofs;
                        }
                        if (ofs == vcount.size()) {
                            vcount.push_back(1);
                            vids.push_back(id);
                        }
                        else ++vcount[ofs];
                    }
                }
                for (INTTYPE j = 0; j < vcount.size(); ++j) 
                    if (vcount[j] == 1) vertexes[vids[j]].border = 1;
            }
        }
    }

    void compactMesh() {
        INTTYPE dst = 0;
        for (auto &v:vertexes) v.tcount = 0;
        for (INTTYPE i = 0; i < faces.size(); ++i)
            if (!faces[i].deleted) {
                faces[dst++] = faces[i];
                for (INTTYPE j = 0; j < 3; ++j)
                    vertexes[faces[i].vertexes[j]].tcount = 1;
            }
        faces.resize(dst);
        dst = 0;
        for (INTTYPE i = 0; i < vertexes.size(); ++i)
            if (vertexes[i].tcount) {
                vertexes[i].tstart = dst;
                vertexes[dst].coord = vertexes[i].coord;
                ++dst;
            }
        for (auto &f: faces)
            for (INTTYPE j = 0; j < 3; ++j)
                f.vertexes[j] = vertexes[f.vertexes[j]].tstart;
        vertexes.resize(dst);
    }


public:
    ImprovedQuadricMethod(): insertSwitch(1) { }
    
    void switchOff() {
        insertSwitch = 0;
    }
    
    void simplify(const std::vector<REALTYPE> samplingThresholds,
                  std::vector<Result>& results) {
        assert(insertSwitch == 0);
        if (!samplingThresholds.size()) return;
        INTTYPE numDeleted = 0;
        std::vector<INTTYPE> deleted0, deleted1;
        INTTYPE numFaces = faces.size();

        auto th = samplingThresholds.begin();
        for (INTTYPE iteration = 0; ; ++iteration) {
            
            if (iteration % 5 == 0) updateMesh(iteration);

            if (REALTYPE(numFaces - numDeleted) / numFaces <= *th) {
                sample(results);
                if (++th == samplingThresholds.end()) break;
            }

            
            for (auto &f: faces) f.dirty = 0;

            REALTYPE agressiveness = 4;
            REALTYPE threshold = 1e-9 * pow(REALTYPE(iteration + 3), agressiveness);

            for (auto &f: faces) {
                if (f.error[3] > threshold) continue;
                if (f.deleted) continue;
                if (f.dirty) continue;

                
                for (INTTYPE i = 0; i < 3; ++i) 
                    if (f.error[i] < threshold) {
                        INTTYPE index0 = f.vertexes[i];
                        Vertex& v0 = vertexes[index0];
                        INTTYPE index1 = f.vertexes[(i + 1) % 3];
                        Vertex& v1 = vertexes[index1];

                        if (v0.border != v1.border) continue;

                        Vector pos;
                        calcError(v0, v1, pos);
                        static int k = 0;
                        deleted0.resize(v0.tcount);
                        deleted1.resize(v1.tcount);

                        if (flipped(pos, index0, index1, v0, v1, deleted0)) continue;
                        if (flipped(pos, index1, index0, v1, v0, deleted1)) continue;
                        
                        v0.coord = pos;
                        v0.Q += v1.Q;
                        INTTYPE tstart = refs.size();

                        updateFaces(index0, v0, deleted0, numDeleted);
                        updateFaces(index0, v1, deleted1, numDeleted);

                        INTTYPE tcount = refs.size() - tstart;

                        if (tcount <= v0.tcount)
                            //memcpy(&refs[v0.tstart], &refs[tstart], tcount * sizeof(Ref));
                            std::copy(refs.begin() + tstart, refs.begin() + tstart + tcount, refs.begin() + v0.tstart);
                        else
                            v0.tstart = tstart;

                        v0.tcount = tcount;
                        break;
                    }
                if (REALTYPE(numFaces - numDeleted) / numFaces <= *th) break;
            }
        }

        compactMesh();
    }

    void insertVertex(const REALTYPE x, const REALTYPE y, const REALTYPE z) {
        assert(insertSwitch);
        assert(!isnan(x));
        assert(!isnan(y));
        assert(!isnan(z));
        assert(!isinf(x));
        assert(!isinf(y));
        assert(!isinf(z));
        vertexes.push_back(Vertex(x, y, z));
    }

    void insertFace(const INTTYPE v0, const INTTYPE v1, const INTTYPE v2) {
        assert(insertSwitch);
        faces.push_back(Face(v0, v1, v2));
   }
    
    void sample(std::vector<Result>& results) const {
        assert(!insertSwitch);
        std::vector<REALTYPE> fs;
        std::vector<REALTYPE> fns;

        for (auto const &f: faces)
            if (!f.deleted) {
                fs.push_back(vertexes[f.vertexes[0]].coord[0]);
                fs.push_back(vertexes[f.vertexes[0]].coord[1]);
                fs.push_back(vertexes[f.vertexes[0]].coord[2]);
                fs.push_back(vertexes[f.vertexes[1]].coord[0]);
                fs.push_back(vertexes[f.vertexes[1]].coord[1]);
                fs.push_back(vertexes[f.vertexes[1]].coord[2]);
                fs.push_back(vertexes[f.vertexes[2]].coord[0]);
                fs.push_back(vertexes[f.vertexes[2]].coord[1]);
                fs.push_back(vertexes[f.vertexes[2]].coord[2]);

                fns.push_back(f.normal[0]);
                fns.push_back(f.normal[1]);
                fns.push_back(f.normal[2]);
                fns.push_back(f.normal[0]);
                fns.push_back(f.normal[1]);
                fns.push_back(f.normal[2]);   
                fns.push_back(f.normal[0]);
                fns.push_back(f.normal[1]);
                fns.push_back(f.normal[2]);   
            }
        results.push_back(std::make_pair(fs, fns));
    }
    
};


#endif /* IMPROVED_QUADRIC_H */
