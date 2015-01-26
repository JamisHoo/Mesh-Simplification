/******************************************************************************
 *  Copyright (c) 2014. All rights reserved.
 *
 *  Project: Mesh Simplification
 *  Filename: main.cc 
 *  Version: 1.0
 *  Author: Jinming Hu
 *  E-mail: hjm211324@gmail.com
 *  Date: Jun. 24, 2014
 *  Time: 14:39:01
 *  Description: 
 *****************************************************************************/

#include "common.h"
#include "parser.h"
#include "quadric.h"
#include "improved_quadric.h"
#include "preview.h"
#include <iostream>
#include <fstream>
#include <iterator>
#include <thread>
#include <chrono>

int main(int argc, char** argv) {
    using namespace std;
    using namespace Preview;
    using namespace MeshSimplification;
    

    // argc == 2, GUI mode
    // argc == 4, CLI mode
    assert(argc == 2 || argc == 4);

    // input obj file
    string file;
    ifstream fin(argv[1]);
    file.clear();
    char buff;
    while (fin.get(buff)) file += buff;
    fin.close();
    
    if (argc == 2) {
        
        QuadricMethod quadric;
        ImprovedQuadricMethod improvedQuadric;

        auto vertexCallback = [&quadric, &improvedQuadric](const REALTYPE x, const REALTYPE y, const REALTYPE z) { 
            quadric.insertVertex(x, y, z); 
            improvedQuadric.insertVertex(x, y, z);
        };
        auto faceCallback = [&quadric, &improvedQuadric](const int v0, const int v1, const int v2) { 
            quadric.insertFace(v0, v1, v2); 
            improvedQuadric.insertFace(v0, v1, v2);
            ++totalNumFaces; 
        };

        // parse obj file
        ObjParser objParser(file, lowerBound, upperBound, vertexCallback, faceCallback);

        // calculate view point position
        objCenter = (lowerBound + upperBound) * 0.5;
        REALTYPE viewPointDistance = std::max(std::max(upperBound[0] - lowerBound[0], 
                                                       upperBound[1] - lowerBound[1]),
                                              upperBound[2] - lowerBound[2]) * 2.3;
        viewPoint = objCenter + Vector(viewPointDistance, 0, 0);

        // cannot insert any vertex or face after this
        quadric.switchOff();
        improvedQuadric.switchOff();

        // sampling points
        // must be in sorted from largest to smallest
        // linear
        //for (INTTYPE i = numSamples; i > 0; --i) samplingThresholds.push_back(REALTYPE(1) / numSamples * i);
        // exponential
        for (INTTYPE i = 0; i < numSamples; ++i) samplingThresholds.push_back(exp(-0.06 * i));

        // simplify the mesh with multithreads
        auto simplifyFunc0 = [&] { quadric.simplify(samplingThresholds, results[0]); };
        auto simplifyFunc1 = [&] { improvedQuadric.simplify(samplingThresholds, results[1]); };
        thread simplifyThread0(simplifyFunc0);
        thread simplifyThread1(simplifyFunc1);

        // show result
        thread previewThread(initGL, argc, argv);
        previewThread.join();
        // if the display thread exits, ignore the computing threads
        // This may lead to crash when exiting.
        simplifyThread0.detach();
        simplifyThread1.detach();
        //simplifyThread0.join();
        //simplifyThread1.join();
    }
    else if (argc == 4) {

        ImprovedQuadricMethod improvedQuadric;

        auto vertexCallback = [&improvedQuadric](const REALTYPE x, const REALTYPE y, const REALTYPE z) { 
            improvedQuadric.insertVertex(x, y, z);
        };
        auto faceCallback = [&improvedQuadric](const int v0, const int v1, const int v2) { 
            improvedQuadric.insertFace(v0, v1, v2);
            ++totalNumFaces; 
        };

        // parse obj file
        ObjParser objParser(file, lowerBound, upperBound, vertexCallback, faceCallback);

        // calculate view point position
        objCenter = (lowerBound + upperBound) * 0.5;
        REALTYPE viewPointDistance = std::max(std::max(upperBound[0] - lowerBound[0], 
                                                       upperBound[1] - lowerBound[1]),
                                              upperBound[2] - lowerBound[2]) * 2.3;
        viewPoint = objCenter + Vector(viewPointDistance, 0, 0);

        // cannot insert any vertex or face after this
        improvedQuadric.switchOff();

        
        samplingThresholds.push_back(stold(argv[3]));
        
        std::chrono::high_resolution_clock clock;
        auto t1 = clock.now();
        // simplify the mesh
        improvedQuadric.simplify(samplingThresholds, results[1]);
        auto t2 = clock.now();
        cout << "Simplified in " << (t2 - t1).count() / 1e6 << " ms." << endl;

        // output results
        ofstream fout(argv[2]);
        objParser.output(results[1][0].first, fout);
        fout.close();

    }
    
}
