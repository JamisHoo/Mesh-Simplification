/******************************************************************************
 *  Copyright (c) 2014. All rights reserved.
 *
 *  Project: Mesh Simplification
 *  Filename: preview.h 
 *  Version: 1.0
 *  Author: Jinming Hu
 *  E-mail: hjm211324@gmail.com
 *  Date: Jun. 27, 2014
 *  Time: 15:34:27
 *  Description: functions to display results with OpenGL
 *****************************************************************************/
#ifndef PREVIEW_H
#define PREVIEW_H

#include <GL/freeglut.h>
#include "common.h"

namespace Preview {
    // numELE(vertexes) == numELE(normals) == 9 * numFacs
    // numELE: number of elements
    double* vertexes[2] = { nullptr, nullptr };
    double* normals[2] = { nullptr, nullptr };
    int numFaces[2] = { 0, 0 };
    int subWindow[2];
    double g_rotation = 0;

    MeshSimplification::Vector viewPoint, objCenter;
    MeshSimplification::Vector lowerBound, upperBound;
    std::vector< MeshSimplification::Result > results[2];
    std::vector< MeshSimplification::REALTYPE > samplingThresholds;

    int onDisplayNum[2] = { -1, -1 };

    // results.size == numSamples
    // results[0] is original mesh
    int numSamples = 100;
    int totalNumFaces = 0;

    // window size
    int winWidth = 1024;
    int winHeight = 640;


    // FIXME
    // States below aren't declared in freeglut_std.h version 2.0.
    // They work fine on my machine.
    // This may lead to some problems in the future,
    // but I don't have a better idea.
    const int GLUT_WHEEL_UP_ = 3;
    const int GLUT_WHEEL_DOWN_ = 4;

    void resizeCallBack(int, int) {
        // fix window size
        glutReshapeWindow(winWidth, winHeight);
    }

    void keyboardCallBack(unsigned char key, int, int) {
        if (key == 27) glutLeaveMainLoop();
    }

    void updateDisplayNum(const int i, bool up_down) { //up == 0 down == 1
        if (up_down == 0) // up
            if (onDisplayNum[i] > 0) --onDisplayNum[i];

        if (up_down == 1) // down
            if (onDisplayNum[i] + 1 < results[i].size()) ++onDisplayNum[i];


        // if typeid(MeshSimplification::REALTYPE) != typeid(double)
        /*
        delete[] vertexes;
        delete[] normals;
        vertexes = new double[results[onDisplayNum].first.size()];
        normals = new double[results[onDisplayNum].second.size()];

        std::copy(results[onDisplayNum].first.begin(), 
                  results[onDisplayNum].first.end(),
                  vertexes);
        std::copy(results[onDisplayNum].second.begin(),
                  results[onDisplayNum].second.end(),
                  normals);
        */
        // else
        if (onDisplayNum[i] == -1) 
            vertexes[i] = normals[i] = nullptr, numFaces[i] = 0;
        else {
            vertexes[i] = &results[i][onDisplayNum[i]].first[0];
            normals[i] = &results[i][onDisplayNum[i]].second[0];
            numFaces[i] = results[i][onDisplayNum[i]].first.size() / 9;
        }

    }

    void mouseCallBack(int button, int state, int, int) {
        if (state == GLUT_UP && (button == GLUT_WHEEL_UP_ || button == GLUT_WHEEL_DOWN_)) {
            if (button == GLUT_WHEEL_UP_) {
                updateDisplayNum(0, 0);
                updateDisplayNum(1, 0);
            }
            else if (button == GLUT_WHEEL_DOWN_) {
                updateDisplayNum(0, 1);
                updateDisplayNum(1, 1);
            }
            else assert(0);
        }
    }

    void display() {
        for (int i = 0; i < 2; ++i) {
            glutSetWindow(subWindow[i]);
            glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
            glLoadIdentity();
            gluLookAt(viewPoint[0], viewPoint[1], viewPoint[2],
                      objCenter[0], objCenter[1], objCenter[2],
                      0, 1, 0);
            glPushMatrix();
            
            glTranslatef(objCenter[0], objCenter[1], objCenter[2]);
            glRotatef(g_rotation, 0, 1, 0);
            glRotatef(90, 0, 1, 0);
            g_rotation += -0.75;
            glTranslatef(-objCenter[0], -objCenter[1], -objCenter[2]);
            
            if (onDisplayNum[i] == -1) 
                updateDisplayNum(i, 1);
            
            glEnableClientState(GL_VERTEX_ARRAY);
            glEnableClientState(GL_NORMAL_ARRAY);

            glVertexPointer(3, GL_DOUBLE, 0, vertexes[i]);
            glNormalPointer(GL_DOUBLE, 0, normals[i]);

            glDrawArrays(GL_TRIANGLES, 0, numFaces[i] * 3);

            glDisableClientState(GL_VERTEX_ARRAY);
            glDisableClientState(GL_NORMAL_ARRAY);
            
            glPopMatrix();

            double spacing = std::max(std::max(upperBound[0] - lowerBound[0],
                                               upperBound[1] - lowerBound[1]),
                                      upperBound[2] - lowerBound[2]);

            glRasterPos3f(upperBound[0], lowerBound[1] - spacing * 0.1, 0);
            std::string prompt = std::to_string(double(numFaces[i]) / totalNumFaces * 100);
            prompt = prompt.substr(0, 5);
            prompt += "%";
            for (auto c: prompt)
                glutBitmapCharacter(GLUT_BITMAP_HELVETICA_18, c);

            glRasterPos3f(upperBound[0], lowerBound[1] - spacing * 0.14, 0);
            prompt = std::to_string(numFaces[i]) + " faces";
            for (auto c: prompt)
                glutBitmapCharacter(GLUT_BITMAP_HELVETICA_18, c);

            glutSwapBuffers();
        }
    }


    void setWindowEnvironment(const int width, const int height) {
        int field_of_view_angle = 45;
        double z_near = 0.1;
        double z_far = 500.0;

        glutDisplayFunc(display);
        glutKeyboardFunc(keyboardCallBack);
        glutMouseFunc(mouseCallBack);
        glutReshapeFunc(resizeCallBack);
        glutIdleFunc(display);

        glMatrixMode(GL_PROJECTION);
        glViewport(0, 0, width / 2, height);
        GLfloat aspect = GLfloat(double(width / 2) / height);
        glMatrixMode(GL_PROJECTION);
        glLoadIdentity();
        gluPerspective(field_of_view_angle, aspect, z_near, z_far);
        glMatrixMode(GL_MODELVIEW);
        glShadeModel(GL_SMOOTH);
        glClearColor(0.379, 0.926, 0.242, 0.5);
        glClearDepth(1.0);
        glEnable(GL_DEPTH_TEST);
        glDepthFunc(GL_LEQUAL);
        glHint(GL_PERSPECTIVE_CORRECTION_HINT, GL_NICEST);
        GLfloat amb_light[] = { 0.01, 0.01, 0.01, 1.0 };
        GLfloat diffuse[] = { 0.3, 0.3, 0.99, 1.0 };
        GLfloat specular[] = {0.9, 0.9, 0.9, 1.0 };
        glLightModelfv(GL_LIGHT_MODEL_AMBIENT, amb_light);
        glLightfv(GL_LIGHT0, GL_DIFFUSE, diffuse);
        glLightfv(GL_LIGHT0, GL_SPECULAR, specular);
        glEnable(GL_LIGHT0);
        glEnable(GL_COLOR_MATERIAL);
        glShadeModel(GL_SMOOTH);
        glLightModeli(GL_LIGHT_MODEL_TWO_SIDE, GL_FALSE);
        glDepthFunc(GL_LEQUAL);
        glEnable(GL_DEPTH_TEST);
        glEnable(GL_LIGHTING);
        glEnable(GL_LIGHT0);
    }

    void initGL(int argc, char** argv) {
        // modify these variables to fit different objs
        int width = winWidth;
        int height = winHeight;
        std::string title = "Mesh Simplification";

        glutInit(&argc, argv);
        glutSetOption(GLUT_ACTION_ON_WINDOW_CLOSE, GLUT_ACTION_CONTINUE_EXECUTION);
        glutInitDisplayMode(GLUT_RGB | GLUT_DOUBLE | GLUT_DEPTH);
        glutInitWindowSize(width, height);

        int mainWindow = glutCreateWindow(title.c_str());

        //setWindowEnvironment(width, height, display);

        subWindow[0] = glutCreateSubWindow(mainWindow, 0, 0, winWidth / 2, winHeight);
        setWindowEnvironment(width, height);

        subWindow[1] = glutCreateSubWindow(mainWindow, winWidth / 2, 0, winWidth / 2, winHeight);
        setWindowEnvironment(width, height);

        glutMainLoop();
    }
}

#endif /* PREVIEW_H */
