//
// Created by Дмитрий on 23.12.2022.
//
#pragma once
#include <cmath>
#include <vector>
#include <chrono>
#include <iostream>
#include "StartConditions.h"

void optimizedSolution(int n, int m, double Nmax, double eps) {

    double a = -1, b = 1;
    double c = -1, d = 1;
    double h = (b - a) / n;
    double k = (d - c) / m;
    double paramX = 1 / pow(h, 2);
    double paramY = 1 / pow(k, 2);
    double centerParam = -2 * (paramX + paramY);
    std::vector<double> x(n + 1);
    std::vector<double> y(m + 1);
    for (int i = 0; i < x.size(); i++) {
        x[i] = a + i * h;
    }

    for (int i = 0; i < y.size(); i++) {
        y[i] = c + i * k;
    }

    std::vector <std::vector <double>> V;
    V.assign((m+1), std::vector <double>(n + 1, 0));
    for (int i = 0; i < m + 1; i++) {
        V[i][0] = exactSolution(a, y[i]);
        V[i][n] = exactSolution(b, y[i]);
    }
    for (int i = 1; i < n; i++) {
        V[0][i] = exactSolution(x[i], c);
        V[m][i] = exactSolution(x[i], d);
    }

    std::vector<std::vector<double>> fLaplas;
    fLaplas.assign((m + 1), std::vector <double>(n + 1, 0));

    for (int i = 1; i < m; i++) {
        for (int j = 1; j < n; j++) {
            fLaplas[i][j] = laplasian(x[j], y[i]);
        }
    }
    int numStep = 0;
    double epsMax, currentEps;
    double vOld, vNew;
    auto start = std::chrono::high_resolution_clock::now();
    while (true) {
        epsMax = 0.0;
        for (int i = 1; i < m; i++) {
            for (int j = 1; j < n; j++) {
                vOld = V[i][j];
                vNew = centerParam * vOld + paramX * (V[i][j - 1] + V[i][j + 1]) + paramY * (V[i - 1][j] + V[i + 1][j]);
                V[i][j] = vOld - (vNew + fLaplas[i][j]) / centerParam;
                currentEps = fabs(V[i][j] - vOld);
                if (currentEps > epsMax) {
                    epsMax = currentEps;
                }
            }
        }
        numStep++;
        if ((epsMax < eps) || (numStep >= Nmax)) {
            break;
        }
    }
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration <double> duration = end - start;
    double runtime = duration.count();

    std::vector <std::vector <double>> discr;
    discr.assign((m + 1), std::vector <double>(n + 1, 0));

    double discrMax = 0.0, temp;
    for (int i = 1; i < m; i++) {
        for (int j = 1; j < n; j++) {
            discr[i][j] = centerParam * V[i][j] + paramX * (V[i][j - 1] + V[i][j + 1]) + paramY * (V[i - 1][j] + V[i + 1][j]) + fLaplas[i][j];
            temp = fabs(discr[i][j]);
            if (temp > discrMax) {
                discrMax = temp;
            }
        }
    }

    double schematicErrorMax = 0.0;
    for (int i = 0; i < V.size(); i++) {
        for (int j = 0; j < V[0].size(); j++) {
            temp = fabs(V[i][j] - exactSolution(x[j], y[i]));
            if (temp > schematicErrorMax) {
                schematicErrorMax = temp;
            }
        }
    }
    std::cout << "End slow algorithm calculation" << std::endl;
    std::cout << "Total iterations: " << numStep << std::endl;
    std::cout << "Achievement accuracy: " << epsMax << std::endl;
    std::cout << "Schematic Error: " << schematicErrorMax << std::endl;
    std::cout << "Discrepancy norma: " << discrMax << std::endl;
    std::cout << "Total time: " << runtime << std::endl;
}
