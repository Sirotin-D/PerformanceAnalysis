//
// Created by Дмитрий on 23.12.2022.
//
#pragma once
#include <cmath>
#include <chrono>
#include <iostream>
#include <vector>
#include "StartConditions.h"

void slowSolution(int n, int m, double Nmax, double eps) {

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

    std::vector<std::vector<double>> V;
    V.assign((m + 1), std::vector<double>(n + 1, 0));

    for (int i = 0; i < m + 1; i++) {
        V[i][0] = exactSolution(a, y[i]);
        V[i][n] = exactSolution(b, y[i]);
    }

    for (int i = 1; i < n; i++) {
        V[0][i] = exactSolution(x[i], c);
        V[m][i] = exactSolution(x[i], d);
    }

    std::vector<std::vector<double>> A;
    A.assign((n - 1) * (m - 1), std::vector<double>((n - 1) * (m - 1), 0));

    std::vector<double> fLaplas((m - 1) * (n - 1), 0);

    for (int i = 1, row = 0; i < m; i++) {
        for (int j = 1; j < n; j++, row++) {

            bool isMatrixFilled = true;

            if (j == 1) {
                A[row][(i - 1) * (n - 1) + j - 1] = centerParam;
                fLaplas[row] = -laplasian(x[j], y[i]);
                fLaplas[row] -= paramX * V[i][j - 1];
                A[row][(i - 1) * (n - 1) + j + 1 - 1] = paramX;
                if (i == 1) {
                    A[row][(i + 1 - 1) * (n - 1) + j - 1] = paramY;
                    fLaplas[row] -= paramY * V[i - 1][j];
                }
                else {
                    if (i == m - 1) {
                        A[row][(i - 1 - 1) * (n - 1) + j - 1] = paramY;
                        fLaplas[row] -= paramY * V[i + 1][j];
                    }
                    else {
                        A[row][(i + 1 - 1) * (n - 1) + j - 1] = paramY;
                        A[row][(i - 1 - 1) * (n - 1) + j - 1] = paramY;
                    }
                }
                isMatrixFilled = false;
            }
            if (j == (n - 1)) {
                A[row][(i - 1) * (n - 1) + j - 1] = centerParam;
                fLaplas[row] = -laplasian(x[j], y[i]);
                fLaplas[row] -= paramX * V[i][j + 1];
                A[row][(i - 1) * (n - 1) + j - 1 - 1] = paramX;
                if (i == 1) {
                    A[row][(i + 1 - 1) * (n - 1) + j - 1] = paramY;
                    fLaplas[row] -= paramY * V[i - 1][j];
                }
                else {
                    if (i == m - 1) {
                        A[row][(i - 1 - 1) * (n - 1) + j - 1] = paramY;
                        fLaplas[row] -= paramY * V[i + 1][j];
                    }
                    else {
                        A[row][(i + 1 - 1) * (n - 1) + j - 1] = paramY;
                        A[row][(i - 1 - 1) * (n - 1) + j - 1] = paramY;
                    }
                }
                isMatrixFilled = false;
            }
            if ((i == 1) && ((j > 1) && (j < n - 1))) {
                A[row][(i - 1) * (n - 1) + j - 1] = centerParam;
                fLaplas[row] = -laplasian(x[j], y[i]);
                fLaplas[row] -= paramY * V[i - 1][j];
                A[row][(i + 1 - 1) * (n - 1) + j - 1] = paramY;
                A[row][(i - 1) * (n - 1) + j - 1 - 1] = paramX;
                A[row][(i - 1) * (n - 1) + j + 1 - 1] = paramX;
                isMatrixFilled = false;
            }
            if ((i == (m - 1)) && ((j > 1) && (j < n - 1))) {
                A[row][(i - 1) * (n - 1) + j - 1] = centerParam;
                fLaplas[row] = -laplasian(x[j], y[i]);
                fLaplas[row] -= paramY * V[i + 1][j];
                A[row][(i - 1 - 1) * (n - 1) + j - 1] = paramY;
                A[row][(i - 1) * (n - 1) + j - 1 - 1] = paramX;
                A[row][(i - 1) * (n - 1) + j + 1 - 1] = paramX;
                isMatrixFilled = false;

            }
            if (isMatrixFilled) {
                A[row][(i - 1) * (n - 1) + j - 1] = centerParam;
                A[row][(i + 1 - 1) * (n - 1) + j - 1] = paramY;
                A[row][(i - 1 - 1) * (n - 1) + j - 1] = paramY;
                A[row][(i - 1) * (n - 1) + j + 1 - 1] = paramX;
                A[row][(i - 1) * (n - 1) + j - 1 - 1] = paramX;
                fLaplas[row] = -laplasian(x[j], y[i]);
            }
        }
    }

    for (int i = 0; i < A.size(); i++) {
        for (int j = 0; j < A[0].size(); j++) {
            A[i][j] = -A[i][j];
        }
    }

    for (int i = 0; i < fLaplas.size(); i++) {
        fLaplas[i] = -fLaplas[i];
    }

    int numStep = 0;
    double epsMax, currentEpr;
    std::vector<double> V_internal((n - 1) * (m - 1), 0);
    double vOld, vNew;
    auto start = std::chrono::high_resolution_clock::now();
    while (true) {
        epsMax = 0.0;
        for (int i = 0; i < V_internal.size(); i++) {
            vOld = V_internal[i];
            vNew = fLaplas[i];
            for (int j = 0; j < V_internal.size(); j++)
                if (j != i) {
                    vNew -= A[i][j] * V_internal[j];
                }
            vNew /= A[i][i];
            currentEpr = fabs(vOld - vNew);
            if (currentEpr > epsMax) {
                epsMax = currentEpr;
            }
            V_internal[i] = vNew;
        }
        numStep++;
        if ((epsMax < eps) || (numStep >= Nmax)) {
            break;
        }
    }
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration <double> duration = end - start;
    double runtime = duration.count();

    for (int i = 1, s = 0; i < m; i++) {
        for (int j = 1; j < n; j++, s++) {
            V[i][j] = V_internal[s];
        }
    }

    std::vector<double> discr(V_internal.size(), 0);

    double discrMax = 0.0;
    for (int i = 0; i < V_internal.size(); i++) {
        for (int j = 0; j < V_internal.size(); j++) {
            discr[i] += A[i][j] * V_internal[j];
        }
        discr[i] -= fLaplas[i];
        if (fabs(discr[i]) > discrMax) {
            discrMax = fabs(discr[i]);
        }
    }

    double schematicErrorMax = 0.0;
    for (int i = 0; i < V.size(); i++) {
        for (int j = 0; j < V[0].size(); j++) {
            double tmp = fabs(V[i][j] - exactSolution(x[j], y[i]));
            if (tmp > schematicErrorMax) {
                schematicErrorMax = tmp;
            }
        }
    }

    std::cout << "End calculation" << std::endl;
    std::cout << "Total iterations: " << numStep << std::endl;
    std::cout << "Achievement accuracy: " << epsMax << std::endl;
    std::cout << "Schematic Error: " << schematicErrorMax << std::endl;
    std::cout << "Discrepancy norma: " << discrMax << std::endl;
    std::cout << "Total time: " << runtime << std::endl;
}