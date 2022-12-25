//
// Created by Дмитрий on 23.12.2022.
//
#pragma once
#include <cmath>

static double laplasian(double _x, double _y) {
    return 4 * (1 - pow(_x, 2) - pow(_y, 2)) * exp(1 - pow(_x, 2) - pow(_y, 2));
}

static double exactSolution(double _x, double _y) {
    return exp(1 - pow(_x, 2) - pow(_y, 2));
}
