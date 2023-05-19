/*
 *  basé sur fem.c
 *  Library for LEPL1110 : Finite Elements for dummies
 *
 *  Copyright (C) 2021 UCL-IMMC : Vincent Legat
 *  All rights reserved.
 *
 */

#include "fem.h"

int calculateBandwidth(double **A, int size) {
    int bandwidth = 0;    
    for (int i = 0; i < size; i++) {
        for (int j = 0; j < size; j++) {
            if (A[i][j] != 0) {
                int distance = abs(j - i);
                if (distance > bandwidth) {
                    bandwidth = distance;
                }
            }
        }
    }
    
    return bandwidth;
}