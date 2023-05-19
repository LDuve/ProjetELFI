/*
 *  basé sur fem.c
 *  Library for LEPL1110 : Finite Elements for dummies
 *
 *  Copyright (C) 2021 UCL-IMMC : Vincent Legat
 *  All rights reserved.
 *
 */
#include "fem.h"

void  femFullSystemConstrain(femFullSystem *mySystem,  int myNode, double myValue) 
{
    double  **A, *B;
    int     i, size;
    
    A    = mySystem->A;
    B    = mySystem->B;
    size = mySystem->size;
    
    for (i=0; i < size; i++) {
        B[i] -= myValue * A[i][myNode];
        A[i][myNode] = 0; }
    
    for (i=0; i < size; i++) 
        A[myNode][i] = 0; 
    
    A[myNode][myNode] = 1;
    B[myNode] = myValue;
}

void  femFullSystemConstrainNeumann(femFullSystem *mySystem, int myNode, double myValue) 
{
    double  *B;        

    B    = mySystem->B;
    
    B[myNode] += myValue;
}



void femElasticityAddBoundaryCondition(femProblem *theProblem, char *nameDomain, femBoundaryType type, double value)
{
    int iDomain = geoGetDomain(nameDomain);
    if (iDomain == -1)  Error("Undefined domain :-(");

    femBoundaryCondition* theBoundary = malloc(sizeof(femBoundaryCondition));
    theBoundary->domain = theProblem->geometry->theDomains[iDomain];
    theBoundary->value = value;
    theBoundary->type = type;  
    theProblem->nBoundaryConditions++;
    int size = theProblem->nBoundaryConditions; 

    
    if (theProblem->conditions == NULL)
        theProblem->conditions = malloc(size*sizeof(femBoundaryCondition*));
    else 
        theProblem->conditions = realloc(theProblem->conditions, size*sizeof(femBoundaryCondition*));
    theProblem->conditions[size-1] = theBoundary;
    
    int *elem = theBoundary->domain->elem;
    int nElem = theBoundary->domain->nElem; 

    theBoundary->domain->elemUnique = malloc(nElem * sizeof(femElem*));
    for (int i = 0; i < nElem; i++) {
        theBoundary->domain->elemUnique[i].elemIndex = elem[i];
    }

    if(type == DIRICHLET_X || type == DIRICHLET_Y ){
        int shift;
        if (type == DIRICHLET_X)  shift = 0;   
        if (type == DIRICHLET_Y)  shift = 1;    
        
        if (type == DIRICHLET_X || type == DIRICHLET_Y){
            for (int e=0; e<nElem; e++) {
                for (int i=0; i<2; i++) {  
                    int node = theBoundary->domain->mesh->elem[2*elem[e]+i];  
                    theProblem->constrainedNodes[2*node+shift] = size-1;  
                    }
                }    
            }
        }
    
    if(type == NEUMANN_X || type == NEUMANN_Y ){
        double val = theBoundary->value ;
        int shiftN;
        if (type == NEUMANN_X)  shiftN = 0;   
        if (type == NEUMANN_Y)  shiftN = 1;

        for (int e=0; e<nElem; e++) {
            for (int i=0; i<2; i++) {                    
                    int node = theBoundary->domain->mesh->elem[2*elem[e]+i]; 
                    theProblem->constrainedNodes[2*node+shiftN] = size-1;    
                }
            
            // Calcul de la différence de distance pour chaque élément
            double x1 = theProblem->geometry->theNodes->X[theBoundary->domain->mesh->elem[2 * elem[e]]];
            double y1 = theProblem->geometry->theNodes->Y[theBoundary->domain->mesh->elem[2 * elem[e]]];
            double x2 = theProblem->geometry->theNodes->X[theBoundary->domain->mesh->elem[2 * elem[e] + 1]];
            double y2 = theProblem->geometry->theNodes->Y[theBoundary->domain->mesh->elem[2 * elem[e] + 1]];
            double distance = sqrt(pow(x2 - x1, 2) + pow(y2 - y1, 2));
            
            theBoundary->domain->elemUnique[e].value = val * distance * 0.5;
            }
        }      
}
