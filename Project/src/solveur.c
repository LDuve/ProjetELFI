/*
 *  basé sur fem.c
 *  Library for LEPL1110 : Finite Elements for dummies
 *
 *  Copyright (C) 2021 UCL-IMMC : Vincent Legat
 *  All rights reserved.
 *
 */

#include "fem.h"


double* femFullSystemEliminateBande(femFullSystem *mySystem)
{
    double  **A, *B, factor;
    int     i, j, k, jend, size, band;
    A    = mySystem->A;
    B    = mySystem->B;
    size = mySystem->size;
    band = calculateBandwidth(A, size);

    /* Incomplete Cholesky factorization */ 

    for (k=0; k < size; k++) {
        if ( fabs(A[k][k]) <= 1e-4 ) {
            Error("Cannot eleminate with such a pivot"); }
        jend = fmin(k + band,size);
        for (i = k+1 ; i <  jend; i++) {
            factor = A[k][i] / A[k][k];
            for (j = i ; j < jend; j++) 
                A[i][j] = A[i][j] - A[k][j] * factor;
            B[i] = B[i] - B[k] * factor; }}
        
    /* Back-substitution */

    for (i = (size-1); i >= 0 ; i--) {
        factor = 0;
        jend = fmin(i + band,size);
        for (j = i+1 ; j < jend; j++)
            factor += A[i][j] * B[j];
        B[i] = ( B[i] - factor)/A[i][i]; }

    return(mySystem->B);
}

double* femFullSystemEliminate(femFullSystem *mySystem)
{
    double  **A, *B, factor;
    int     i, j, k, size;
    
    A    = mySystem->A;
    B    = mySystem->B;
    size = mySystem->size;
    
    /* Gauss elimination */
    
    for (k=0; k < size; k++) {
        if ( fabs(A[k][k]) <= 1e-16 ) {
            printf("Pivot index %d  ",k);
            printf("Pivot value %e  ",A[k][k]);
            Error("Cannot eliminate with such a pivot"); }
        for (i = k+1 ; i <  size; i++) {
            factor = A[i][k] / A[k][k];
            for (j = k+1 ; j < size; j++) 
                A[i][j] = A[i][j] - A[k][j] * factor;
            B[i] = B[i] - B[k] * factor; }}
    
    /* Back-substitution */
    
    for (i = size-1; i >= 0 ; i--) {
        factor = 0;
        for (j = i+1 ; j < size; j++)
            factor += A[i][j] * B[j];
        B[i] = ( B[i] - factor)/A[i][i]; }
    
    return(mySystem->B);    
}

double* femFullSystemEliminateFrontal(femFullSystem *mySystem)
{
    double  **A, *B, factor;
    int     i, j, k, size;
    
    A    = mySystem->A;
    B    = mySystem->B;
    size = mySystem->size;

    /* Frontal elimination */

    for (k = 0; k < size; k++) {
        if (fabs(A[k][k]) <= 1e-16) {
            printf("Pivot index %d  ", k);
            printf("Pivot value %e  ", A[k][k]);
            Error("Cannot eliminate with such a pivot");
        }

        for (i = k + 1; i < size; i++) {
            if (A[i][k] != 0) {
                factor = A[i][k] / A[k][k];

                for (j = k + 1; j < size; j++)
                    A[i][j] -= factor * A[k][j];

                B[i] -= factor * B[k];
            }
        }
    }

    /* Back-substitution */

    for (i = size - 1; i >= 0; i--) {
        factor = 0;
        for (j = i + 1; j < size; j++)
            factor += A[i][j] * B[j];
        B[i] = (B[i] - factor) / A[i][i];
    }

    return mySystem->B;
}

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

double *femElasticitySolve(femProblem *theProblem, int solver)
{

    femFullSystem  *theSystem = theProblem->system;
    femIntegration *theRule = theProblem->rule;
    femDiscrete    *theSpace = theProblem->space;
    femGeo         *theGeometry = theProblem->geometry;
    femNodes       *theNodes = theGeometry->theNodes;
    femMesh        *theMesh = theGeometry->theElements;
    
    
    double x[4],y[4],phi[4],dphidxsi[4],dphideta[4],dphidx[4],dphidy[4];
    int iElem,iInteg,iEdge,i,j,d,map[4],mapX[4],mapY[4];
    
    int nLocal = theMesh->nLocalNode;

    double a   = theProblem->A;
    double b   = theProblem->B;
    double c   = theProblem->C;      
    double rho = theProblem->rho;
    double g   = theProblem->g;
    double **A = theSystem->A;
    double *B  = theSystem->B;
    
    
    for (iElem = 0; iElem < theMesh->nElem; iElem++) {
        for (j=0; j < nLocal; j++) {
            map[j]  = theMesh->elem[iElem*nLocal+j];
            mapX[j] = 2*map[j];
            mapY[j] = 2*map[j] + 1;
            x[j]    = theNodes->X[map[j]];
            y[j]    = theNodes->Y[map[j]];} 
        
        for (iInteg=0; iInteg < theRule->n; iInteg++) {    
            double xsi    = theRule->xsi[iInteg];
            double eta    = theRule->eta[iInteg];
            double weight = theRule->weight[iInteg];  
            femDiscretePhi2(theSpace,xsi,eta,phi);
            femDiscreteDphi2(theSpace,xsi,eta,dphidxsi,dphideta);
            
            double dxdxsi = 0.0;
            double dxdeta = 0.0;
            double dydxsi = 0.0; 
            double dydeta = 0.0;
            for (i = 0; i < theSpace->n; i++) {  
                dxdxsi += x[i]*dphidxsi[i];       
                dxdeta += x[i]*dphideta[i];   
                dydxsi += y[i]*dphidxsi[i];   
                dydeta += y[i]*dphideta[i]; }
            double jac = fabs(dxdxsi * dydeta - dxdeta * dydxsi);
            
            for (i = 0; i < theSpace->n; i++) {    
                dphidx[i] = (dphidxsi[i] * dydeta - dphideta[i] * dydxsi) / jac;       
                dphidy[i] = (dphideta[i] * dxdxsi - dphidxsi[i] * dxdeta) / jac; }            
            for (i = 0; i < theSpace->n; i++) { 
                for(j = 0; j < theSpace->n; j++) {
                    A[mapX[i]][mapX[j]] += (dphidx[i] * a * dphidx[j] + dphidy[i] * c * dphidy[j]) * jac * weight;                                                                                            
                    A[mapX[i]][mapY[j]] += (dphidx[i] * b * dphidy[j] + dphidy[i] * c * dphidx[j]) * jac * weight;                                                                                           
                    A[mapY[i]][mapX[j]] += (dphidy[i] * b * dphidx[j] + dphidx[i] * c * dphidy[j]) * jac * weight;                                                                                            
                    A[mapY[i]][mapY[j]] += (dphidy[i] * a * dphidy[j] + dphidx[i] * c * dphidx[j]) * jac * weight; }}
             for (i = 0; i < theSpace->n; i++) {
                B[mapY[i]] -= phi[i] * g * rho * jac * weight; }}} 
    

    int *theConstrainedNodes = theProblem->constrainedNodes;

    for (int i=0; i < theSystem->size; i++) {
        
        if (theConstrainedNodes[i] != -1  && (theProblem->conditions[theConstrainedNodes[i]]->type == DIRICHLET_X ||  theProblem->conditions[theConstrainedNodes[i]]->type == DIRICHLET_Y)) { 
            double value = theProblem->conditions[theConstrainedNodes[i]]->value;
            femFullSystemConstrain(theSystem,i,value); }
        
        
        if (theConstrainedNodes[i] != -1 && (theProblem->conditions[theConstrainedNodes[i]]->type == NEUMANN_X ||  theProblem->conditions[theConstrainedNodes[i]]->type == NEUMANN_Y)) { 
            double GoodValue = theProblem->conditions[theConstrainedNodes[i]]->domain->elemUnique->value ;
            femFullSystemConstrainNeumann(theSystem,i,GoodValue);
            
            }
        }
    
    if(solver==1){
        return femFullSystemEliminate(theSystem);
    }
    if(solver==2){
        return femFullSystemEliminateBande(theSystem);
    }
    if(solver ==3){
        return femFullSystemEliminateFrontal(theSystem);
    }
}

double *femElasticitySolveSym(femProblem *theProblem, int solver)
{

    femFullSystem  *theSystem = theProblem->system;
    femIntegration *theRule = theProblem->rule;
    femDiscrete    *theSpace = theProblem->space;
    femGeo         *theGeometry = theProblem->geometry;
    femNodes       *theNodes = theGeometry->theNodes;
    femMesh        *theMesh = theGeometry->theElements;
    
    
    double x[4],y[4],phi[4],dphidxsi[4],dphideta[4],dphidx[4],dphidy[4];
    int iElem,iInteg,iEdge,i,j,d,map[4],mapX[4],mapY[4];
    
    int nLocal = theMesh->nLocalNode;
    double r = theGeometry->rayon;
    double a   = theProblem->A;
    double b   = theProblem->B;
    double c   = theProblem->C;      
    double rho = theProblem->rho;
    double g   = theProblem->g;
    double **A = theSystem->A;
    double *B  = theSystem->B;
    
    
    for (iElem = 0; iElem < theMesh->nElem; iElem++) {
        for (j=0; j < nLocal; j++) {
            map[j]  = theMesh->elem[iElem*nLocal+j];
            mapX[j] = 2*map[j];
            mapY[j] = 2*map[j] + 1;
            x[j]    = theNodes->X[map[j]];
            y[j]    = theNodes->Y[map[j]];} 
        
        for (iInteg=0; iInteg < theRule->n; iInteg++) {    
            double xsi    = theRule->xsi[iInteg];
            double eta    = theRule->eta[iInteg];
            double weight = theRule->weight[iInteg];  
            femDiscretePhi2(theSpace,xsi,eta,phi);
            femDiscreteDphi2(theSpace,xsi,eta,dphidxsi,dphideta);
            
            double dxdxsi = 0.0;
            double dxdeta = 0.0;
            double dydxsi = 0.0; 
            double dydeta = 0.0;

            for (i = 0; i < theSpace->n; i++) {  
                dxdxsi += x[i]*dphidxsi[i];       
                dxdeta += x[i]*dphideta[i];   
                dydxsi += y[i]*dphidxsi[i];   
                dydeta += y[i]*dphideta[i]; }
            double jac = fabs(dxdxsi * dydeta - dxdeta * dydxsi);
            
            for (i = 0; i < theSpace->n; i++) {    
                dphidx[i] = (dphidxsi[i] * dydeta - dphideta[i] * dydxsi) / jac;       
                dphidy[i] = (dphideta[i] * dxdxsi - dphidxsi[i] * dxdeta) / jac; }            
            for (i = 0; i < theSpace->n; i++) { 
                for(j = 0; j < theSpace->n; j++) {
                    A[mapX[i]][mapX[j]] += (dphidx[i] * a * dphidx[j]*r + dphidy[i] * c * dphidy[j]*r+ b  *dphidx[i] * phi[j] + phi[i] *(b*dphidx[j]+a*phi[j]/r)) * jac * weight*3.14;                                                                                            
                    A[mapX[i]][mapY[j]] += (dphidx[i] * b * dphidy[j]*r + dphidy[i] * c * dphidx[j]*r+b*phi[i]*dphidy[j]) * jac * weight*3.14;                                                                                           
                    A[mapY[i]][mapX[j]] += (dphidy[i] * b * dphidx[j] *r+ dphidx[i] * c * dphidy[j]*r + b* dphidy[i]*phi[j]) * jac * weight*3.14;                                                                                            
                    A[mapY[i]][mapY[j]] += (dphidy[i] * a * dphidy[j] *r+ dphidx[i] * c * dphidx[j]*r) * jac * weight*3.14; }}
             for (i = 0; i < theSpace->n; i++) {
                B[mapY[i]] -= phi[i] * g * rho * jac * weight; }}} 

    int *theConstrainedNodes = theProblem->constrainedNodes;
    for (int i=0; i < theSystem->size; i++) {
        
        if (theConstrainedNodes[i] != -1  && (theProblem->conditions[theConstrainedNodes[i]]->type == DIRICHLET_X ||  theProblem->conditions[theConstrainedNodes[i]]->type == DIRICHLET_Y)) { //correspondant au nœud actuel que nous examinons. si different de -1 alors il y a une contraite 
            printf("RENTRE  Dirichlet!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n");
            double value = theProblem->conditions[theConstrainedNodes[i]]->value;
            femFullSystemConstrain(theSystem,i,value); }
        
        
        if (theConstrainedNodes[i] != -1 && (theProblem->conditions[theConstrainedNodes[i]]->type == NEUMANN_X ||  theProblem->conditions[theConstrainedNodes[i]]->type == NEUMANN_Y)) { 
            printf("RENTRE Neumannnn!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n");
            double GoodValue = theProblem->conditions[theConstrainedNodes[i]]->domain->elemUnique->value ;
            femFullSystemConstrainNeumann(theSystem,i,GoodValue);            
            }
        }
                 
    if(solver==1){
        return femFullSystemEliminate(theSystem);
    }
    if(solver==2){
        return femFullSystemEliminateBande(theSystem);
    }
    if(solver ==3){
        return femFullSystemEliminateFrontal(theSystem);
    }
}
