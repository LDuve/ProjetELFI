/*
 *  main.c
 *  Projet 2022-2023
 *  Elasticite lineaire plane
 *
 *  Code de calcul
 *
 *  Copyright (C) 2023 UCL-IMMC : Vincent Legat
 *  All rights reserved.
 *
 */
 
#include "fem.h"

int main(void)
{  

    double Lx = 1.0;
    double Ly = 1.0;
    double Bx = 0.5;
    double By = 0.0;    

    femGeo* theGeometry = geoGetGeometry();   

    theGeometry->LxPlate     =  Lx;
    theGeometry->LyPlate     =  Ly;
    theGeometry->LxBar       =  Bx;
    theGeometry->LyBar       =  By;
    theGeometry->h           =  Lx * 0.25;
    theGeometry->elementType = FEM_TRIANGLE;  //FEM_TRIANGLE  FEM_QUAD

    geoMeshRead("../data/mesh.txt");    
    
    geoSetDomainName(1,"Sym");
    geoSetDomainName(2,"Bot");

    femProblem* theProblem = femElasticityRead(theGeometry,"../data/problem.txt");    
    
    femElasticityAddBoundaryCondition(theProblem,"Sym",DIRICHLET_Y,0.0);
    femElasticityAddBoundaryCondition(theProblem,"Bot",DIRICHLET_X,0.0);

    //Numéro de solveur 1:Full, 2:Bande, 3:Frontal
    int solvertype =2;
    //Numéro de problème 5:normal, 6:axysimétrique
    int symétrie = 5;
    
    double* theSoluce;
    femElasticityPrint(theProblem);
    if(symétrie==5){
        theSoluce = femElasticitySolve(theProblem, solvertype);
    }
    else{theSoluce = femElasticitySolveSym(theProblem, solvertype);}

    femNodes *theNodes = theGeometry->theNodes;


    femFieldWrite(theNodes->nNodes,2,&theSoluce[0],"../data/U.txt");
    femFieldWrite(theNodes->nNodes,2,&theSoluce[1],"../data/V.txt");
    femElasticityFree(theProblem); 
    geoFree();

    return 0;  
}
