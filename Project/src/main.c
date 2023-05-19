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
    femGeo* theGeometry = geoGetGeometry();   
    geoMeshRead("../data/mesh.txt");
    femProblem* theProblem = femElasticityRead(theGeometry,"../data/problem.txt");
    femElasticityPrint(theProblem);
    double *theSoluce = femElasticitySolve(theProblem); 
    femNodes *theNodes = theGeometry->theNodes;
    femFieldWrite(theNodes->nNodes,2,&theSoluce[0],"../data/U.txt");
    femFieldWrite(theNodes->nNodes,2,&theSoluce[1],"../data/V.txt");
    femElasticityFree(theProblem); 
    geoFree();



    double Lx = 1.0;
    double Ly = 1.0;
    double Bx = 0.5;
    double By = 0.0;    
      
    geoInitialize();
    femGeo* theGeometry = geoGetGeometry();
    
    theGeometry->LxPlate     =  Lx;
    theGeometry->LyPlate     =  Ly; 
    theGeometry->LxBar       =  Bx;
    theGeometry->LyBar       =  By;        
    theGeometry->h           =  Lx * 0.25;    
    theGeometry->elementType = FEM_TRIANGLE;  //FEM_TRIANGLE  FEM_QUAD
  
//    geoMeshGenerateGeo();
//    geoMeshGenerate();
    geoMeshGenerateGeoProjet();
//    geoMeshGenerateProjecteur();

//  geoMeshGenerateBar();   
    

    geoMeshImport();
//    geoSetDomainName(0,"gauche");  //DIRICHLET_X
//    geoSetDomainName(7,"Bottom");    //DIRICHLET_Y


   //pour geoMeshGenerateGeoProjet();
    geoSetDomainName(0,"petitGauche");
    geoSetDomainName(1,"platGauche");
    geoSetDomainName(2,"gauche");  //DIRICHLET_X    /!\ -1
    geoSetDomainName(3,"dessus");    //DIRICHLET_Y
    geoSetDomainName(4,"droite");    //DIRICHLET_X
    geoSetDomainName(5,"platDroite");  //DIRICHLET_X    /!\ -1
    geoSetDomainName(6,"petitDroite");    //DIRICHLET_Y
    geoSetDomainName(7,"dessous");    //DIRICHLET_X



/*
//juste la bar
    geoSetDomainName(0,"gauche");  //DIRICHLET_X    /!\ -1
    geoSetDomainName(1,"dessus");    //DIRICHLET_Y
    geoSetDomainName(2,"droite");    //DIRICHLET_X    
  
  */  
        
//
//  -2- Creation probleme 
//
    
    double E   = 211.e10;
    double nu  = 0.3;
    double rho = 7.85e3; 
    double g   = 9.81;
    femProblem* theProblem = femElasticityCreate(theGeometry,E,nu,rho,g,PLANAR_STRESS);  // PLANAR_STRAIN PLANAR_STRESS AXISYM
    //femElasticityAddBoundaryCondition(theProblem,"petitGauche",DIRICHLET_X,0.0);              //DIRICHLET_X,DIRICHLET_Y,NEUMANN_X,NEUMANN_Y,TANGENTIEL,NORMAL
    //femElasticityAddBoundaryCondition(theProblem,"platGauche", DIRICHLET_Y ,0.0);
    femElasticityAddBoundaryCondition(theProblem,"gauche",DIRICHLET_X,0.0);
    femElasticityAddBoundaryCondition(theProblem,"dessus", DIRICHLET_Y ,0.0);
    femElasticityAddBoundaryCondition(theProblem,"droite",DIRICHLET_X,0.0);              //DIRICHLET_X,DIRICHLET_Y,NEUMANN_X,NEUMANN_Y,TANGENTIEL,NORMAL
    //femElasticityAddBoundaryCondition(theProblem,"platDroite", DIRICHLET_Y ,0.0);
    //femElasticityAddBoundaryCondition(theProblem,"petitDroite",NEUMANN_Y,500000.0);
    //femElasticityAddBoundaryCondition(theProblem,"dessous", NEUMANN_X , 900.0);
    femElasticityPrint(theProblem);
    double *theSoluce = femElasticitySolve(theProblem); 







    return 0;  
}

 
