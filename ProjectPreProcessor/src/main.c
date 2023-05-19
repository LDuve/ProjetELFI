/*
 *  main.c
 *  Projet 2022-2023
 *  Elasticite lineaire plane
 *
 *  Preprocesseur
 *
 *  Copyright (C) 2023 UCL-IMMC : Vincent Legat
 *  All rights reserved.
 *
 */
 
#include "glfem.h"

int main(void)
{  

//
//  -1- Construction de la geometrie 
//

    double Lx = 1.0;
    double Ly = 1.0;   
    double LxBx = 0.5;
    double LyBy = 2.0;  
    geoInitialize();
    femGeo* theGeometry = geoGetGeometry();
    
    theGeometry->LxPlate     =  Lx;
    theGeometry->LyPlate     =  Ly;  
    theGeometry->LxBar       =  LxBx;
    theGeometry->LyBar       =  LyBy;   
    theGeometry->h           =  Lx * 0.05;    
    theGeometry->elementType = FEM_QUAD;
  
//    geoMeshGenerate();      // Utilisation de OpenCascade
    
//  geoMeshGenerateGeo();   // Utilisation de outils de GMSH  
                            // Attention : les entit�s sont diff�rentes !
                            // On a aussi invers� la g�omtrie pour rire !
                            
    geoMeshGenerateGeoFile("../data/mesh.geo");   // Lecture fichier geo
  
    geoMeshImport();
 
   //pour geoMeshGenerateGeoProjet
    geoSetDomainName(0,"petitGauche");
    geoSetDomainName(1,"platGauche");
    geoSetDomainName(2,"gauche"); 
    geoSetDomainName(3,"dessus");    
    geoSetDomainName(4,"droite");   
    geoSetDomainName(5,"platDroite"); 
    geoSetDomainName(6,"petitDroite");    
    geoSetDomainName(7,"dessous");   

/*
    //pour geoMeshGenerateBar
    geoSetDomainName(0,"gauche");
    geoSetDomainName(1,"dessus");
    geoSetDomainName(2,"droite"); 
    geoSetDomainName(3,"dessus");    
*/
/*
   //pour geoMeshGenerateProjecteur
    geoSetDomainName(0,"petitGauche");
    geoSetDomainName(1,"dessusGauche");
    geoSetDomainName(2,"milieu"); 
    geoSetDomainName(3,"dessusDroite");    
    geoSetDomainName(4,"petitDroite");   
    geoSetDomainName(5,"dessous"); 
*/
 

    geoMeshWrite("../data/mesh.txt");
          
//
//  -2- Definition du probleme
//
    
    double E   = 211.e9;
    double nu  = 0.3;
    double rho = 7.85e3; 
    double g   = 9.81;
    femProblem* theProblem = femElasticityCreate(theGeometry,E,nu,rho,g,PLANAR_STRAIN);
    //femElasticityAddBoundaryCondition(theProblem,"petitGauche",DIRICHLET_X,0.0);              //DIRICHLET_X,DIRICHLET_Y,NEUMANN_X,NEUMANN_Y,TANGENTIEL,NORMAL
    //femElasticityAddBoundaryCondition(theProblem,"platGauche", DIRICHLET_Y ,0.0);
    femElasticityAddBoundaryCondition(theProblem,"gauche",DIRICHLET_X,0.0);
    femElasticityAddBoundaryCondition(theProblem,"dessus", DIRICHLET_Y ,0.0);
    femElasticityAddBoundaryCondition(theProblem,"droite",DIRICHLET_X,0.0);              //DIRICHLET_X,DIRICHLET_Y,NEUMANN_X,NEUMANN_Y,TANGENTIEL,NORMAL
    //femElasticityAddBoundaryCondition(theProblem,"platDroite", DIRICHLET_Y ,0.0);
    //femElasticityAddBoundaryCondition(theProblem,"petitDroite",NEUMANN_Y,500000.0);
    //femElasticityAddBoundaryCondition(theProblem,"dessous", NEUMANN_X , 900.0);
    femElasticityPrint(theProblem);
    femElasticityWrite(theProblem,"../data/problem.txt");
 

//
//  -3- Champ de la taille de r�f�rence du maillage
//

    double *meshSizeField = malloc(theGeometry->theNodes->nNodes*sizeof(double));
    femNodes *theNodes = theGeometry->theNodes;
    for(int i=0; i < theNodes->nNodes; ++i)
        meshSizeField[i] = theGeometry->geoSize(theNodes->X[i], theNodes->Y[i]);
    double hMin = femMin(meshSizeField,theNodes->nNodes);  
    double hMax = femMax(meshSizeField,theNodes->nNodes);  
    printf(" ==== Global requested h : %14.7e \n",theGeometry->h);
    printf(" ==== Minimum h          : %14.7e \n",hMin);
    printf(" ==== Maximum h          : %14.7e \n",hMax);
    
//
//  -4- Visualisation 
//  
    
    int mode = 1; 
    int domain = 0;
    int freezingButton = FALSE;
    double t, told = 0;
    char theMessage[MAXNAME];
   
 
    GLFWwindow* window = glfemInit("EPL1110 : Project 2022-23 ");
    glfwMakeContextCurrent(window);

    do {
        int w,h;
        glfwGetFramebufferSize(window,&w,&h);
        glfemReshapeWindows(theGeometry->theNodes,w,h);

        t = glfwGetTime();  
        if (glfwGetKey(window,'D') == GLFW_PRESS) { mode = 0;}
        if (glfwGetKey(window,'V') == GLFW_PRESS) { mode = 1;}
        if (glfwGetKey(window,'N') == GLFW_PRESS && freezingButton == FALSE) { domain++; freezingButton = TRUE; told = t;}
        
        if (t-told > 0.5) {freezingButton = FALSE; }
        if (mode == 1) {
            glfemPlotField(theGeometry->theElements,meshSizeField);
            glfemPlotMesh(theGeometry->theElements); 
            sprintf(theMessage, "Number of elements : %d ",theGeometry->theElements->nElem);
            glColor3f(1.0,0.0,0.0); glfemMessage(theMessage); }
        if (mode == 0) {
            domain = domain % theGeometry->nDomains;
            glfemPlotDomain( theGeometry->theDomains[domain]); 
            sprintf(theMessage, "%s : %d ",theGeometry->theDomains[domain]->name,domain);
             glColor3f(1.0,0.0,0.0); glfemMessage(theMessage);  }
            
         glfwSwapBuffers(window);
         glfwPollEvents();
    } while( glfwGetKey(window,GLFW_KEY_ESCAPE) != GLFW_PRESS &&
             glfwWindowShouldClose(window) != 1 );
            
    // Check if the ESC key was pressed or the window was closed

    free(meshSizeField);
    femElasticityFree(theProblem) ; 
    geoFree();
    glfwTerminate(); 
    
    exit(EXIT_SUCCESS);
    return 0;  
}


 
