#include "stru.h"
#include <stdio.h>
#include "fem.h"

#include <GLFW/glfw3.h>

void geoMeshGenerate() {
 
    int ierr;
/*
x1: Coordonnée x du coin inférieur gauche du rectangle.
y1: Coordonnée y du coin inférieur gauche du rectangle.
z1: Coordonnée z du coin inférieur gauche du rectangle (zéro par défaut).
x2: Coordonnée x du coin supérieur droit du rectangle.
y2: Coordonnée y du coin supérieur droit du rectangle.
z2: Coordonnée z du coin supérieur droit du rectangle (zéro par défaut).
tag: Étiquette à donner à l'entité créée (zéro par défaut).
makeSolid: Si vrai, crée une surface solide plutôt qu'une surface (faux par défaut).
*/
    double l = 1.;
    double ha = 1.;

    double px1 = 0.0;
    double py1 = 0.0;
    int idRect = gmshModelOccAddRectangle(px1,py1,0.0,l,3*ha,-1,0.0,&ierr); 
    
    double px2 = 0.0;
    double py2 = ha;
    double l2 = 2.0 * l / 5.0;
    int idRect2 = gmshModelOccAddRectangle(px2,py2,0.0,l2,2*ha,-1,0.0,&ierr);    

    double x3 = 3.0 * l / 5.0;
    double y3 = ha;
    double l3 = 2.0 * l / 5.0;
    int idRect3 = gmshModelOccAddRectangle(x3,y3,0.0,l3,2*ha,-1,0.0,&ierr);  
    
    int rect[] = {2,idRect};
    int rect2[] = {2,idRect2};
    int rect3[] = {2,idRect3};

    gmshModelOccCut(rect,2,rect2,2,NULL,NULL,NULL,NULL,NULL,-1,1,1,&ierr);
    gmshModelOccCut(rect,2,rect3,2,NULL,NULL,NULL,NULL,NULL,-1,1,1,&ierr);
    

    gmshModelOccSynchronize(&ierr); 

    return;
}





int main(void)
{
    double Lx = 1.0;
    double Ly = 2.0;

    int ierr;
    geoInitialize();
    femGeo* theGeometry = geoGetGeometry();

    theGeometry->LxPlate =  Lx;
    theGeometry->LyPlate =  Ly;
    theGeometry->xPlate  =  0.0;
    theGeometry->yPlate  =  0.0;
    theGeometry->xHole   =  Lx / 4.0;
    theGeometry->yHole   =  Ly / 4.0;
    theGeometry->rHole   =  Lx / 8.0;
    theGeometry->xNotch  = -Lx / 2.0;
    theGeometry->yNotch  = -Ly / 2.0;
    theGeometry->rNotch  =  Lx / 2.0;

    theGeometry->h       =  Lx * 0.1;    // 0.1 c'est plus joli :-)
    theGeometry->hHole   =  theGeometry->h * 0.2;
    theGeometry->hNotch  =  theGeometry->h * 0.05;
    theGeometry->dHole   =  theGeometry->h * 1.0;
    theGeometry->dNotch  =  theGeometry->h * 4.0;

    geoMeshGenerate();
    geoMeshImport();
    
    geoSetDomainName(0,"Outer Disk");
    geoSetDomainName(1,"Bottom");
    geoSetDomainName(2,"Left");
    geoSetDomainName(3,"Right");
    geoSetDomainName(4,"Top");
    geoSetDomainName(5,"Inner Disk");

    char filename[] = "../data/mesh.txt";
    geoMeshWrite(filename);


    double *meshSizeField = malloc(theGeometry->theNodes->nNodes*sizeof(double));
    femNodes *theNodes = theGeometry->theNodes;
    for(int i=0; i < theNodes->nNodes; ++i)
     meshSizeField[i] = geoSize(theNodes->X[i], theNodes->Y[i]);
    double hMin = femMin(meshSizeField,theNodes->nNodes);  
    double hMax = femMax(meshSizeField,theNodes->nNodes); 

    printf(" ==== Global requested h : %14.7e \n",theGeometry->h);
    printf(" ==== Minimum h          : %14.7e \n",hMin);
    printf(" ==== Maximum h          : %14.7e \n",hMax);

        
    int mode = 1; // Change mode by pressing "j", "k", "l"
    int domain = 0;
    int freezingButton = FALSE;
    double t, told = 0;
    char theMessage[256];
    double pos[2] = {20,460};
 
 
    GLFWwindow* window = glfemInit("EPL1110 : Mesh generation ");
    glfwMakeContextCurrent(window);

    do {
        int w,h;
    
    
        glfwGetFramebufferSize(window,&w,&h);
        glfemReshapeWindows(theGeometry->theNodes,w,h);

        t = glfwGetTime();  
    //    glfemChangeState(&mode, theMeshes->nMesh);
        if (glfwGetKey(window,'D') == GLFW_PRESS) { mode = 0;}
        if (glfwGetKey(window,'V') == GLFW_PRESS) { mode = 1;}
        if (glfwGetKey(window,'N') == GLFW_PRESS && freezingButton == FALSE) { domain++; freezingButton = TRUE; told = t;}

        
        if (t-told > 0.5) {freezingButton = FALSE; }
            
        
        
         
        if (mode == 1) {
            glfemPlotField(theGeometry->theElements, meshSizeField);
            glfemPlotMesh(theGeometry->theElements); 
            sprintf(theMessage, "Number of elements : %d ",theGeometry->theElements->nElem);
 
            
            glColor3f(1.0,0.0,0.0); glfemMessage(theMessage);
 
            
            
            }
        if (mode == 0) {
            domain = domain % theGeometry->nDomains;
            glfemPlotDomain( theGeometry->theDomains[domain]); 
            
            
            
            sprintf(theMessage, "%s : %d ",theGeometry->theDomains[domain]->name,domain);
 
            
            glColor3f(1.0,0.0,0.0); glfemMessage(theMessage);
            }
            
         glfwSwapBuffers(window);
         glfwPollEvents();
    } while( glfwGetKey(window,GLFW_KEY_ESCAPE) != GLFW_PRESS &&
             glfwWindowShouldClose(window) != 1 );
            
    // Check if the ESC key was pressed or the window was closed

    free(meshSizeField);  
    geoFinalize();
    glfwTerminate(); 
    
    exit(EXIT_SUCCESS);
    return 0;  
}
