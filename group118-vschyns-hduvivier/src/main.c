#include "stru.h"
#include <stdio.h>


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

    double x1 = 0.0;
    double y1 = 0.0;
    int idRect = gmshModelOccAddRectangle(x1,y1,0.0,l,3*ha,-1,0.0,&ierr); 
    
    double x2 = 0.0;
    double y2 = ha;
    double l2 = 2.0 * l / 5.0;
    int idRect2 = gmshModelOccAddRectangle(x2,y2,0.0,l2,2*ha,-1,0.0,&ierr);    

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
    printf("Ici commence le projet");
    int x = 0;
    double y = 1.2;


    gmshInitialize();
    geoMeshGenerate();
    
      //Ajouter ces lignes
    gmshModelMesh(2);
    gmshOptionSetNumber("Mesh.SurfaceFaces", 1);
    gmshOptionSetNumber("Mesh.Points", 1);
    gmshFltkInitialize();
    gmshFltkRun();
    
    gmshFinalize();
    return 0;
}
