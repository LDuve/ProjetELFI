#include "stru.h"
#include <stdio.h>



#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include "D:\Documents\EPL\Q6\LEPL1110\Projet\ProjetELFI\group118-vschyns-hduvivier\gmsh\gmsh-4.11.1-Windows64-sdk\include\gmshc.h"


void geoMeshGenerate() {
 
    double l = 5.;
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


    //
//  -2- D�finition de la fonction callback pour la taille de r�f�rence
//      Synchronisation de OpenCascade avec gmsh
//      G�n�ration du maillage (avec l'option Mesh.SaveAll :-)
                  
   
    //geoSetSizeCallback(geoSize);  //donner la carte de pour le maillage //servais a reduire la taille des morceaux
                                  
    gmshModelOccSynchronize(&ierr); //envoyer les informations de géométrie à gmsh      
    gmshOptionSetNumber("Mesh.SaveAll", 1, &ierr);
    gmshModelMeshGenerate(2, &ierr);  //genere le maillage
    //gmshFltkRun(&ierr);  //pour visualiser dans gmsh mais pas utile

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
