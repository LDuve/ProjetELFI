#include "fem.h"

//
// Ici, vous pouvez définir votre géométrie :-)
//  (1) Raffiner intelligemment.... (yes )
//  (2) Construire la geometrie avec OpenCascade 
//  (3) Construire la geometrie avec les outils de GMSH
//  (4) Obtenir la geometrie en lisant un fichier .geo de GMSH


double geoSize(double x, double y) {

    femGeo* theGeometry = geoGetGeometry();
    return theGeometry->h * (1.0 - 0.5 * x);
}



//contour avec les lignes  8 points
void geoMeshGenerateGeoProjet() {

    femGeo* theGeometry = geoGetGeometry();
    geoSetSizeCallback(geoSize); 
   
    int ierr;
    double LxBar = theGeometry->LxBar ;
    double LyBar = theGeometry->LyBar;
    double Lx = theGeometry->LxPlate;  
    double Ly = theGeometry->LyPlate; 
    double hTot = Ly+LyBar;
    double lc = theGeometry->h;     

    int p1 = gmshModelGeoAddPoint( -2.5*Lx, 0., 0., lc, 1, &ierr);
    int p2 = gmshModelGeoAddPoint(-2.5*Lx,  Ly, 0., lc, 2, &ierr);
    int p3 = gmshModelGeoAddPoint( -LxBar/2, Ly, 0., lc, 3, &ierr);
    int p4 = gmshModelGeoAddPoint(-LxBar/2, hTot, 0., lc, 4, &ierr);
    int p5 = gmshModelGeoAddPoint(LxBar/2,  hTot  , 0., lc, 5, &ierr);
    int p6 = gmshModelGeoAddPoint(LxBar/2,  Ly  , 0., lc, 6, &ierr);
    int p7 = gmshModelGeoAddPoint(2.5*Lx, Ly  , 0., lc, 7, &ierr);
    int p8 = gmshModelGeoAddPoint( 2.5*Lx, 0.  , 0., lc, 8, &ierr);


    int l1 = gmshModelGeoAddLine(p1, p2, 1, &ierr);
    int l2 = gmshModelGeoAddLine(p2, p3, 2, &ierr);
    int l3 = gmshModelGeoAddLine(p3, p4, 3, &ierr);
    int l4 = gmshModelGeoAddLine(p4, p5, 4, &ierr);
    int l5 = gmshModelGeoAddLine(p5, p6, 5, &ierr);
    int l6 = gmshModelGeoAddLine(p6, p7, 6, &ierr);
    int l7 = gmshModelGeoAddLine(p7, p8, 7, &ierr);
    int l8 = gmshModelGeoAddLine(p8, p1, 8, &ierr);

    
    int lTags[] = {l1, l2, l3, l4, l5, l6, l7, l8 };  
    int c1[] = {1};
    c1[0] = gmshModelGeoAddCurveLoop(lTags, 8, 1, 0, &ierr);  
    int s1 = gmshModelGeoAddPlaneSurface(c1, 1, 1, &ierr);

    gmshModelGeoSynchronize(&ierr);


    if (theGeometry->elementType == FEM_QUAD) {
        gmshOptionSetNumber("Mesh.SaveAll",1,&ierr);
        gmshOptionSetNumber("Mesh.RecombineAll",1,&ierr);
        gmshOptionSetNumber("Mesh.Algorithm",8,&ierr);  
        gmshOptionSetNumber("Mesh.RecombinationAlgorithm",1.0,&ierr); 
        gmshModelGeoMeshSetRecombine(2,1,45,&ierr);  
        gmshModelMeshGenerate(2,&ierr);  }
  
    if (theGeometry->elementType == FEM_TRIANGLE) {
        gmshOptionSetNumber("Mesh.SaveAll",1,&ierr);
        gmshModelMeshGenerate(2,&ierr);  }

    gmshFltkRun(&ierr);
}


//juste la barre           4 points
void geoMeshGenerateBar() {
    
    femGeo* theGeometry = geoGetGeometry();
    geoSetSizeCallback(geoSize);   


    int ierr;
    double LxBar = theGeometry->LxBar ;  
    double h = theGeometry->LyBar + theGeometry->LyPlate; 
    double lc = theGeometry->h;     



    int p1 = gmshModelGeoAddPoint(-LxBar/2, -h/2, 0., lc, 1, &ierr);
    int p2 = gmshModelGeoAddPoint(-LxBar/2,  h/2, 0., lc, 2, &ierr);
    int p3 = gmshModelGeoAddPoint( LxBar/2,  h/2, 0., lc, 3, &ierr);
    int p4 = gmshModelGeoAddPoint( LxBar/2, -h/2, 0., lc, 4, &ierr);
    
    


    int l1 = gmshModelGeoAddLine(p1, p2, 1, &ierr);
    int l2 = gmshModelGeoAddLine(p2, p3, 2, &ierr);
    int l3 = gmshModelGeoAddLine(p3, p4, 3, &ierr);
    int l4 = gmshModelGeoAddLine(p4, p1, 4, &ierr);


    int lTags[] = {l1, l2, l3, l4};  
    int c1[] = {1};
    c1[0] = gmshModelGeoAddCurveLoop(lTags, 4, 1, 0, &ierr);  
    int s1 = gmshModelGeoAddPlaneSurface(c1, 1, 1, &ierr);
    gmshModelGeoSynchronize(&ierr);


    if (theGeometry->elementType == FEM_QUAD) {
        gmshOptionSetNumber("Mesh.SaveAll",1,&ierr);
        gmshOptionSetNumber("Mesh.RecombineAll",1,&ierr);
        gmshOptionSetNumber("Mesh.Algorithm",4,&ierr);  
        gmshOptionSetNumber("Mesh.RecombinationAlgorithm",1.0,&ierr); 
        gmshModelGeoMeshSetRecombine(2,1,45,&ierr);  
        gmshModelMeshGenerate(2,&ierr);  }
  
    if (theGeometry->elementType == FEM_TRIANGLE) {
        gmshOptionSetNumber("Mesh.SaveAll",1,&ierr);
        gmshModelMeshGenerate(2,&ierr);  }

    //gmshFltkRun(&ierr);
}


//juste le projecteur      6 points
void geoMeshGenerateProjecteur() {

    femGeo* theGeometry = geoGetGeometry();
    geoSetSizeCallback(geoSize); 
   
    int ierr;
    double Lx = theGeometry->LxPlate;  
    double LxBar = theGeometry->LxBar ;  
    double Ly = theGeometry->LyPlate; 
    double lc = theGeometry->h;     

    int p1 = gmshModelGeoAddPoint( -2.5*Lx, 0., 0., lc, 1, &ierr);
    int p2 = gmshModelGeoAddPoint( -2.5*Lx, Ly, 0., lc, 2, &ierr);
    int p3 = gmshModelGeoAddPoint(-LxBar/2, Ly, 0., lc, 3, &ierr);
    int p4 = gmshModelGeoAddPoint( LxBar/2, Ly, 0., lc, 4, &ierr);
    int p5 = gmshModelGeoAddPoint(2.5*Lx,   Ly, 0., lc, 5, &ierr);
    int p6 = gmshModelGeoAddPoint(2.5*Lx,   0., 0., lc, 6, &ierr);


    int l1 = gmshModelGeoAddLine(p1, p2, 1, &ierr);
    int l2 = gmshModelGeoAddLine(p2, p3, 2, &ierr);
    int l3 = gmshModelGeoAddLine(p3, p4, 3, &ierr);
    int l4 = gmshModelGeoAddLine(p4, p5, 4, &ierr);
    int l5 = gmshModelGeoAddLine(p5, p6, 5, &ierr);
    int l6 = gmshModelGeoAddLine(p6, p1, 6, &ierr);



    int lTags[] = {l1, l2, l3, l4, l5, l6}; 
    int c1[] = {1};
    c1[0] = gmshModelGeoAddCurveLoop(lTags, 6, 1, 0, &ierr);  
    int s1 = gmshModelGeoAddPlaneSurface(c1, 1, 1, &ierr);
    gmshModelGeoSynchronize(&ierr);


    if (theGeometry->elementType == FEM_QUAD) {
        gmshOptionSetNumber("Mesh.SaveAll",1,&ierr);
        gmshOptionSetNumber("Mesh.RecombineAll",1,&ierr);
        gmshOptionSetNumber("Mesh.Algorithm",6,&ierr);  
        gmshOptionSetNumber("Mesh.RecombinationAlgorithm",1.0,&ierr); 
        gmshModelGeoMeshSetRecombine(2,1,45,&ierr);  
        gmshModelMeshGenerate(2,&ierr);  }
  
    if (theGeometry->elementType == FEM_TRIANGLE) {
        gmshOptionSetNumber("Mesh.SaveAll",1,&ierr);
        gmshModelMeshGenerate(2,&ierr);  }
    //gmshFltkRun(&ierr);
}



void geoMeshGenerateGeoFile(const char *filename){
    femGeo* theGeometry = geoGetGeometry();
    int ierr;
    gmshOpen(filename, &ierr); ErrorGmsh(ierr);
    if (theGeometry->elementType == FEM_QUAD) {
        gmshOptionSetNumber("Mesh.SaveAll",1,&ierr);
        gmshOptionSetNumber("Mesh.RecombineAll",1,&ierr);
        gmshOptionSetNumber("Mesh.Algorithm",8,&ierr);   
        gmshOptionSetNumber("Mesh.RecombinationAlgorithm",1.0,&ierr); 
        gmshModelGeoMeshSetRecombine(2,1,45,&ierr); 
        gmshModelMeshGenerate(2,&ierr);  }
 
    if (theGeometry->elementType == FEM_TRIANGLE) {
        gmshOptionSetNumber("Mesh.SaveAll",1,&ierr);
        gmshModelMeshGenerate(2,&ierr);  }
    return;
}
