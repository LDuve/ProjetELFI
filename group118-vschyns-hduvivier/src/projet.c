#include "fem.h"




double geoSize(double x, double y){

    femGeo* theGeometry = geoGetGeometry();
    
    double h = theGeometry->h;
    double x0 = theGeometry->xNotch;
    double y0 = theGeometry->yNotch;
    double r0 = theGeometry->rNotch;
    double h0 = theGeometry->hNotch;
    double d0 = theGeometry->dNotch;
  
    
    double x1 = theGeometry->xHole;
    double y1 = theGeometry->yHole;
    double r1 = theGeometry->rHole;
    double h1 = theGeometry->hHole;
    double d1 = theGeometry->dHole;


    
    double hfinal = h;
        
    return hfinal;
    

}


void geoMeshGenerate() {
 
    int ierr;

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
    
    geoSetSizeCallback(geoSize);   
    gmshModelOccSynchronize(&ierr);  
    gmshOptionSetNumber("Mesh.SaveAll", 1, &ierr);
    gmshModelMeshGenerate(2, &ierr); 
    return;
}