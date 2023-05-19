 #include "fem.h"


double geoSize(double x, double y){

    femGeo* theGeometry = geoGetGeometry();
    
    double h = theGeometry->h; //Taille de référence globale


    return h;
}
//enlever les carré
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

    gmshFltkRun(&ierr);

    return;
}

//contour avec les lignes
void geoMeshGenerateGeoProjet() {

    femGeo* theGeometry = geoGetGeometry();
    geoSetSizeCallback(geoSize); 
   
    int ierr;
    double l = theGeometry->LxPlate;  // w=1
    double h = theGeometry->LyPlate; // h=1
    double lc = theGeometry->h;     //0.05

    int p1 = gmshModelGeoAddPoint(-2.5*l, 0, 0., lc, 1, &ierr);
    int p2 = gmshModelGeoAddPoint( -2.5*l, h, 0., lc, 2, &ierr);
    int p3 = gmshModelGeoAddPoint( -0.5*l ,  h, 0., lc, 3, &ierr);
    int p4 = gmshModelGeoAddPoint(-0.5*l,  3*h, 0., lc, 4, &ierr);
    int p5 = gmshModelGeoAddPoint(0.5*l,    3*h, 0., lc, 5, &ierr);
    int p6 = gmshModelGeoAddPoint(0.5*l,      h, 0., lc, 6, &ierr);
    int p7 = gmshModelGeoAddPoint(2.5*l,     h, 0., lc, 7, &ierr);
    int p8 = gmshModelGeoAddPoint(2.5*l,   0, 0., lc, 8, &ierr);


    int l1 = gmshModelGeoAddLine(p1, p2, 1, &ierr);
    int l2 = gmshModelGeoAddLine(p2, p3, 2, &ierr);
    int l3 = gmshModelGeoAddLine(p3, p4, 3, &ierr);
    int l4 = gmshModelGeoAddLine(p4, p5, 4, &ierr);
    int l5 = gmshModelGeoAddLine(p5, p6, 5, &ierr);
    int l6 = gmshModelGeoAddLine(p6, p7, 6, &ierr);
    int l7 = gmshModelGeoAddLine(p7, p8, 7, &ierr);
    int l8 = gmshModelGeoAddLine(p8, p1, 8, &ierr);

    //int l9 = gmshModelGeoAddLine(p3, p6, 9, &ierr);
    
    int lTags[] = {l1, l2, l3, l4, l5, l6, l7, l8 };  //,l9
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

    //gmshFltkRun(&ierr);
}


//deux carré dans mm fct essai
void geoMeshGenerateBarProjet() {

    femGeo* theGeometry = geoGetGeometry();
    geoSetSizeCallback(geoSize); 
   
    int ierr;
    double l = theGeometry->LxPlate;  // w=1
    double h = theGeometry->LyPlate; // h=1
    double lc = theGeometry->h;     //0.05

    int p1 = gmshModelGeoAddPoint(-2.5*l, 0, 0., lc, 1, &ierr);
    int p2 = gmshModelGeoAddPoint( -2.5*l, h, 0., lc, 2, &ierr);
    int p3 = gmshModelGeoAddPoint( -0.5*l ,  h, 0., lc, 3, &ierr);
    int p4 = gmshModelGeoAddPoint(-0.5*l,  3*h, 0., lc, 4, &ierr);
    int p5 = gmshModelGeoAddPoint(0.5*l,    3*h, 0., lc, 5, &ierr);
    int p6 = gmshModelGeoAddPoint(0.5*l,      h, 0., lc, 6, &ierr);
    int p7 = gmshModelGeoAddPoint(2.5*l,     h, 0., lc, 7, &ierr);
    int p8 = gmshModelGeoAddPoint(2.5*l,   0, 0., lc, 8, &ierr);


    int l1 = gmshModelGeoAddLine(p1, p2, 1, &ierr);
    int l2 = gmshModelGeoAddLine(p2, p3, 2, &ierr);
    int l3 = gmshModelGeoAddLine(p3, p4, 3, &ierr);
    int l4 = gmshModelGeoAddLine(p4, p5, 4, &ierr);
    int l5 = gmshModelGeoAddLine(p5, p6, 5, &ierr);
    int l6 = gmshModelGeoAddLine(p6, p7, 6, &ierr);
    int l7 = gmshModelGeoAddLine(p7, p8, 7, &ierr);
    int l8 = gmshModelGeoAddLine(p8, p1, 8, &ierr);

    int l9 = gmshModelGeoAddLine(p3, p6, 9, &ierr);
    
    int lTags[] = {l1, l2, l3, l4, l5, l6, l7, l8 ,l9};  //
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

    //gmshFltkRun(&ierr);
}

//juste la barre
void geoMeshGenerateBar() {
    
    femGeo* theGeometry = geoGetGeometry();
    geoSetSizeCallback(geoSize);   


    int ierr;
    double w = theGeometry->LxBar ;  // w=1
    double h = theGeometry->LyBar + theGeometry->LyPlate; // h=1
    double r = w/4;
    double lc = theGeometry->h;     //0.05


    int p1 = gmshModelGeoAddPoint(-w/2, -h/2, 0., lc, 1, &ierr);
    int p2 = gmshModelGeoAddPoint(-w/2,  h/2, 0., lc, 4, &ierr);
    int p3 = gmshModelGeoAddPoint( w/2,  h/2, 0., lc, 3, &ierr);
    int p4 = gmshModelGeoAddPoint( w/2, -h/2, 0., lc, 2, &ierr);
    
    


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
//juste le projecteur
void geoMeshGenerateProjecteur() {

    femGeo* theGeometry = geoGetGeometry();
    geoSetSizeCallback(geoSize); 
   
    int ierr;
    double l = theGeometry->LxPlate;  // w=1
    double w = theGeometry->LxBar ;  // w=1
    double h = theGeometry->LyPlate; // h=1
    double lc = theGeometry->h;     //0.05

    int p1 = gmshModelGeoAddPoint(-2.5*l, 0, 0., lc, 1, &ierr);
    int p2 = gmshModelGeoAddPoint( -2.5*l, h, 0., lc, 2, &ierr);
    int p3 = gmshModelGeoAddPoint(-w/2,  h/2, 0., lc, 4, &ierr);
    int p4 = gmshModelGeoAddPoint( w/2,  h/2, 0., lc, 3, &ierr);
    int p7 = gmshModelGeoAddPoint(2.5*l,     h, 0., lc, 7, &ierr);
    int p8 = gmshModelGeoAddPoint(2.5*l,   0, 0., lc, 8, &ierr);


    int l1 = gmshModelGeoAddLine(p1, p2, 1, &ierr);
    int l2 = gmshModelGeoAddLine(p2, p3, 2, &ierr);
    int l3 = gmshModelGeoAddLine(p3, p4, 2, &ierr);
    int l4 = gmshModelGeoAddLine(p4, p8, 3, &ierr);
    int l5 = gmshModelGeoAddLine(p8, p1, 4, &ierr);



    int lTags[] = {l1, l2, l3, l4,l5}; 
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






double *femElasticitySolve(femProblem *theProblem)
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
        
        if (theConstrainedNodes[i] != -1  && (theProblem->conditions[theConstrainedNodes[i]]->type == DIRICHLET_X ||  theProblem->conditions[theConstrainedNodes[i]]->type == DIRICHLET_Y)) { //correspondant au nœud actuel que nous examinons. si different de -1 alors il y a une contraite 
            //printf("RENTRE  Dirichlet!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n");
            double value = theProblem->conditions[theConstrainedNodes[i]]->value;
            femFullSystemConstrain(theSystem,i,value); }
        
        
        if (theConstrainedNodes[i] != -1 && (theProblem->conditions[theConstrainedNodes[i]]->type == NEUMANN_X ||  theProblem->conditions[theConstrainedNodes[i]]->type == NEUMANN_Y)) { 
            //printf("RENTRE Neumannnn!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n");
            double GoodValue = theProblem->conditions[theConstrainedNodes[i]]->domain->elemUnique->value ;
            femFullSystemConstrainNeumann(theSystem,i,GoodValue);
            
            }
        }
                 
    return femFullSystemEliminateFrontal(theSystem);
}

double *femElasticitySolveSym(femProblem *theProblem)
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
    
    //printMatrices(theSystem);
    //printBoundaryCondition(theProblem, theProblem->conditions[0]); 
    //printBoundaryCondition(theProblem ,theProblem->conditions[1]);
    //printBoundaryCondition(theProblem ,theProblem->conditions[2]);

    int *theConstrainedNodes = theProblem->constrainedNodes;
    /*
    for(int i=0; i < theSystem->size; i++) {
        if (theConstrainedNodes[i] != -1){printf("node === %d \n" , i);}
        printf("theConstrainedNodes[i] == %d \n", theConstrainedNodes[i]);}
    */  
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
                 
    return femFullSystemEliminate(theSystem);
}



