//
//  controlPoint.cpp
//  Assignment1
//
//  Created by YangJialin on 16/4/4.
//
//
#include "Eigen/SVD"
#include "controlPoint.hpp"
#include "myMath.hpp"
using Eigen::MatrixXd;

void chainInit(myChain *chain){
    if (pN<1) {
        printf("control point initialization fialed");
        return;
    }
    double b=100;
    
    //chain->ks = pow(b, 2);
    //chain->kd = 2*b;
    chain->ks = 100;
    chain->kd = 10;
    
    chain->dt = 1/100.0;
    for (int i=0; i<pN; i++) {
        chain->m[i] = 1000;
    }
    for (int i=0; i<pN; i++) {
        chain->f[i].x = 0;
        chain->f[i].y = -1.0;
        chain->f[i].z = 0.0;
    }
    for (int i=0; i<pN; i++) {
        chain->v[i].x = 0.0;
        chain->v[i].y = 0.0;
        chain->v[i].z = 0.0;
    }
    chain->p[0].x = 0.0;
    chain->p[0].y = 1.0;
    chain->p[0].z = 0.0;
    
    for (int i=1; i<pN; i++) {
        chain->p[i].x = 0.0;
        chain->p[i].y = 1.0-(2*i)/(double)(pN-1);
        chain->p[i].z = 0.0;
    }
    chain->r = 2.0/(pN-1);
    /*
    for (int i=1; i<(pN+1)/2; i++) {
        chain->p[i].x = 0;
        chain->p[i].y = 1-(2*i)/(float)(pN-1);
        chain->p[i].z = 0;
    }
    
    for (int i=(pN+1)/2; i<pN; i++) {
        chain->p[i].x = 2*(1+i-(pN+1)/2)/(float)(pN-1);
        chain->p[i].y = 0;
        chain->p[i].z = 0;
    }*/
}


void computeAcceleration2(struct myChain * chain, struct point a[pN]){
    //row and col for J
    int row = pN; //number of constrains
    int col = (pN-1)*2; //number of x1, y1,..., xn, yn
    double Radius = chain->r;
    
    //q is the vector to store the position(x, y) of each point.
    MatrixXd q(col, 1);
    //qp is the vector to store the velosity(x, y) of each point.
    MatrixXd qp(col, 1);
    //C is the vector to define the constrian of each point.
    MatrixXd C(row, 1);
    //Cp is the vector to define the constrain of each point's velocity.
    MatrixXd Cp(row, 1);
    //MatrixXd Cp1(row,1);
    //J is the matrix of $C/$q, called the Jacobian of C.
    MatrixXd J(row, col);
    //JT is the transpos of matrix J
    MatrixXd JT(col, row);
    //Jp is the time derivative of the Matrix J.
    MatrixXd Jp(row, col);
    //Q is the matrix of globle force vector.
    MatrixXd Q(col, 1);
    //Qc is the matrix of constrain force.
    MatrixXd Qc(col, 1);
    //Qm is the matrix of mouse force.
    MatrixXd Qm(col, 1);
    //lamda is a vector with the dimention of C.
    MatrixXd lamda(row, 1);
    //define matrix to calculate the lamda
    MatrixXd JJT(row, row);
    MatrixXd right(row, 1);
    //q initialization
    for(int i=0; i<pN-1; i++){
        q(i*2,0) = chain->p[i+1].x;
        q(i*2+1,0) = chain->p[i+1].y;
    }
    //qp initialization
    for(int i=0; i<pN-1; i++){
        qp(i*2,0) = chain->v[i+1].x;
        qp(i*2+1,0) = chain->v[i+1].y;
    }
    
    //C(0,0) = pow((chain->p[1].x-0.0), 2) + pow((chain->p[1].y-1.0), 2) - pow(Radius, 2);
    //the constrain of inner point
    for(int i=0; i<row-1; i++){
        C(i,0) = pow((chain->p[i+1].x-chain->p[i].x),2) + pow((chain->p[i+1].y-chain->p[i].y),2)-pow(Radius,2);
    }
    //the constrain of last point
    C(row-1,0) = pow(chain->p[pN-1].x, 2) + pow(chain->p[pN-1].y, 2) - 1.0;
    if(C(row-1,0)>0.000001||C(row-1,0)<-0.000001)
    printf("%f\n", C(row-1,0));

    //J initialization
    
    for(int j=0; j<col; j++){
        if(j==0){
            J(0,j) = 2*(chain->p[1].x-chain->p[0].x);
        }else if(j==1){
            J(0,j) = 2*(chain->p[1].y-chain->p[0].y);
        }else{
            J(0,j) = 0;
        }
    }
    
    
    //inner row
    for(int i=1; i<row-1; i++){
        for(int j=0; j<col; j++){
            if(j == 2*(i-1)){
                J(i,j) = -2*(chain->p[i+1].x-chain->p[i].x);
            }else if(j == 2*(i-1)+1){
                J(i,j) = -2*(chain->p[i+1].y-chain->p[i].y);
            }else if(j == 2*(i-1)+2){
                J(i,j) = 2*(chain->p[i+1].x-chain->p[i].x);
            }else if(j == 2*(i-1)+3){
                J(i,j) = 2*(chain->p[i+1].y-chain->p[i].y);
            }else{
                J(i,j) = 0;
            }
        }
    }
    
    //last row
    for(int j=0; j<col; j++){
        if(j==col-2){
            J(row-1,j) = 2*chain->p[pN-1].x;
        }else if(j==col-1){
            J(row-1,j) = 2*chain->p[pN-1].y;
        }else{
            J(row-1,j) = 0;
        }
    }
    
    /*------------------------------------------------------------*/
    for(int j=0; j<col; j++){
        if(j==0){
            Jp(0,j) = 2*chain->v[1].x;
        }else if(j==1){
            Jp(0,j) = 2*chain->v[1].y;
        }else{
            Jp(0,j) = 0;
        }
    }
    //Jp initialization
    for(int i=1; i<row-1; i++){
        for(int j=0; j<col; j++){
            if(j == 2*(i-1)){
                Jp(i,j) = -2*(chain->v[i+1].x-chain->v[i].x);
            }else if(j == 2*(i-1)+1){
                Jp(i,j) = -2*(chain->v[i+1].y-chain->v[i].y);
            }else if(j == 2*(i-1)+2){
                Jp(i,j) = 2*(chain->v[i+1].x-chain->v[i].x);
            }else if(j == 2*(i-1)+3){
                Jp(i,j) = 2*(chain->v[i+1].y-chain->v[i].y);
            }else{
                Jp(i,j) = 0;
            }
        }
    }
    
    //last row
    for(int j=0; j<col; j++){
        if(j==col-2){
            Jp(row-1,j) = 2*chain->v[pN-1].x;
        }else if(j==col-1){
            Jp(row-1,j) = 2*chain->v[pN-1].y;
        }else{
            Jp(row-1,j) = 0;
        }
    }
    
    Cp = J*qp;
    
    //initializing Q
    for(int i=0; i<pN-1; i++){
        Q(i*2, 0) = chain->f[i+1].x;
        Q(i*2+1, 0) = chain->f[i+1].y;
    }
    
    //initializing Qm
    for(int i=0; i<pN-1; i++){
        Qm(i*2, 0) = 0;
        Qm(i*2+1, 0) = 0;
        if(g_iLeftMouseButton){
            if(i==pointselect-1){
                if(pointselect<pN){
                    Qm(i*2,0) = msForceHa.x+msForceUp.x;
                    Qm(i*2+1,0) = msForceHa.y+msForceUp.y;
                }
            }else if(pointselect == 0){
                Qm((pN-2)*2,0) = msForceHa.x+msForceUp.x;
                Qm((pN-2)*2+1,0) = msForceHa.y+msForceUp.y;
            }
        }
    }
    
    //calculate lamda
    JT = J.transpose();
    JJT = J*JT;
    right = -(Jp*qp) - (J*(Q+Qm)) - (C*chain->ks) - (Cp*chain->kd);
    Eigen::JacobiSVD<MatrixXd> svd(JJT, Eigen::ComputeThinU | Eigen::ComputeThinV);
    svd.setThreshold(1E-6);
    lamda = svd.solve(right);
    
    //calculate constrain force
    Qc = JT*lamda;
    a[0].x=0;
    a[0].y=0;
    a[0].z=0;
    for (int i=1; i<pN; i++) {
        a[i].x = Q((i-1)*2,0) + Qc((i-1)*2,0) + Qm((i-1)*2,0);
        a[i].y = Q((i-1)*2+1,0) + Qc((i-1)*2+1,0) + Qm((i-1)*2+1,0);
        a[i].z = 0;
    }
    
    
}


void myEuler(struct myChain * chain)
{
    point a[pN];
    
    computeAcceleration2(chain, a);
    
    for (int i=1; i<pN; i++)
        
    {
        chain->v[i].z = 0.0;
        chain->v[i].x += chain->dt * a[i].x;
        chain->v[i].y += chain->dt * a[i].y;

        chain->p[i].z = 0.0;
        chain->p[i].x += chain->dt * chain->v[i].x;
        chain->p[i].y += chain->dt * chain->v[i].y;
    }

}


void myRK4(struct myChain * chain)
{
    point F1p[pN], F1v[pN],
    F2p[pN], F2v[pN],
    F3p[pN], F3v[pN],
    F4p[pN], F4v[pN];
    
    point a[pN];
    
    
    struct myChain buffer;
    
    int i;
    
    buffer = *chain; // make a copy of jello
    
    computeAcceleration2(chain, a);
    
    for (i=0; i<pN; i++)
            {
                pMULTIPLY(chain->v[i],chain->dt,F1p[i]);
                pMULTIPLY(a[i],chain->dt,F1v[i]);
                pMULTIPLY(F1p[i],0.5,buffer.p[i]);
                pMULTIPLY(F1v[i],0.5,buffer.v[i]);
                pSUM(chain->p[i],buffer.p[i],buffer.p[i]);
                pSUM(chain->v[i],buffer.v[i],buffer.v[i]);
            }
    
    computeAcceleration2(&buffer, a);
    
    for (i=0; i<pN; i++)
            {
                // F2p = dt * buffer.v;
                pMULTIPLY(buffer.v[i],chain->dt,F2p[i]);
                // F2v = dt * a(buffer.p,buffer.v);
                pMULTIPLY(a[i],chain->dt,F2v[i]);
                pMULTIPLY(F2p[i],0.5,buffer.p[i]);
                pMULTIPLY(F2v[i],0.5,buffer.v[i]);
                pSUM(chain->p[i],buffer.p[i],buffer.p[i]);
                pSUM(chain->v[i],buffer.v[i],buffer.v[i]);
            }
    
    computeAcceleration2(&buffer, a);
    
    for (i=0; i<pN; i++)
            {
                // F3p = dt * buffer.v;
                pMULTIPLY(buffer.v[i],chain->dt,F3p[i]);
                // F3v = dt * a(buffer.p,buffer.v);
                pMULTIPLY(a[i],chain->dt,F3v[i]);
                pMULTIPLY(F3p[i],0.5,buffer.p[i]);
                pMULTIPLY(F3v[i],0.5,buffer.v[i]);
                pSUM(chain->p[i],buffer.p[i],buffer.p[i]);
                pSUM(chain->v[i],buffer.v[i],buffer.v[i]);
            }
    
    computeAcceleration2(&buffer, a);
    
    
    for (i=0; i<pN; i++)
            {
                // F3p = dt * buffer.v;
                pMULTIPLY(buffer.v[i],chain->dt,F4p[i]);
                // F3v = dt * a(buffer.p,buffer.v);
                pMULTIPLY(a[i],chain->dt,F4v[i]);
                
                pMULTIPLY(F2p[i],2,buffer.p[i]);
                pMULTIPLY(F3p[i],2,buffer.v[i]);
                pSUM(buffer.p[i],buffer.v[i],buffer.p[i]);
                pSUM(buffer.p[i],F1p[i],buffer.p[i]);
                pSUM(buffer.p[i],F4p[i],buffer.p[i]);
                pMULTIPLY(buffer.p[i],1.0 / 6,buffer.p[i]);
                pSUM(buffer.p[i],chain->p[i],chain->p[i]);
                
                pMULTIPLY(F2v[i],2,buffer.p[i]);
                pMULTIPLY(F3v[i],2,buffer.v[i]);
                pSUM(buffer.p[i],buffer.v[i],buffer.p[i]);
                pSUM(buffer.p[i],F1v[i],buffer.p[i]);
                pSUM(buffer.p[i],F4v[i],buffer.p[i]);
                pMULTIPLY(buffer.p[i],1.0 / 6,buffer.p[i]);
                pSUM(buffer.p[i],chain->v[i],chain->v[i]);
            }
    
    return;  
}
