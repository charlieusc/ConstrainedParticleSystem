/*

  USC/Viterbi/Computer Science
  "Jello Cube" Assignment 1 starter code

*/

#include "physics.h"

//spring rest lenth
const double Rlen = 1.0/7.0;
const double RlenD = Rlen * sqrt(2.0);
const double RlenT = Rlen * sqrt(3.0);
double DiagScal = 1.5;
point xp, xm, yp, ym, zp, zm, ilp;
point statPoint = {0, 0, 0};
double ilPara;



void initPoint(point &a){
    a.x = 0;
    a.y = 0;
    a.z = 0;
}
void initCollisionPoint(point p){
    xp.x = 2;
    xp.y = p.y;
    xp.z = p.z;
    xm.x = -2;
    xm.y = p.y;
    xm.z = p.z;
    
    yp.x = p.x;
    yp.y = 2;
    yp.z = p.z;
    ym.x = p.x;
    ym.y = -2;
    ym.z = p.z;
    
    zp.x = p.x;
    zp.y = p.y;
    zp.z = 2;
    zm.x = p.x;
    zm.y = p.y;
    zm.z = -2;
}


double inlinedPlanePara(world *jello, int i, int j, int k){
    return jello->a * jello->p[i][j][k].x + jello->b * jello->p[i][j][k].y +jello->c * jello->p[i][j][k].z +jello->d;
}
void initInLinedCollisionPoint(world *jello, int i, int j, int k){

    ilp.x = jello->p[i][j][k].x - jello->a/normalizer(jello->a, jello->b, jello->c) * ilPara/normalizer(jello->a, jello->b, jello->c);
    ilp.y = jello->p[i][j][k].y - jello->b/normalizer(jello->a, jello->b, jello->c) * ilPara/normalizer(jello->a, jello->b, jello->c);
    ilp.z = jello->p[i][j][k].z - jello->c/normalizer(jello->a, jello->b, jello->c) * ilPara/normalizer(jello->a, jello->b, jello->c);
}

/* Computes acceleration to every control point of the jello cube, 
   which is in state given by 'jello'.
   Returns result in array 'a'. */

void computeAcceleration(struct world * jello, struct point a[8][8][8])
{
  /* for you to implement ... */
    //F=ma, a=F/m
    int i,j,k;
    point Fh, Fd, Ft;
    double Kh = jello->kElastic;
    double Kd = jello->dElastic;
    double CKh = jello->kCollision;
    double CKd = jello->dCollision;

    for (i=0; i<=7; i++)
        for (j=0; j<=7; j++)
            for (k=0; k<=7; k++)
            {
                initPoint(Fh);
                initPoint(Fd);
                initPoint(Ft);
                //Structural spring force
                //hook & damping
                
                if (i+1<=7) {
                    calculateHookForce(Kh, jello->p[i][j][k], jello->p[i+1][j][k],Rlen , Fh);
                    calculateDampingForce(Kd, jello->p[i][j][k], jello->p[i+1][j][k], jello->v[i][j][k], jello->v[i+1][j][k], Fd);
                }
                if (i-1>=0) {
                    calculateHookForce(Kh, jello->p[i][j][k], jello->p[i-1][j][k],Rlen , Fh);
                    calculateDampingForce(Kd, jello->p[i][j][k], jello->p[i-1][j][k], jello->v[i][j][k], jello->v[i-1][j][k], Fd);
                }
                if (j+1<=7) {
                    calculateHookForce(Kh, jello->p[i][j][k], jello->p[i][j+1][k],Rlen , Fh);
                    calculateDampingForce(Kd, jello->p[i][j][k], jello->p[i][j+1][k], jello->v[i][j][k], jello->v[i][j+1][k], Fd);
                }
                if (j-1>=0) {
                    calculateHookForce(Kh, jello->p[i][j][k], jello->p[i][j-1][k],Rlen , Fh);
                    calculateDampingForce(Kd, jello->p[i][j][k], jello->p[i][j-1][k], jello->v[i][j][k], jello->v[i][j-1][k], Fd);
                }
                if (k+1<=7) {
                    calculateHookForce(Kh, jello->p[i][j][k], jello->p[i][j][k+1],Rlen , Fh);
                    calculateDampingForce(Kd, jello->p[i][j][k], jello->p[i][j][k+1], jello->v[i][j][k], jello->v[i][j][k+1], Fd);
                }
                if (k-1>=0) {
                    calculateHookForce(Kh, jello->p[i][j][k], jello->p[i][j][k-1],Rlen , Fh);
                    calculateDampingForce(Kd, jello->p[i][j][k], jello->p[i][j][k-1], jello->v[i][j][k], jello->v[i][j][k-1], Fd);
                }
                
                //bend spring force
                //hook & damping
                
                if (i+2<=7) {
                    calculateHookForce(Kh, jello->p[i][j][k], jello->p[i+2][j][k],2*Rlen , Fh);
                    calculateDampingForce(Kd, jello->p[i][j][k], jello->p[i+2][j][k], jello->v[i][j][k], jello->v[i+2][j][k], Fd);
                }
                if (i-2>=0) {
                    calculateHookForce(Kh, jello->p[i][j][k], jello->p[i-2][j][k],2*Rlen , Fh);
                    calculateDampingForce(Kd, jello->p[i][j][k], jello->p[i-2][j][k], jello->v[i][j][k], jello->v[i-2][j][k], Fd);
                }
                if (j+2<=7) {
                    calculateHookForce(Kh, jello->p[i][j][k], jello->p[i][j+2][k],2*Rlen , Fh);
                    calculateDampingForce(Kd, jello->p[i][j][k], jello->p[i][j+2][k], jello->v[i][j][k], jello->v[i][j+2][k], Fd);
                }
                if (j-2>=0) {
                    calculateHookForce(Kh, jello->p[i][j][k], jello->p[i][j-2][k],2*Rlen , Fh);
                    calculateDampingForce(Kd, jello->p[i][j][k], jello->p[i][j-2][k], jello->v[i][j][k], jello->v[i][j-2][k], Fd);
                }
                if (k+2<=7) {
                    calculateHookForce(Kh, jello->p[i][j][k], jello->p[i][j][k+2],2*Rlen , Fh);
                    calculateDampingForce(Kd, jello->p[i][j][k], jello->p[i][j][k+2], jello->v[i][j][k], jello->v[i][j][k+2], Fd);
                }
                if (k-2>=0) {
                    calculateHookForce(Kh, jello->p[i][j][k], jello->p[i][j][k-2],2*Rlen , Fh);
                    calculateDampingForce(Kd, jello->p[i][j][k], jello->p[i][j][k-2], jello->v[i][j][k], jello->v[i][j][k-2], Fd);
                }

                //Shear spring force
                //i,j
                if (i+1<=7 && j+1<=7) {
                    calculateHookForce(Kh, jello->p[i][j][k], jello->p[i+1][j+1][k],RlenD , Fh);
                    calculateDampingForce(Kd, jello->p[i][j][k], jello->p[i+1][j+1][k], jello->v[i][j][k], jello->v[i+1][j+1][k], Fd);
                    if (k+1<=7) {
                        calculateHookForce(DiagScal*Kh, jello->p[i][j][k], jello->p[i+1][j+1][k+1],RlenT , Fh);
                        calculateDampingForce(DiagScal*Kd, jello->p[i][j][k], jello->p[i+1][j+1][k+1], jello->v[i][j][k], jello->v[i+1][j+1][k+1], Fd);
                    }
                    if (k-1>=0) {
                        calculateHookForce(DiagScal*Kh, jello->p[i][j][k], jello->p[i+1][j+1][k-1],RlenT , Fh);
                        calculateDampingForce(DiagScal*Kd, jello->p[i][j][k], jello->p[i+1][j+1][k-1], jello->v[i][j][k], jello->v[i+1][j+1][k-1], Fd);
                    }
                }
                if (i+1<=7 && j-1>=0) {
                    calculateHookForce(Kh, jello->p[i][j][k], jello->p[i+1][j-1][k],RlenD , Fh);
                    calculateDampingForce(Kd, jello->p[i][j][k], jello->p[i+1][j-1][k], jello->v[i][j][k], jello->v[i+1][j-1][k], Fd);
                    if (k+1<=7) {
                        calculateHookForce(DiagScal*Kh, jello->p[i][j][k], jello->p[i+1][j-1][k+1],RlenT , Fh);
                        calculateDampingForce(DiagScal*Kd, jello->p[i][j][k], jello->p[i+1][j-1][k+1], jello->v[i][j][k], jello->v[i+1][j-1][k+1], Fd);
                    }
                    if (k-1>=0) {
                        calculateHookForce(DiagScal*Kh, jello->p[i][j][k], jello->p[i+1][j-1][k-1],RlenT , Fh);
                        calculateDampingForce(DiagScal*Kd, jello->p[i][j][k], jello->p[i+1][j-1][k-1], jello->v[i][j][k], jello->v[i+1][j-1][k-1], Fd);
                    }
                }
                if (i-1>=0 && j+1<=7) {
                    calculateHookForce(Kh, jello->p[i][j][k], jello->p[i-1][j+1][k],RlenD , Fh);
                    calculateDampingForce(Kd, jello->p[i][j][k], jello->p[i-1][j+1][k], jello->v[i][j][k], jello->v[i-1][j+1][k], Fd);
                    if (k+1<=7) {
                        calculateHookForce(DiagScal*Kh, jello->p[i][j][k], jello->p[i-1][j+1][k+1],RlenT , Fh);
                        calculateDampingForce(DiagScal*Kd, jello->p[i][j][k], jello->p[i-1][j+1][k+1], jello->v[i][j][k], jello->v[i-1][j+1][k+1], Fd);
                    }
                    if (k-1>=0) {
                        calculateHookForce(DiagScal*Kh, jello->p[i][j][k], jello->p[i-1][j+1][k-1],RlenT , Fh);
                        calculateDampingForce(DiagScal*Kd, jello->p[i][j][k], jello->p[i-1][j+1][k-1], jello->v[i][j][k], jello->v[i-1][j+1][k-1], Fd);
                    }
                }
                if (i-1>=0 && j-1>=0) {
                    calculateHookForce(Kh, jello->p[i][j][k], jello->p[i-1][j-1][k],RlenD , Fh);
                    calculateDampingForce(Kd, jello->p[i][j][k], jello->p[i-1][j-1][k], jello->v[i][j][k], jello->v[i-1][j-1][k], Fd);
                    if (k+1<=7) {
                        calculateHookForce(DiagScal*Kh, jello->p[i][j][k], jello->p[i-1][j-1][k+1],RlenT , Fh);
                        calculateDampingForce(DiagScal*Kd, jello->p[i][j][k], jello->p[i-1][j-1][k+1], jello->v[i][j][k], jello->v[i-1][j-1][k+1], Fd);
                    }
                    if (k-1>=0) {
                        calculateHookForce(DiagScal*Kh, jello->p[i][j][k], jello->p[i-1][j-1][k-1],RlenT , Fh);
                        calculateDampingForce(DiagScal*Kd, jello->p[i][j][k], jello->p[i-1][j-1][k-1], jello->v[i][j][k], jello->v[i-1][j-1][k-1], Fd);
                    }
                    
                }
                //i,k
                if (i+1<=7 && k+1<=7) {
                    calculateHookForce(Kh, jello->p[i][j][k], jello->p[i+1][j][k+1],RlenD , Fh);
                    calculateDampingForce(Kd, jello->p[i][j][k], jello->p[i+1][j][k+1], jello->v[i][j][k], jello->v[i+1][j][k+1], Fd);
                }
                if (i+1<=7 && k-1>=0) {
                    calculateHookForce(Kh, jello->p[i][j][k], jello->p[i+1][j][k-1],RlenD , Fh);
                    calculateDampingForce(Kd, jello->p[i][j][k], jello->p[i+1][j][k-1], jello->v[i][j][k], jello->v[i+1][j][k-1], Fd);
                }
                if (i-1>=0 && k+1<=7) {
                    calculateHookForce(Kh, jello->p[i][j][k], jello->p[i-1][j][k+1],RlenD , Fh);
                    calculateDampingForce(Kd, jello->p[i][j][k], jello->p[i-1][j][k+1], jello->v[i][j][k], jello->v[i-1][j][k+1], Fd);
                }
                if (i-1>=0 && k-1>=0) {
                    calculateHookForce(Kh, jello->p[i][j][k], jello->p[i-1][j][k-1],RlenD , Fh);
                    calculateDampingForce(Kd, jello->p[i][j][k], jello->p[i-1][j][k-1], jello->v[i][j][k], jello->v[i-1][j][k-1], Fd);
                }
                //j,k
                if (j+1<=7 && k+1<=7) {
                    calculateHookForce(Kh, jello->p[i][j][k], jello->p[i][j+1][k+1],RlenD , Fh);
                    calculateDampingForce(Kd, jello->p[i][j][k], jello->p[i][j+1][k+1], jello->v[i][j][k], jello->v[i][j+1][k+1], Fd);
                }
                if (j+1<=7 && k-1>=0) {
                    calculateHookForce(Kh, jello->p[i][j][k], jello->p[i][j+1][k-1],RlenD , Fh);
                    calculateDampingForce(Kd, jello->p[i][j][k], jello->p[i][j+1][k-1], jello->v[i][j][k], jello->v[i][j+1][k-1], Fd);
                }
                if (j-1>=0 && k+1<=7) {
                    calculateHookForce(Kh, jello->p[i][j][k], jello->p[i][j-1][k+1],RlenD , Fh);
                    calculateDampingForce(Kd, jello->p[i][j][k], jello->p[i][j-1][k+1], jello->v[i][j][k], jello->v[i][j-1][k+1], Fd);
                }
                if (j-1>=0 && k-1>=0) {
                    calculateHookForce(Kh, jello->p[i][j][k], jello->p[i][j-1][k-1],RlenD , Fh);
                    calculateDampingForce(Kd, jello->p[i][j][k], jello->p[i][j-1][k-1], jello->v[i][j][k], jello->v[i][j-1][k-1], Fd);
                }
                
                //collision force
                initCollisionPoint(jello->p[i][j][k]);
                if (jello->p[i][j][k].x>2) {
                    calculateHookForce(CKh * (jello->p[i][j][k].x-2), jello->p[i][j][k], xp,0 , Fh);
                    calculateDampingForce(CKd * (jello->p[i][j][k].x-2), jello->p[i][j][k], xp, jello->v[i][j][k], statPoint, Fd);
                }
                if (jello->p[i][j][k].x<-2) {
                    calculateHookForce(CKh * (-2-jello->p[i][j][k].x), jello->p[i][j][k], xm,0 , Fh);
                    calculateDampingForce(CKd * (-2-jello->p[i][j][k].x), jello->p[i][j][k], xm, jello->v[i][j][k], statPoint, Fd);
                }
                if (jello->p[i][j][k].y>2) {
                    calculateHookForce(CKh * (jello->p[i][j][k].y-2) , jello->p[i][j][k], yp,0 , Fh);
                    calculateDampingForce(CKd * (jello->p[i][j][k].y-2), jello->p[i][j][k], yp, jello->v[i][j][k], statPoint, Fd);
                }
                if (jello->p[i][j][k].y<-2) {
                    calculateHookForce(CKh * (-2-jello->p[i][j][k].y), jello->p[i][j][k], ym,0 , Fh);
                    calculateDampingForce(CKd * (-2-jello->p[i][j][k].y), jello->p[i][j][k], ym, jello->v[i][j][k], statPoint, Fd);
                }
                if (jello->p[i][j][k].z>2) {
                    calculateHookForce(CKh * (jello->p[i][j][k].z-2), jello->p[i][j][k], zp,0 , Fh);
                    calculateDampingForce(CKd * (jello->p[i][j][k].z-2), jello->p[i][j][k], zp, jello->v[i][j][k], statPoint, Fd);
                }
                if (jello->p[i][j][k].z<-2) {
                    calculateHookForce(CKh * (-2-jello->p[i][j][k].z), jello->p[i][j][k], zm,0 , Fh);
                    calculateDampingForce(CKd * (-2-jello->p[i][j][k].z), jello->p[i][j][k], zm, jello->v[i][j][k], statPoint, Fd);
                }
                
                //Inlined Plane collision
                /*
                initInLinedCollisionPoint(jello->p[i][j][k]);
                //1.414 *(jello->p[i][j][k].y-jello->p[i][j][k].z)
                if(jello->p[i][j][k].z < jello->p[i][j][k].y){
                    calculateHookForce(CKh, jello->p[i][j][k], ilp,0 , Fh);
                    calculateDampingForce(CKd, jello->p[i][j][k], ilp, jello->v[i][j][k], statPoint, Fd);
                }
                 */
                ilPara = inlinedPlanePara(jello, i, j, k);
                initInLinedCollisionPoint(jello, i, j, k);
                //1.414 *(jello->p[i][j][k].y-jello->p[i][j][k].z)
                if(ilPara > 0){
                    calculateHookForce(CKh, jello->p[i][j][k], ilp,0 , Fh);
                    calculateDampingForce(CKd, jello->p[i][j][k], ilp, jello->v[i][j][k], statPoint, Fd);
                }
                
                
                
                Ft.x =Fh.x + Fd.x + jello->forceField[i * jello->resolution * jello->resolution + j * jello->resolution + k].x;
                Ft.y =Fh.y + Fd.y + jello->forceField[i * jello->resolution * jello->resolution + j * jello->resolution + k].y;
                Ft.z =Fh.z + Fd.z + jello->forceField[i * jello->resolution * jello->resolution + j * jello->resolution + k].z;
                
                
                
                Ft.x+=msForceUp.x+msForceHa.x;
                Ft.y+=msForceUp.y+msForceHa.y;
                Ft.z+=msForceUp.z+msForceHa.z;
                //printf("%d: %f\n",i * jello->resolution * jello->resolution + j * jello->resolution + k,jello->forceField[i * jello->resolution * jello->resolution + j * jello->resolution + k].x);
                /*
                Ft.x = Fh.x + Fd.x;
                Ft.y = Fh.y + Fd.y;
                Ft.z = Fh.z + Fd.z;
                */
                a[i][j][k].x = Ft.x/jello->mass;
                a[i][j][k].y = Ft.y/jello->mass;
                a[i][j][k].z = Ft.z/jello->mass;
                //printf("%d,%d,%d x:%f/n", i,j,k,a[i][j][k].x);
                //printf("%d,%d,%d y:%f/n", i,j,k,a[i][j][k].y);
                //printf("%d,%d,%d z:%f/n", i,j,k,a[i][j][k].z);
                
            }
    msForceUp = {0,0,0};
    msForceHa = {0,0,0};
}

/* performs one step of Euler Integration */
/* as a result, updates the jello structure */
void Euler(struct world * jello)
{
  int i,j,k;
  point a[8][8][8];

  computeAcceleration(jello, a);
  
  for (i=0; i<=7; i++)
    for (j=0; j<=7; j++)
      for (k=0; k<=7; k++)
      {
        jello->p[i][j][k].x += jello->dt * jello->v[i][j][k].x;
        jello->p[i][j][k].y += jello->dt * jello->v[i][j][k].y;
        jello->p[i][j][k].z += jello->dt * jello->v[i][j][k].z;
        jello->v[i][j][k].x += jello->dt * a[i][j][k].x;
        jello->v[i][j][k].y += jello->dt * a[i][j][k].y;
        jello->v[i][j][k].z += jello->dt * a[i][j][k].z;

      }
}

/* performs one step of RK4 Integration */
/* as a result, updates the jello structure */
void RK4(struct world * jello)
{
  point F1p[8][8][8], F1v[8][8][8], 
        F2p[8][8][8], F2v[8][8][8],
        F3p[8][8][8], F3v[8][8][8],
        F4p[8][8][8], F4v[8][8][8];

  point a[8][8][8];


  struct world buffer;

  int i,j,k;

  buffer = *jello; // make a copy of jello

  computeAcceleration(jello, a);

  for (i=0; i<=7; i++)
    for (j=0; j<=7; j++)
      for (k=0; k<=7; k++)
      {
         pMULTIPLY(jello->v[i][j][k],jello->dt,F1p[i][j][k]);
         pMULTIPLY(a[i][j][k],jello->dt,F1v[i][j][k]);
         pMULTIPLY(F1p[i][j][k],0.5,buffer.p[i][j][k]);
         pMULTIPLY(F1v[i][j][k],0.5,buffer.v[i][j][k]);
         pSUM(jello->p[i][j][k],buffer.p[i][j][k],buffer.p[i][j][k]);
         pSUM(jello->v[i][j][k],buffer.v[i][j][k],buffer.v[i][j][k]);
      }

  computeAcceleration(&buffer, a);

  for (i=0; i<=7; i++)
    for (j=0; j<=7; j++)
      for (k=0; k<=7; k++)
      {
         // F2p = dt * buffer.v;
         pMULTIPLY(buffer.v[i][j][k],jello->dt,F2p[i][j][k]);
         // F2v = dt * a(buffer.p,buffer.v);     
         pMULTIPLY(a[i][j][k],jello->dt,F2v[i][j][k]);
         pMULTIPLY(F2p[i][j][k],0.5,buffer.p[i][j][k]);
         pMULTIPLY(F2v[i][j][k],0.5,buffer.v[i][j][k]);
         pSUM(jello->p[i][j][k],buffer.p[i][j][k],buffer.p[i][j][k]);
         pSUM(jello->v[i][j][k],buffer.v[i][j][k],buffer.v[i][j][k]);
      }

  computeAcceleration(&buffer, a);

  for (i=0; i<=7; i++)
    for (j=0; j<=7; j++)
      for (k=0; k<=7; k++)
      {
         // F3p = dt * buffer.v;
         pMULTIPLY(buffer.v[i][j][k],jello->dt,F3p[i][j][k]);
         // F3v = dt * a(buffer.p,buffer.v);     
         pMULTIPLY(a[i][j][k],jello->dt,F3v[i][j][k]);
         pMULTIPLY(F3p[i][j][k],0.5,buffer.p[i][j][k]);
         pMULTIPLY(F3v[i][j][k],0.5,buffer.v[i][j][k]);
         pSUM(jello->p[i][j][k],buffer.p[i][j][k],buffer.p[i][j][k]);
         pSUM(jello->v[i][j][k],buffer.v[i][j][k],buffer.v[i][j][k]);
      }
         
  computeAcceleration(&buffer, a);


  for (i=0; i<=7; i++)
    for (j=0; j<=7; j++)
      for (k=0; k<=7; k++)
      {
         // F3p = dt * buffer.v;
         pMULTIPLY(buffer.v[i][j][k],jello->dt,F4p[i][j][k]);
         // F3v = dt * a(buffer.p,buffer.v);     
         pMULTIPLY(a[i][j][k],jello->dt,F4v[i][j][k]);

         pMULTIPLY(F2p[i][j][k],2,buffer.p[i][j][k]);
         pMULTIPLY(F3p[i][j][k],2,buffer.v[i][j][k]);
         pSUM(buffer.p[i][j][k],buffer.v[i][j][k],buffer.p[i][j][k]);
         pSUM(buffer.p[i][j][k],F1p[i][j][k],buffer.p[i][j][k]);
         pSUM(buffer.p[i][j][k],F4p[i][j][k],buffer.p[i][j][k]);
         pMULTIPLY(buffer.p[i][j][k],1.0 / 6,buffer.p[i][j][k]);
         pSUM(buffer.p[i][j][k],jello->p[i][j][k],jello->p[i][j][k]);

         pMULTIPLY(F2v[i][j][k],2,buffer.p[i][j][k]);
         pMULTIPLY(F3v[i][j][k],2,buffer.v[i][j][k]);
         pSUM(buffer.p[i][j][k],buffer.v[i][j][k],buffer.p[i][j][k]);
         pSUM(buffer.p[i][j][k],F1v[i][j][k],buffer.p[i][j][k]);
         pSUM(buffer.p[i][j][k],F4v[i][j][k],buffer.p[i][j][k]);
         pMULTIPLY(buffer.p[i][j][k],1.0 / 6,buffer.p[i][j][k]);
         pSUM(buffer.p[i][j][k],jello->v[i][j][k],jello->v[i][j][k]);
      }

  return;  
}
