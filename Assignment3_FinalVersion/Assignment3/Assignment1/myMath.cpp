//
//  myMath.cpp
//  Assignment1
//
//  Created by YangJialin on 16/2/11.
//
//

#include "myMath.hpp"

point msForceUp = {0,0,0};
point msForceHa = {0,0,0};
int pointselect = 1;
double msforceScale = 0.1;

point crossProduct(point a, point b){
    point c;
    c.x = a.y * b.z - a.z * b.y;
    c.y = a.z * b.x - a.x * b.z;
    c.z = a.x * b.y - a.y * b.x;
    
    return c;
}

double normalizer(double a, double b, double c){
    return sqrt(a*a + b*b + c*c);
}

double dotP(point a, point b){
    return a.x * b.x + a.y * b.y + a.z * b.z;
}

void VecMi(point a, point b, point &r){
    r.x = a.x - b.x;
    r.y = a.y - b.y;
    r.z = a.z - b.z;
}

void calculateHookForce(double k, point A, point B, double R, point &F){
    
    //F = -k * (|A-B| - R)*(A-B)/|A-B|
    
    point L;
    VecMi(A, B, L);
    double La = sqrt(L.x * L.x + L.y * L.y + L.z * L.z);
    F.x += -k * (La - R) * L.x / La;
    F.y += -k * (La - R) * L.y / La;
    F.z += -k * (La - R) * L.z / La;
    
}

void calculateDampingForce(double k,point A, point B, point Va, point Vb, point &F){
    //F = -k * ((Va-Vb)dot(B-A)/|B-A|) * (B-A)/|B-A|
    
    point L,Vab;
    VecMi(A, B, L);
    VecMi(Va, Vb, Vab);
    double La = sqrt(L.x * L.x + L.y * L.y + L.z * L.z);
    double dot = dotP(Vab, L);
    F.x += -k * dot / La * L.x / La;
    F.y += -k * dot / La * L.y / La;
    F.z += -k * dot / La * L.z / La;
    
}

void calculateDrawingPoint(world *jello, point *drawPointset){
    int counter = 0;
    double dx, dy, dz;
    for (double j = -2; j < 2; j+=0.1) {
        for (double k = -2; k < 2; k+=0.1) {
            if (jello->a != 0) {
                dx = (-jello->b * j - jello->c * k - jello->d)/jello->a;
                drawPointset[counter]={dx, j, k};
            }else if(jello->b != 0){
                dy = (- jello->c * k - jello->d)/jello->b;
                drawPointset[counter]={j, dy, k};
            }else if(jello->c != 0){
                dz = (- jello->d)/jello->c;
                drawPointset[counter]={j, k, dz};
            }else{
                drawPointset[counter]={10, 10, 10};
            }
            counter++;
            
        }
    }
    return;
}
