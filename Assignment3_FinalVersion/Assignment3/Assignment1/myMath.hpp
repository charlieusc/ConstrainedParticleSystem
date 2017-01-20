//
//  myMath.hpp
//  Assignment1
//
//  Created by YangJialin on 16/2/11.
//
//

#ifndef myMath_hpp
#define myMath_hpp

#include <stdio.h>
#include "jello.h"
extern point msForceUp,msForceHa;
extern int pointselect;
extern double msforceScale;

point crossProduct(point a, point b);
double normalizer(double a, double b, double c);
double dotP(point a, point b);
void VecMi(point a, point b, point &r);
//calculate hook force between two points
void calculateHookForce(double k, point A, point B, double R, point &F);
//calculate damping force between two points
void calculateDampingForce(double k,point A, point B, point Va, point Vb, point &F);
void calculateDrawingPoint(world *jello, point* drawPointset);
#endif /* myMath_hpp */
