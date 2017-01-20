//
//  controlPoint.hpp
//  Assignment1
//
//  Created by YangJialin on 16/4/4.
//
//

#ifndef controlPoint_hpp
#define controlPoint_hpp

#include <stdio.h>
#include "jello.h"
#endif /* controlPoint_hpp */
//define control points


void chainInit(myChain *chain);
void computeAcceleration2(struct myChain * chain, struct point a[pN]);
void myEuler(struct myChain * chain);
void myRK4(struct myChain * chain);