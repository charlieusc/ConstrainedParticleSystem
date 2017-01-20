/*

  USC/Viterbi/Computer Science
  "Jello Cube" Assignment 1 starter code

*/


#ifndef _SHOWCUBE_H_
#define _SHOWCUBE_H_
#include "myMath.hpp"
#include "controlPoint.hpp"

void showCube(struct world * jello);
void showInlinedPlane(world *jello, point *mySet);
void showBoundingBox();
void showAxis();
void DrawCircle(float cx, float cy, float r, int num_segments);
void DrawPoint(myChain *chain);
void DrawStick(myChain *chain);
void DrawChain(myChain *chain);

#endif
