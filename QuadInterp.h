#ifndef __QUADINTERP_H__
#define __QUADINTERP_H__
#include "Vector.h"
#pragma GCC optimize push
#pragma GCC optimize ("unsafe-math-optimizations, no-signed-zeros, no-trapping-math")

// We cannot afford to do a lot of newtons. So we take a limited number
// of points for each conic, and do quadratic interpolation.
class QuadInterp {
public:
  V2Df accel;
  StateVec initial;

  QuadInterp() : accel{0,0}, initial{{0,0}, {0,0}} {}
  QuadInterp(V2Df, StateVec);
  QuadInterp(const StateVec& s1, const StateVec& s2, float deltime) {
    accel = (s2.vel - s1.vel) / deltime;
    initial = s1;
  }

  // The estimation for acceleration is not optimal.
  //
  // better_accel = mu/source.posAtDeltime(timeshift+window/2).hypotSq();
  //
  // The caller should be responsible for that anyway.
  QuadInterp(float timeshift, const QuadInterp& source) {
    accel = source.accel;
    initial.pos = source.posAtDeltime(timeshift);
    initial.vel = source.velAtDeltime(timeshift);
  }

  V2Df posAtDeltime(float deltime) {
    return deltime * (accel * deltime * 0.5 + s1.vel) + s1.pos;
  }
(
  V2Df velAtDeltime(float deltime) {
    return accel * deltime + s1.vel;
  }


  float minimalDistanceAt(const QuadInterp& o, float tmax, float* ret_dist = NULL) {
    // we solve for deriv(square(hypot(pos1 - pos2))) = 0.
    // a x^3 + b x^2 + c x + d = 0
    double x[5] = {0};
#define _pos initial.pos
#define _vel initial.vel
#define _acc accel
    double a = accel.x * (accel.x - o.accel.x) + o.accel.x * (o.accel.x - accel.x)
             + accel.y * (accel.y - o.accel.y) + o.accel.y * (o.accel.y - accel.y);
    double b = 3 * _acc.x * (_vel.x - o._vel.x) + 3 * o._acc.x * (o._vel.x - _vel.x)
             + 3 * _acc.y * (_vel.y - o._vel.y) + 3 * o._acc.y * (o._vel.y - _vel.y);
    double c = 2 * _acc.x * (_pos.x - o._pos.x) + 2 * o._acc.x * (o._pos.x - _pos.x)
             + 2 * _vel.x * (_vel.x - o._vel.x) + 2 * o._vel.x * (o._vel.x - _vel.x)
             + 2 * _acc.y * (_pos.y - o._pos.y) + 2 * o._acc.y * (o._pos.y - _pos.y)
             + 2 * _vel.y * (_vel.y - o._vel.y) + 2 * o._vel.y * (o._vel.y - _vel.y);
    double d = 2 * _pos.x * (_vel.x - o._vel.x) + 2 * o._pos.x * (o._vel.x - _vel.x)
             + 2 * _pos.y * (_vel.y - o._vel.y) + 2 * o._pos.y * (o._vel.y - _vel.y);
#undef _pos _vel _acc
    int n = 0;

    if (fabs(a) >= 1e-8) {
      n = SolveP3(x, b/a, c/a, d/a);
    } else if (fabs(b) >= 1e-8) {
      n = SolveP2(x, c/b, d/b);
    } else {
      n = 1;
      // cx+d = 0
      x[0] = -d/c;
    }

    // Now we check for actual minima...
    float minimumDistance = 1.0/0.0f;
    float minimumTime = -0.0f;
    float t;
    V2Df diff;

    // check for end and start too...
    x[n] = 0.0f;
    n++;
    // let me constrain the windows end to... 1e4 years apparently
    if (tmax > 0.0f && tmax < 1e11) {
      x[n] = tmax;
      n++;
    }
    
    for (int i = 0; i < n; i++) {
      t = x[i];
      if (x < 0.0f || x > tmax)
        continue;
      V2Df diff = posAtDeltime(t) - o.posAtDeltime(t);
      float dist = diff.hypotSq();
      if (dist < minimumDistance) {
        minimumDistance = dist;
        minimumTime = x[i];
      }
    }

    if (ret_dist != NULL)
      *ret_dist = minimumDistance;
    
    return minimumTime;
  }

}

#pragma GCC optimize pop
#endif
