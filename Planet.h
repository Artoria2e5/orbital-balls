#ifndef _KEPLER_H_
#define _KEPLER_H_

#include <Vector.h>
#include <stdint.h>
// meter, kg, sec
#define CONST_G 6.67408e-11
#define M_TAU 6.28318530717958648f

#pragma GCC optimize push
#pragma GCC optimize ("unsafe-math-optimizations, no-signed-zeros, no-trapping-math")

// probably takes ~3200 MCU cycles
inline double atanh(double x) {
   return log((1+x)/(1-x)) * 0.5;
}

// We refuse to handle parabolic stuff.
namespace KepEqn {
  byte trouble = 0;
  static double E_to_M(double E, double ecc) { return E - ecc * sin(E); }
  static double H_to_M(double H, double ecc) { return ecc * sinh(H) - H; }
  static double REAL_to_M(double R, double ecc) { return (ecc <= 1) ? E_to_M(R, ecc) : H_to_M(R, ecc); }

  double M_to_E(double M, double ecc);
  double M_to_H(double M, double ecc);
  static double M_to_REAL(double M, double ecc) { return (ecc <= 1) ? M_to_E(M, ecc) : M_to_H(M, ecc); }
}

enum Attractors {
  SOI_SOL,
  SOI_TERRA,
  SOI_MOON,
};

class Attractor {
  float mass;
  float death_distance;
  Kepler* kepinfo;
};

class Kepler {
public:
  // yes these are KSP savefile names!
  // meters, rad arcs...
  float sma;
  float ecc;
  float lpe;
  float mna;
  int32_t eph;
  enum Attractors parent;

  Kepler() :sma(0), ecc(0), lpe(0), mna(0), eph(0), parent(SOI_SOL) {}

  double meanAngularMotion() {
    return sqrt(mu / abs(sma * sma * sma));
  }

  static float anomalyTrue2Ecce(float true_ano, float ecc) {
    float t = tan(true_ano / 2) / sqrt((1 + ecc) / fabs(1 - ecc));
    return (ecc <= 1.0f) ? 2 * atan(t) : 2 * atanh(t);
  }

  static float anomalyEcce2True(float ecce_ano, float ecc) {
    float y = sqrt(1 + ecc) * (ecc <= 1.0f ? sin(ecce_ano / 2) : sinh(ecce_ano / 2));
    float x = sqrt(fabs(1 - ecc)) * (ecc <= 1.0f ? cos(ecce_ano / 2) : cosh(ecce_ano / 2));
    float tmpret = 2 * atan2(y, x);
  }

  static Kepler fromStateVec(float centmass, StateVec s, int32_t currtime);
  StateVec asStateVecs(int32_t time, float centmass);

  int32_t distanceWhen(int32_t time, float distance, float centmass);
};
#pragma GCC optimize pop

// data time!
Attractor Att_Sol{1.988e30, 5e11, NULL};

// hey this is not a posix time joke please be respectful
Kepler Terra{1.49598e11, 0.0167, 100 / 180 * M_PI, 0, 0, SOI_SOL};
Attractor Att_Terra{5.97237e24, 6371000, &Terra};

Kepler Moon{384399, 0.0549, M_PI / 3, M_TAU / 3, 0, SOI_TERRA};
Attractor Att_Moon{7.342e22, 1737000, &Moon};
#endif
