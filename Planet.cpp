#include "Planet.h"

#pragma GCC optimize push
#pragma GCC optimize ("unsafe-math-optimizations, no-signed-zeros, no-trapping-math")

namespace KepEqn {
  // 978-3-662-11187-1, 4 ed.
  double M_to_E(double M, double ecc) {
    double E = e < 0.8 ? M : M_PI;
    double f;

    // Newton
    // Probably ~5000 cycle/iter
    int i = 0;
    do {
      E -= (f = E - ecc * sin(E) - M) / (1.0 - ecc * cos(E));
      if (++i > 15) {
        trouble |= 0x01;
        break;
      }
    } while (fabs(f) > DBL_EPSILON * 1e4)
    return E;
  }

  double M_to_H(double M, double ecc) {
    double H = ((M < 0) ? -1 : 1) * log(2.0 + fabs(M) / ecc + 1.8);
    double f;

    // Newton
    // ~7000 cycles/iter
    int i = 0;
    do {
      E -= (f = ecc * sinh(H) - H - M) / (ecc * cosh(E) - 1.0);
      if (++i > 15) {
        trouble |= 0x02;
        break;
      }
    } while (fabs(f) > DBL_EPSILON * 1e4 * (1.0 + fabs(H + M)))
    return E;
  }
}

Kepler Kepler::fromStateVec(StateVec s, int32_t currtime, float centmass) {
  float mu = CONST_G * centmass;
  float h = s.pos.cross(s.vel);
  V2Df e = s.vel.cross(h) / mu - s.pos / s.pos.length();
  V2Df n{-h.y, h.x};
  
  float true_ano = arccos(e.dot(s.pos) / e.length() / s.pos.length());
  if (s.pos.dot(s.vel) < 0) true_ano = M_TAU - true_ano;

  float ecce_ano = Kepler::anomalyTrue2Ecce(true_ano, e.length());

  float argu_peri = acos(n.dot(e) / n.length() / e.length());
  float mna = KepEqn::REAL_to_M(ecce_ano, e.length());
  float sma = 1 / (2 / s.pos.length() - square(s.vel.length()) / mu);

  float ecc = e.length();
  if (ecc == 1.0f) ecc -= FLT_EPSILON;

  return Kepler{sma, ecc, argu_peri, mna, currtime};
}

StateVec Kepler::asStateVecs(int32_t time, float centmass) {
  double n = meanAngularMotion();
  float period = 2 * pi / n;
  int32_t normalizedTime = ecc <= 1 ? fmod(time, period) : time;
  float men_ano = normalizedTime * n + mna;
  float ecc_ano = KepEqn::M_to_REAL(men_ano, ecc);
  float tru_ano = anomalyEcce2True(ecc_ano, ecc);

  float r;
  if (ecc <= 1)
    r = sma * (1 - ecc * cos(ecc_ano));
  else
    r = sma * (1 - ecc * ecc) / (1 + ecc * cos(tru_ano));

  Vector2D velocity;
  if (ecc <= 1)
    velocity = sqrt(CONST_G * centmass / r) * Vector2D{- sin(ecc_ano), sqrt(1 - ecc * ecc) * sin(ecc_ano)}.rotate(lpe);
  else {
    float phi = atan2(ecc * sin(tru_ano), 1 + ecc * cos(tru_ano));
    float vel = sqrt(CONST_G * centmass * (2 / r - 1 / sma));
    velocity = vel * Vector2D{phi + lpe};
  }

  return StateVec{r * Vector2D{tru_ano + lpe}, velocity};
}

int32_t Kepler::distanceNextWhen(int32_t time, float distance, float centmass) {
  // the inverse of asStateVecs; given r find rest
  float cosEcc, cosTru;
  float tmpAno;
  float candidateM0, candidateM1;
  if (ecc <= 1) {
    if (distance > sma) {
      return -1;
    }
    
    cosEcc = (sma - distance) / (sma * ecc);
    tmpAno = acos(cosEcc);
    candidateM0 = E_to_M(tmpAno, ecc);
  } else {
    cosTrue = - (ecc * ecc * sma + distance - sma) / (ecc * distance);
    tmpAno = acos(cosTrue);
    tmpAno = anomalyTrue2Ecce(tmpAno);
    candidateM0 = H_to_M(tmpAno, ecc);
  }

  // By symmetry
  candidateM1 = M_TAU - candidateM0;

  if (candidateM0 > candidateM1) {
    std::swap(candidateM0, candidateM1);
  }

  // Now we should somehow figure out which M comes first...
  double n = meanAngularMotion();
  float period = 2 * pi / n;
  int32_t normalizedTime = ecc <= 1 ? fmod(time, period) : time;
  float men_ano = normalizedTime * n + mna;
  if (ecc <= 1) {
    men_ano = fmod(men_ano, M_TAU);
  }

  if (men_ano <= candidateM0) {
    // M0
    return time + (candidateM0 - men_ano) / n;
  } else if (men_ano <= candidateM1) {
    // M1
    return time + (candidateM1 - men_ano) / n;
  } else if (ecc <= 1) {
    // M0
    return time + period + (candidateM1 - men_ano) / n;
  } else {
    // now i see we are keeping it signed so it's all gonna crash in 2038
    return -1;
  }
}

#pragma GCC optimize pop
