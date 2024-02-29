#pragma once
#include <stdint.h>
#include <string.h>
#include <cmath>

namespace eric {

#ifndef M_PI
#define M_PI 3.141592654
#endif

#ifdef PBRT_FLOAT_AS_DOUBLE
  typedef double Float;
#else
  typedef float Float;
#endif  // PBRT_FLOAT_AS_DOUBLE

static constexpr Float MaxFloat = std::numeric_limits<Float>::max();
static constexpr Float Infinity = std::numeric_limits<Float>::infinity();
static constexpr Float MachineEpsilon =
std::numeric_limits<Float>::epsilon() * 0.5;

inline Float gamma(int n) {
    return (n * MachineEpsilon) / (1 - n * MachineEpsilon);
}

inline uint32_t FloatToBits(float number) {
    uint32_t bits;
    memcpy(&bits, &number, sizeof(float));
    return bits;
}

inline float BitsToFloat(uint32_t bits) {
    float f;
    memcpy(&f, &bits, sizeof(uint32_t));
    return f;
}

inline uint64_t FloatToBits(double number) {
    uint64_t bits;
    memcpy(&bits, &number, sizeof(double));
    return bits;
}

inline double BitsToFloat(uint64_t bits) {
    double d;
    memcpy(&d, &bits, sizeof(uint64_t));
    return d;
}

inline float NextFloatUp(float v) {
    if (std::isinf(v) && v > 0.f) {
        return v;
    } 
    if (v == -0.f) {
        v = 0.f;
    }

    uint32_t ui = FloatToBits(v);
    if (v > 0) {
        ui++;
    } else {
        ui--;
    }
    return BitsToFloat(ui);
}

inline double NextFloatUp(double v) {
    if (std::isinf(v) && v > 0.f) {
        return v;
    }
    if (v == -0.f) {
        v = 0.f;
    }
    uint64_t ui = FloatToBits(v);
    if (v > 0.f) {
        ui++;
    } else {
        ui--;
    }
    return BitsToFloat(ui);
}

inline float NextFLoatDown(float v) {
    if (std::isinf(v) && v < 0.f) {
        return v;
    }
    if (v == 0.f) {
        v = -0.f;
    }
    uint32_t ui = FloatToBits(v);
    if (v > 0.f) {
        ui--;
    } else {
        ui++;
    }
    return BitsToFloat(ui);
}

inline double NextFLoatDown(double v) {
    if (std::isinf(v) && v < 0.f) {
        return v;
    }
    if (v == 0.f) {
        v = -0.f;
    }
    uint64_t ui = FloatToBits(v);
    if (v > 0.f) {
        ui--;
    } else {
        ui++;
    }
    return BitsToFloat(ui);
}

class EFloat {
    public:
        Float v, err;
        Float high, low;
    public: 
        EFloat():v(0.f), err(0.f), high(0.f), low(0.f){};
        EFloat(Float v, Float err = 0.f);
        
        bool operator==(const EFloat ef) const {
            return ef.v == v;
        };

        Float UpperBound() const {
            return high;
        }

        Float LowerBound() const {
            return low;
        }

        operator double() const {
            return v;
        };

        operator float() const {
            return v;
        };
        
        EFloat operator-() const {
            return EFloat(-v, err);
        };

        EFloat operator+(const EFloat ef) const;
        EFloat operator-(const EFloat ef) const;
        EFloat operator+(Float f) const;
        EFloat operator-(Float f) const;
        EFloat operator*(const EFloat ef) const;
        EFloat operator/(const EFloat ef) const;

        EFloat operator+=(const EFloat ef);
        EFloat operator*=(const EFloat ef);
        EFloat operator-=(const EFloat ef);
        EFloat operator/=(const EFloat ef);

        Float GetAbsoluteError();

};

inline void Swap(EFloat *e1, EFloat *e2) {
    Float v = (Float)e1->v;
    Float err = (Float)e1->err;

    e1->v = e2->v;
    e1->err = e2->err;
    e2->v = v;
    e2->err = err;
}

inline EFloat operator*(float f, const EFloat ef) {
    return EFloat(f) * ef;
}

inline EFloat operator/(float f, const EFloat ef) {
    return EFloat(f) / ef;
}

inline EFloat operator+(float f, const EFloat ef) {
    return EFloat(f) + ef;
}

inline EFloat operator-(float f, const EFloat ef) {
    return EFloat(f) - ef;
}

/*
** delta = b^2 - 4ac
** t1 = (-b + sqrt(delta)) / (2 * a)   t2 = (-b - sqrt(delta)) / (2 * a);
*/
inline bool Quadratic(EFloat a, EFloat b, EFloat c, EFloat *t1, EFloat *t2) {
    double discrim = double(b) * double(b) - 4 * double(a) * double(c);
    if (discrim < 0.f) {
        return false;
    }

    double rootDiscrim = sqrt(discrim);
    EFloat floatRootDiscrim(rootDiscrim, MachineEpsilon * rootDiscrim);
    EFloat q;
    if (float(b) < 0.f) {
        q = -0.5f * (b - q);        
    } else {
        q = -0.5f * (b + q);
    }
    *t1 = q / a;
    *t2 = c / q;

    if ((Float)*t1 > (Float)*t2) {
        Swap(t1, t2);
    }
    
    return true;
};

};