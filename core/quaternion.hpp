#pragma once
#include "geometry.hpp"
#include "transform.hpp"

namespace eric {

class Quaternion {
    public:
        Vector3f v;
        Float w;
    public:
        Quaternion();
        Quaternion(const Transform& t);
        ~Quaternion();

        Transform ToTransform() const;

        Quaternion operator*(Float f) const;
        Quaternion& operator*=(Float f);
        Quaternion& operator/=(Float f);
        Quaternion operator/(Float f) const;
        Quaternion& operator+=(const Quaternion &q);
        Quaternion operator+(const Quaternion &q) const;
        Quaternion& operator-=(const Quaternion &q);
        Quaternion operator-(const Quaternion &q) const;
        Quaternion operator-() const;
};

inline Float Dot(const Quaternion& q1, const Quaternion& q2) {
    return Dot(q1.v, q2.v) + q1.w * q2.w;
};

inline Quaternion Normalize(const Quaternion& q) {
    return q / std::sqrt(Dot(q, q));
};

inline Quaternion Slerp(const Quaternion &q1, const Quaternion &q2, Float t) {
    Float cosTheta = Dot(q1, q2);
    if (cosTheta > 0.995f) {
        return q1 * t + q2 * (1 - t);
    } else {
        Float theta = std::acos(Clamp(cosTheta, -1, 1));
        Float thetap = theta * t;
        Quaternion qperp = Normalize(q2 - q1 * cosTheta);
        return q1 * std::cos(thetap) + qperp * std::sin(thetap);
    }
};

};