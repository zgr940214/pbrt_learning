#include "quaternion.hpp"
#include <cmath>
#include <algorithm>

namespace eric {
Quaternion::Quaternion():v(0.f, 0.f, 0.f), w(1) {};

Quaternion::Quaternion(const Transform& trans) {
    Transform t(trans);
    w = sqrt(1 + t.m.m[0][0] + t.m.m[1][1] + t.m.m[2][2]) / 2;
    v.x = sqrt(1 + t.m.m[0][0] - t.m.m[1][1] - t.m.m[2][2]) / 2;
    v.y = sqrt(1 - t.m.m[0][0] + t.m.m[1][1] - t.m.m[2][2]) / 2;
    v.z = sqrt(1 - t.m.m[0][0] - t.m.m[1][1] + t.m.m[2][2]) / 2;
    v.x = copysign(v.x, t.m.m[2][1] - t.m.m[1][2]);
    v.y = copysign(v.y, t.m.m[0][2] - t.m.m[2][0]);
    v.z = copysign(v.z, t.m.m[1][0] - t.m.m[0][1]);
    Normalize(*this);
};

Quaternion::~Quaternion() {};

Transform Quaternion::ToTransform() const {
    Transform t;
    Float xx = v.x * v.x;
    Float yy = v.y * v.y;
    Float zz = v.z * v.z;
    Float xy = v.x * v.y;
    Float xz = v.x * v.z;
    Float yz = v.y * v.z;
    Float wx = w * v.x;
    Float wy = w * v.y;
    Float wz = w * v.z;

    t.m.m[0][0] = 1 - 2 * (yy + zz);
    t.m.m[0][1] = 2 * (xy + wz);
    t.m.m[0][2] = 2 * (xz - wy);
    t.m.m[1][0] = 2 * (xy - wz);
    t.m.m[1][1] = 1 - 2 * (xx + zz);
    t.m.m[1][2] = 2 * (yz + wx);
    t.m.m[2][0] = 2 * (xz + wy);
    t.m.m[2][1] = 2 * (yz - wx);
    t.m.m[2][2] = 1 - 2 * (xx + yy);

    return t;
};

Quaternion& Quaternion::operator+=(const Quaternion &q) {
    v += q.v;
    w += q.w;
    return *this;
};

Quaternion& Quaternion::operator/=(Float f) {
    v /= f;
    w /= f;
    return *this;
};

Quaternion Quaternion::operator/(Float f) const {
    Quaternion q(*this);
    q /= f;
    return q;
};

Quaternion& Quaternion::operator*=(Float f) {
    v *= f;
    w *= f;
    return *this;
};

Quaternion Quaternion::operator*(Float f) const {
    Quaternion q(*this);
    q *= f;
    return q;
};

Quaternion Quaternion::operator+(const Quaternion &q) const {
    Quaternion qq(*this);
    qq += q;
    return qq;
};

Quaternion& Quaternion::operator-=(const Quaternion &q) {
    v -= q.v;
    w -= q.w;
    return *this;
};

Quaternion Quaternion::operator-(const Quaternion &q) const {
    Quaternion qq(*this);
    qq -= q;
    return qq;
};

Quaternion Quaternion::operator-() const{
    Quaternion q(*this);
    q.v = -q.v;
    q.w = -q.w;
    return q;
};

};