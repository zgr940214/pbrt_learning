#include "transform.hpp"
#include <cassert>

namespace eric {

template <typename U>
Matrix4x4<U> Inverse(const Matrix4x4<U> &m) {
    Matrix4x4<U> inv;
    U det = m.determinant();
    if (det == 0) {
        throw std::runtime_error("Matrix is singular and cannot be inverted.");
    }
    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++) {
            inv.m[j][i] = m.cofactor(i, j) / det;
        }
    }
    return inv;
};

template <typename U>
Matrix4x4<U> Transpose(const Matrix4x4<U> &) {
    return Matrix4x4<U>(m[0][0], m[1][0], m[2][0], m[3][0],
                m[0][1], m[1][1], m[2][1], m[3][1],
                m[0][2], m[1][2], m[2][2], m[3][2],
                m[0][3], m[1][3], m[2][3], m[3][3]);
};

Transform Inverse(const Transform &t) {
    return Transform(t.mInv, t.m);
}

Transform Transpose(const Transform &t) {
    return Transform(Transpose(t.m), Transpose(t.mInv));
}

bool IsIdentity(const Transform &t) {
    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++) {
            if (i == j) {
                if (t.m.m[i][j] != 1) return false;
            } else {
                if (t.m.m[i][j] != 0) return false;
            }
        }
    }
    return true;
};

bool IsEqual(const Transform &t1, const Transform &t2) {
    return t1 == t2;
};

void AnimatedTransform::Decompose(const Transform* t, Vector3f *T, Quaternion *R, Matrix4f *S) {
    assert(t != NULL);
    Matrix4f m(t->m);
    // T
    T->x = m.m[0][3];
    T->y = m.m[1][3];
    T->z = m.m[2][3];
    
    for (int i = 0; i < 3; i++) {
        m.m[i][3] = 0.f;
    }
    m.m[3][3] = 1.f;
    //R
    
    Matrix4f Rt = m, Rit, Rt_next;
    int count = 0;
    Float norm = 0.f;
    do {
        Rit = Inverse(Transpose(Rt));
        Rt_next = (Rt + Rit) * 0.5f;

        {
            Float n;
            for (int i = 0; i < 3; i++) {
                n = 0.f;
                n += abs(Rt.m[i][0] - Rt_next.m[i][0]);
                n += abs(Rt.m[i][1] - Rt_next.m[i][1]);
                n += abs(Rt.m[i][2] - Rt_next.m[i][2]);
                norm = std::max(norm , n);
            } 
        }
        Rt = Rt_next;
    }while(++count < 100 && norm > 0.001f);

    *R = Quaternion(Rt);

    // S   m = RS;  S = R^-1 * m;
    Matrix4f RInv = Inverse(Rt);
    *S = Matrix4f::Mul(RInv, m);
};

void AnimatedTransform::Interpolate(Float time, Transform* t) const{
    Transform T, R, S;

    if (!actuallyAnimated || time < startTime) {
        *t = *startTransform;
        return;
    } 
    if (time > endTime) {
        *t = *endTransform;
        return;
    }

    // T linear
    Float dt = (time - startTime) / (endTime - startTime);
    Vector3f Tv = this->T[0] * (1 - dt) + this->T[1] * (dt);
    T.m.m[0][3] = Tv.x;
    T.m.m[1][3] = Tv.y;
    T.m.m[2][3] = Tv.z;

    //S linear
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            S.m.m[i][j] = Lerp(dt, this->S[0].m[i][j], this->S[1].m[i][j]);
        }
    }

    //R Sphere Linear
    Quaternion r1; 
    if (Dot(this->R[0], this->R[1]) < 0) {
        R = Slerp(this->R[0], -this->R[1], dt).ToTransform();
    } else {
        R = Slerp(this->R[0], this->R[1], dt).ToTransform();
    }

    *t = T * R * S;
};

Bounds3f AnimatedTransform::BoundPointMotion(const Point3f& p) const {
    Bounds3f bounds((*startTransform)(p), (*endTransform)(p));
    Float cosTheta = Dot(R[0], R[1]);
    Float theta = std::acos(Clamp(cosTheta, -1, 1));
    Float zeros[4];
    int zeroCount = 0;

    for (int c = 0; c < 3; c++) {
        IntervalFindZeros(c1[c].Eval(p), c2[c].Eval(p), 
            c3[c].Eval(p), c4[c].Eval(p), c5[c].Eval(p), 
            theta, Interval(0.f, 1.0f), zeros, &zeroCount);
        if (zeroCount > 0) {
            for (int i = 0; i < zeroCount; i++) {
                Transform t;
                Interpolate(Lerp(zeros[i], startTime, endTime), &t);
                Point3f pp = t(p);
                bounds = Union(bounds, Bounds3f(pp));
            }
        }
    }
    return bounds;
};

Bounds3f AnimatedTransform::MotionBounds(const Bounds3f& b) const{
    if (!actuallyAnimated) {
        return (*startTransform)(b);
    } else if (!hasRotation) {
        Bounds3f sb = (*startTransform)(b);
        Bounds3f eb = (*endTransform)(b);
        return Union(sb, eb);
    } else { //motion bounds for rotation
        Bounds3f bb(b);
        bb = Union(bb, BoundPointMotion(b.pMin));
        bb = Union(bb, BoundPointMotion(b.pMax));
        bb = Union(bb, BoundPointMotion(Point3f(b.pMin.x , b.pMin.y, b.pMax.z)));
        bb = Union(bb, BoundPointMotion(Point3f(b.pMin.x , b.pMax.y, b.pMax.z)));
        bb = Union(bb, BoundPointMotion(Point3f(b.pMin.x , b.pMax.y, b.pMin.z)));
        bb = Union(bb, BoundPointMotion(Point3f(b.pMax.x , b.pMax.y, b.pMin.z)));
        bb = Union(bb, BoundPointMotion(Point3f(b.pMax.x , b.pMin.y, b.pMax.z)));
        bb = Union(bb, BoundPointMotion(Point3f(b.pMax.x , b.pMin.y, b.pMin.z)));
        return bb;
    }
};

Point3f AnimatedTransform::operator()(Float time, const Point3f &p) {
    Transform f;
    Interpolate(time, &f);
    Point3f pp(p);
    return f(pp);
};
    
Vector3f AnimatedTransform::operator()(Float time, const Vector3f &v) {
    Transform f;
    Interpolate(time, &f);
    Vector3f vv(v);
    return f(vv);
};

Ray AnimatedTransform::operator()(const Ray &r) {
    Transform f;
    Interpolate(r.time, &f);
    Ray rr(r);
    return f(rr);
};

RayDifferential AnimatedTransform::operator()(const RayDifferential &r) {
    Transform f;
    Interpolate(r.time, &f);
    RayDifferential rr(r);
    return f(rr);
};

};