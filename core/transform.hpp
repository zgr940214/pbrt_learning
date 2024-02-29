#pragma once
#include "bounds.hpp"
#include "geometry.hpp"
#include "interaction.hpp"
#include "intervals.hpp"
#include "quaternion.hpp"
#include "ray.hpp"

#if defined(_WITH_SSE_)
#include <emmintrin.h>  // SSE2
#include <smmintrin.h>  // SSE4.1
#endif

namespace eric {

/*
-----------------------------------------------------------------------------------------
** Transformation:
-----------------------------------------------------------------------------------------
*/
template <typename T>
class Matrix4x4 {
 public:
  using data_t = T;
  using ptr_t = T *;
  using row_t = data_t[4];
  union {
    data_t m[4][4];
    row_t row[4];
  };

 public:
  Matrix4x4()
      : m[0][0](1),
        m[0][1](0),
        m[0][2](0),
        m[0][3](0),
        m[1][0](0),
        m[1][1](1),
        m[1][2](0),
        m[1][3](0),
        m[2][0](0),
        m[2][1](0),
        m[2][2](1),
        m[2][3](0),
        m[3][0](0),
        m[3][1](0),
        m[3][2](0),
        m[3][3](1){};
  Matrix4x4(data_t mat[4][4]) {
    std::copy(&mat[0][0], &mat[0][0] + 16, &m[0][0]);
  };
  Matrix4x4(data_t t00, data_t t01, data_t t02, data_t t03, data_t t10,
            data_t t11, data_t t12, data_t t13, data_t t20, data_t t21,
            data_t t22, data_t t23, data_t t30, data_t t31, data_t t32,
            data_t t33)
      : m[0][0](t00),
        m[0][1](t01),
        m[0][2](t02),
        m[0][3](t03),
        m[1][0](t10),
        m[1][1](t11),
        m[1][2](t12),
        m[1][3](t13),
        m[2][0](t20),
        m[2][1](t21),
        m[2][2](t22),
        m[2][3](t23),
        m[3][0](t30),
        m[3][1](t31),
        m[3][2](t32),
        m[3][3](t33){};

  Matrix4x4 transpose() {
    Matrix4x4 result;
#if defined(_WITH_SSE)
    __m128 row1 = _mm_loadu_ps(&m[0][0]);
    __m128 row2 = _mm_loadu_ps(&m[1][0]);
    __m128 row3 = _mm_loadu_ps(&m[2][0]);
    __m128 row4 = _mm_loadu_ps(&m[3][0]);

    _MM_TRANSPOSE4_PS(row1, row2, row3, row4);

    _mm_storeu_ps(&result.m[0][0], row1);
    _mm_storeu_ps(&result.m[1][0], row2);
    _mm_storeu_ps(&result.m[2][0], row3);
    _mm_storeu_ps(&result.m[3][0], row4);

#else
    result = Matrix4x4(m[0][0], m[1][0], m[2][0], m[3][0], m[0][1], m[1][1],
                       m[2][1], m[3][1], m[0][2], m[1][2], m[2][2], m[3][2],
                       m[0][3], m[1][3], m[2][3], m[3][3]);
#endif
    *this = result;
    return *this;
  };

  static Matrix4x4 Mul(const Matrix4x4 &m1, const Matrix4x4 &m2) {
#if defined(_WITH_SSE)
    Matrix4x4 result;
    for (int i = 0; i < 4; ++i) {
      __m128 row = _mm_loadu_ps(&m1.m[i][0]);
      for (int j = 0; j < 4; ++j) {
        __m128 col = _mm_set_ps(m2.m[3][j], m2.m[2][j], m2.m[1][j], m2.m[0][j]);
        __m128 res = _mm_mul_ps(row, col);
        result.m[i][j] = _mm_cvtss_f32(_mm_hadd_ps(_mm_hadd_ps(res, res), res));
      }
    }
    return result;
#else
    Matrix4x4 r;
    for (int i = 0; i < 4; ++i)
      for (int j = 0; j < 4; ++j)
        r.m[i][j] = m1.m[i][0] * m2.m[0][j] + m1.m[i][1] * m2.m[1][j] +
                    m1.m[i][2] * m2.m[2][j] + m1.m[i][3] * m2.m[3][j];
    return r;
#endif
  };

  Matrix4x4 mul(const Matrix4x4 &matrix) const {
#if defined(_WITH_SSE)
    Matrix4x4 result;
    for (int i = 0; i < 4; ++i) {
      __m128 row = _mm_loadu_ps(&m[i][0]);
      for (int j = 0; j < 4; ++j) {
        __m128 col = _mm_set_ps(matrix.m[3][j], matrix.m[2][j], matrix.m[1][j],
                                matrix.m[0][j]);
        __m128 res = _mm_mul_ps(row, col);
        result.m[i][j] = _mm_cvtss_f32(_mm_hadd_ps(_mm_hadd_ps(res, res), res));
      }
    }
    return result;
#else
    Matrix4x4 r;
    for (int i = 0; i < 4; ++i)
      for (int j = 0; j < 4; ++j)
        r.m[i][j] = m[i][0] * matrix.m[0][j] + m[i][1] * matrix.m[1][j] +
                    m[i][2] * matrix.m[2][j] + m[i][3] * matrix.m[3][j];
    return r;
#endif
  }

  T determinant() const {
#if defined(_WITH_SSE)
    __m128 row1 = _mm_loadu_ps(&mat.mat[0][0]);
    __m128 row2 = _mm_loadu_ps(&mat.mat[1][0]);
    __m128 row3 = _mm_loadu_ps(&mat.mat[2][0]);
    __m128 row4 = _mm_loadu_ps(&mat.mat[3][0]);

    // Transpose the matrix
    _MM_TRANSPOSE4_PS(row1, row2, row3, row4);

    // Calculate minors for the first row
    __m128 minor1 = _mm_sub_ps(
        _mm_mul_ps(_mm_shuffle_ps(row3, row3, _MM_SHUFFLE(3, 0, 2, 1)),
                   _mm_shuffle_ps(row4, row4, _MM_SHUFFLE(2, 3, 1, 0))),
        _mm_mul_ps(_mm_shuffle_ps(row3, row3, _MM_SHUFFLE(2, 3, 1, 0)),
                   _mm_shuffle_ps(row4, row4, _MM_SHUFFLE(3, 0, 2, 1))));
    __m128 minor2 = _mm_sub_ps(
        _mm_mul_ps(_mm_shuffle_ps(row2, row2, _MM_SHUFFLE(3, 0, 2, 1)),
                   _mm_shuffle_ps(row4, row4, _MM_SHUFFLE(2, 3, 1, 0))),
        _mm_mul_ps(_mm_shuffle_ps(row2, row2, _MM_SHUFFLE(2, 3, 1, 0)),
                   _mm_shuffle_ps(row4, row4, _MM_SHUFFLE(3, 0, 2, 1))));
    __m128 minor3 = _mm_sub_ps(
        _mm_mul_ps(_mm_shuffle_ps(row2, row2, _MM_SHUFFLE(3, 0, 2, 1)),
                   _mm_shuffle_ps(row3, row3, _MM_SHUFFLE(2, 3, 1, 0))),
        _mm_mul_ps(_mm_shuffle_ps(row2, row2, _MM_SHUFFLE(2, 3, 1, 0)),
                   _mm_shuffle_ps(row3, row3, _MM_SHUFFLE(3, 0, 2, 1))));

    // Calculate determinant
    __m128 det = _mm_mul_ps(row1, _mm_set_ps(-1.0f, 1.0f, -1.0f, 1.0f));
    det = _mm_add_ps(det, _mm_movehl_ps(det, det));
    det = _mm_add_ss(det, _mm_shuffle_ps(det, det, _MM_SHUFFLE(1, 1, 1, 1)));
    det = _mm_mul_ss(det,
                     _mm_shuffle_ps(minor1, minor2, _MM_SHUFFLE(0, 0, 0, 0)));

    return _mm_cvtss_f32(det);
#else
    return mat[0][0] *
               (mat[1][1] * (mat[2][2] * mat[3][3] - mat[2][3] * mat[3][2]) -
                mat[1][2] * (mat[2][1] * mat[3][3] - mat[2][3] * mat[3][1]) +
                mat[1][3] * (mat[2][1] * mat[3][2] - mat[2][2] * mat[3][1])) -
           mat[0][1] *
               (mat[1][0] * (mat[2][2] * mat[3][3] - mat[2][3] * mat[3][2]) -
                mat[1][2] * (mat[2][0] * mat[3][3] - mat[2][3] * mat[3][0]) +
                mat[1][3] * (mat[2][0] * mat[3][2] - mat[2][2] * mat[3][0])) +
           mat[0][2] *
               (mat[1][0] * (mat[2][1] * mat[3][3] - mat[2][3] * mat[3][1]) -
                mat[1][1] * (mat[2][0] * mat[3][3] - mat[2][3] * mat[3][0]) +
                mat[1][3] * (mat[2][0] * mat[3][1] - mat[2][1] * mat[3][0])) -
           mat[0][3] *
               (mat[1][0] * (mat[2][1] * mat[3][2] - mat[2][2] * mat[3][1]) -
                mat[1][1] * (mat[2][0] * mat[3][2] - mat[2][2] * mat[3][0]) +
                mat[1][2] * (mat[2][0] * mat[3][1] - mat[2][1] * mat[3][0]));
#endif
  };

  T cofactor(int i, int j) const {
    T submatrix[3][3];
    // Fill the submatrix by excluding row i and column j
    for (int row = 0, subRow = 0; row < 4; ++row) {
      if (row == i) continue;
      for (int col = 0, subCol = 0; col < 4; ++col) {
        if (col == j) continue;
        submatrix[subRow][subCol] = this->m[row][col];
        ++subCol;
      }
      ++subRow;
    }
    // Calculate the determinant of the submatrix directly
    T det = submatrix[0][0] * (submatrix[1][1] * submatrix[2][2] -
                               submatrix[1][2] * submatrix[2][1]) -
            submatrix[0][1] * (submatrix[1][0] * submatrix[2][2] -
                               submatrix[1][2] * submatrix[2][0]) +
            submatrix[0][2] * (submatrix[1][0] * submatrix[2][1] -
                               submatrix[1][1] * submatrix[2][0]);
    // Return the cofactor
    return ((i + j) % 2 == 0 ? 1 : -1) * det;
  };

  Matrix4x4 &operator+=(const Matrix4x4 &m) {
    for (int i = 0; i < 4; i++) {
      for (int j = 0; j < 4; j++) {
        m[i][j] += m.m[i][j];
      }
    }
    return *this;
  };

  Matrix4x4 operator+(const Matrix4x4 &m) const {
    Matrix4x4 result(*this);
    result += m;
    return result;
  };

  Matrix4x4 &operator-=(const Matrix4x4 &m) {
    for (int i = 0; i < 4; i++) {
      for (int j = 0; j < 4; j++) {
        m[i][j] -= m.m[i][j];
      }
    }
    return *this;
  };

  Matrix4x4 operator-(const Matrix4x4 &m) const {
    Matrix4x4 result(*this);
    result += m;
    return result;
  };

  Matrix4x4 &operator*=(T v) {
    for (int i = 0; i < 4; i++) {
      for (int j = 0; j < 4; j++) {
        m[i][j] *= v;
      }
    }
    return *this;
  };

  Matrix4x4 operator*(T v) const {
    Matrix4x4 m(*this);
    m *= v;
    return m;
  };

  Matrix4x4 &operator/=(T v) {
    for (int i = 0; i < 4; i++) {
      for (int j = 0; j < 4; j++) {
        m[i][j] /= v;
      }
    }
    return *this;
  };

  Matrix4x4 operator/(T v) const {
    Matrix4x4 m(*this);
    m /= v;
    return m;
  };

  template <typename U>
  friend Matrix4x4<U> Inverse(const Matrix4x4<U> &);

  template <typename U>
  friend Matrix4x4<U> Transpose(const Matrix4x4<U> &);
};

template <typename T>
class Matrix3x3 {
 public:
  // TODO: Mat3 is not used currently, implement it later
};

using Matrix3f = Matrix3x3<Float>;
using Matrix3i = Matrix3x3<int>;
using Matrix4f = Matrix4x4<Float>;
using Matrix4i = Matrix4x4<int>;

class Transform {
 public:
  Matrix4f m;
  Matrix4f mInv;

 public:
  Transform(){};
  Transform(const Float mat[4][4]) {
    m = Matrix4f(mat[0][0], mat[0][1], mat[0][2], mat[0][3], mat[1][0],
                 mat[1][1], mat[1][2], mat[1][3], mat[2][0], mat[2][1],
                 mat[2][2], mat[2][3], mat[3][0], mat[3][1], mat[3][2],
                 mat[3][3]);

    mInv = Inverse(m);
  };

  Transform(const Matrix4f &m) : m(m), mInv(Inverse(m)){};

  Transform(const Matrix4f &m, const Matrix4f &mInv) : m(m), mInv(mInv){};

  Vector3f operator()(const Vector3f &v) const {
    Float x = v.x, y = v.y, z = v.z;
    return Vector3f(m.m[0][0] * x + m.m[0][1] * y + m.m[0][2] * z,
                    m.m[1][0] * x + m.m[1][1] * y + m.m[1][2] * z,
                    m.m[2][0] * x + m.m[2][1] * y + m.m[2][2] * z);
  };

  Point3f operator()(const Point3f &p) const {
    Float x = p.x, y = p.y, z = p.z;
    Float w = m.m[3][0] * x + m.m[3][1] * y + m.m[3][2] * z + m.m[3][3];
    Point3f r =
        Point3f(m.m[0][0] * x + m.m[0][1] * y + m.m[0][2] * z + m.m[0][3],
                m.m[1][0] * x + m.m[1][1] * y + m.m[1][2] * z + m.m[1][3],
                m.m[2][0] * x + m.m[2][1] * y + m.m[2][2] * z + m.m[2][3]);
    if (w != 1 && w != 0) {
      r /= w;
    }
    return r;
  };

  Normal3f operator()(const Normal3f &n) const {
    Float x = n.x, y = n.y, z = n.z;
    return Normal3f(mInv.m[0][0] * x + mInv.m[1][0] * y + mInv.m[2][0] * z,
                    mInv.m[0][1] * x + mInv.m[1][1] * y + mInv.m[2][1] * z,
                    mInv.m[0][2] * x + mInv.m[1][2] * y + mInv.m[2][2] * z);
  };

  Ray operator()(const Ray &r) const {
    Vector3f oError;
    // TODO: oError
    Point3f o = (*this)(r.o);
    Vector3f d = (*this)(r.d);
    return Ray(o, d, r.tMax, r.time, r.medium);
  }

  // TODO: deal with float error
  Ray operator()(const Ray &r, Vector3f *oErr, Vector3f *dErr) const {

  };

  Bounds3f operator()(const Bounds3f &b) const {
    Point3f pmin = (*this)(b.pMin), pmax = (*this)(b.pMax);
    return Bounds3f(pmax, pmin);
  };

  Transform operator*(const Transform &t) const {
    return Transform(Matrix4f::Mul(m, t.m), Matrix4f::Mul(t.mInv, mInv));
  }

  // TODO: complete this one when introduce pError to handle floating number
  // error
  SurfaceInteraction operator()(const SurfaceInteraction &si) const {
    SurfaceInteraction ret;
    // Transform p and pError in SurfaceInteraction 229
    // Transform remaining members of SurfaceInteraction
    return ret;
  }

  bool operator==(const Transform &t) const {
    for (int i = 0; i < 4; i++) {
      for (int j = 0; j < 4; j++) {
        if (m.m[i][j] != t.m.m[i][j]) {
          return false;
        }
      }
    }
    return true;
  };

  bool has_scale() const {
    Float la2 = (*this)(Vector3f(1, 0, 0)).length_squared();
    Float lb2 = (*this)(Vector3f(0, 1, 0)).length_squared();
    Float lc2 = (*this)(Vector3f(0, 0, 1)).length_squared();
#define NOT_ONE(x) ((x) < .999f || (x) > 1.001f)
    return (NOT_ONE(la2) || NOT_ONE(lb2) || NOT_ONE(lc2));
#undef NOT_ONE
  }

  bool SwapsHandedness() const {
    Float det = m.m[0][0] * (m.m[1][1] * m.m[2][2] - m.m[1][2] * m.m[2][1]) -
                m.m[0][1] * (m.m[1][0] * m.m[2][2] - m.m[1][2] * m.m[2][0]) +
                m.m[0][2] * (m.m[1][0] * m.m[2][1] - m.m[1][1] * m.m[2][0]);
    return det < 0;
  }

  friend Transform Inverse(const Transform &t);

  friend Transform Transpose(const Transform &t);

  friend bool IsIdentity(const Transform &t);

  friend bool IsEqual(const Transform &t1, const Transform &t2);
};

/*
** Transform public method
*/
Transform Translate(const Vector3f &delta) {
  Matrix4f m(1, 0, 0, delta.x, 0, 1, 0, delta.y, 0, 0, 1, delta.z, 0, 0, 0, 1);
  Matrix4f minv(1, 0, 0, -delta.x, 0, 1, 0, -delta.y, 0, 0, 1, -delta.z, 0, 0,
                0, 1);
  return Transform(m, minv);
};

Transform Scale(Float x, Float y, Float z) {
  Matrix4f m(x, 0, 0, 0, 0, y, 0, 0, 0, 0, z, 0, 0, 0, 0, 1);
  Matrix4f minv(1 / x, 0, 0, 0, 0, 1 / y, 0, 0, 0, 0, 1 / z, 0, 0, 0, 0, 1);
  return Transform(m, minv);
};

Transform RotateX(Float theta) {
  Float sinTheta = sin(Radians(theta));
  Float cosTheta = cos(Radians(theta));
  return Transform(Matrix4f(1.0f, 0.f, 0.f, 0.f, 0.f, cosTheta, -sinTheta, 0.f,
                            0.f, sinTheta, cosTheta, 0.f, 0.f, 0.f, 0.f, 1.f));
};

Transform RotateY(Float theta) {
  Float sinTheta = sin(Radians(theta));
  Float cosTheta = cos(Radians(theta));
  return Transform(Matrix4f(cosTheta, 0.f, sinTheta, 0.f, 0.f, 1.f, 0.f, 0.f,
                            -sinTheta, 0.f, cosTheta, 0.f, 0.f, 0.f, 0.f, 1.f));
};

Transform RotateZ(Float theta) {
  Float sinTheta = sin(Radians(theta));
  Float cosTheta = cos(Radians(theta));
  return Transform(Matrix4f(cosTheta, -sinTheta, 0.f, 0.f, sinTheta, cosTheta,
                            0.f, 0.f, 0.f, 0.f, 1.f, 0.f, 0.f, 0.f, 0.f, 1.f));
};

Transform Rotate(Float theta, Vector3f axis) {
  Float cosTheta = cos(Radians(theta));
  Float sinTheta = sin(Radians(theta));
  Matrix4f m =
      Matrix4f(cosTheta + axis.x * axis.x * (1 - cosTheta),
               axis.x * axis.y * (1 - cosTheta) - axis.z * sinTheta,
               axis.x * axis.z * (1 - cosTheta) + axis.y * sinTheta, 0.f,

               axis.y * axis.x * (1 - cosTheta) + axis.z * sinTheta,
               cosTheta + axis.y * axis.y * (1 - cosTheta),
               axis.y * axis.z * (1 - cosTheta) - axis.x * sinTheta, 0.f,

               axis.z * axis.x * (1 - cosTheta) - axis.y * sinTheta,
               axis.z * axis.y * (1 - cosTheta) + axis.x * sinTheta,
               cosTheta + axis.z * axis.z * (1 - cosTheta), 0.f,

               0.f, 0.f, 0.f, 1.f);
  return Transform(m);
}

Transform LookAt(const Point3f &pos, const Point3f &look, const Vector3f &up) {
  Vector3f dir = Normalize(look - pos);
  Vector3f left = Normalize(Cross(Normalize(up), dir));
  Vector3f newUp = Cross(dir, left);
  Matrix4f camera2Wld =
      Matrix4f(left.x, newUp.x, dir.x, pos.x, left.y, newUp.y, dir.y, pos.y,
               left.z, newUp.z, dir.z, pos.z, 0.f, 0.f, 0.f, 1.f);

  return Transform(camera2Wld);
};

class AnimationTransform {
 private:
  struct DerivativeTerm {
    Float kc, kx, ky, kz;

    DerivativeTerm(){};
    DerivativeTerm(Float c, Float x, Float y, Float z)
        : kc(c), kx(x), ky(y), kz(z){};

    Float Eval(const Point3f &p) const {
      return kc + kx * p.x + ky * p.y + kz * p.z;
    };
  };
  const Transform *startTransform, *endTransform;
  const Float startTime, endTime;
  const bool actuallyAnimated;
  Vector3f T[2];    // translate
  Quaternion R[2];  // rotation
  Matrix4f S[2];    // scale
  bool hasRotation;
  DerivativeTerm c1[3], c2[3], c3[3], c4[3], c5[3];

 public:
  AnimationTransform(const Transform *s, const Transform *e, Float start,
                     Float end)
      : startTransform(s),
        endTransform(e),
        startTime(start),
        endTime(end),
        actuallyAnimated(s != e) {
    Decompose(s, &T[0], &R[0], &S[0]);
    Decompose(e, &T[1], &R[1], &S[1]);
    hasRotation = (Dot(R[0], R[1]) < 0.995f);
    if (hasRotation) {
      Float cosTheta = Dot(R[0], R[1]);
      Float theta = std::acos(Clamp(cosTheta, -1, 1));
      Quaternion qperp = Normalize(R[1] - R[0] * cosTheta);

      Float t0x = T[0].x;
      Float t0y = T[0].y;
      Float t0z = T[0].z;
      Float t1x = T[1].x;
      Float t1y = T[1].y;
      Float t1z = T[1].z;
      Float q0x = R[0].v.x;
      Float q0y = R[0].v.y;
      Float q0z = R[0].v.z;
      Float q0w = R[0].w;
      Float qperpx = qperp.v.x;
      Float qperpy = qperp.v.y;
      Float qperpz = qperp.v.z;
      Float qperpw = qperp.w;
      Float s000 = S[0].m[0][0];
      Float s001 = S[0].m[0][1];
      Float s002 = S[0].m[0][2];
      Float s010 = S[0].m[1][0];
      Float s011 = S[0].m[1][1];
      Float s012 = S[0].m[1][2];
      Float s020 = S[0].m[2][0];
      Float s021 = S[0].m[2][1];
      Float s022 = S[0].m[2][2];
      Float s100 = S[1].m[0][0];
      Float s101 = S[1].m[0][1];
      Float s102 = S[1].m[0][2];
      Float s110 = S[1].m[1][0];
      Float s111 = S[1].m[1][1];
      Float s112 = S[1].m[1][2];
      Float s120 = S[1].m[2][0];
      Float s121 = S[1].m[2][1];
      Float s122 = S[1].m[2][2];

      c1[0] = DerivativeTerm(
          -t0x + t1x,
          (-1 + q0y * q0y + q0z * q0z + qperpy * qperpy + qperpz * qperpz) *
                  s000 +
              q0w * q0z * s010 - qperpx * qperpy * s010 +
              qperpw * qperpz * s010 - q0w * q0y * s020 -
              qperpw * qperpy * s020 - qperpx * qperpz * s020 + s100 -
              q0y * q0y * s100 - q0z * q0z * s100 - qperpy * qperpy * s100 -
              qperpz * qperpz * s100 - q0w * q0z * s110 +
              qperpx * qperpy * s110 - qperpw * qperpz * s110 +
              q0w * q0y * s120 + qperpw * qperpy * s120 +
              qperpx * qperpz * s120 +
              q0x * (-(q0y * s010) - q0z * s020 + q0y * s110 + q0z * s120),
          (-1 + q0y * q0y + q0z * q0z + qperpy * qperpy + qperpz * qperpz) *
                  s001 +
              q0w * q0z * s011 - qperpx * qperpy * s011 +
              qperpw * qperpz * s011 - q0w * q0y * s021 -
              qperpw * qperpy * s021 - qperpx * qperpz * s021 + s101 -
              q0y * q0y * s101 - q0z * q0z * s101 - qperpy * qperpy * s101 -
              qperpz * qperpz * s101 - q0w * q0z * s111 +
              qperpx * qperpy * s111 - qperpw * qperpz * s111 +
              q0w * q0y * s121 + qperpw * qperpy * s121 +
              qperpx * qperpz * s121 +
              q0x * (-(q0y * s011) - q0z * s021 + q0y * s111 + q0z * s121),
          (-1 + q0y * q0y + q0z * q0z + qperpy * qperpy + qperpz * qperpz) *
                  s002 +
              q0w * q0z * s012 - qperpx * qperpy * s012 +
              qperpw * qperpz * s012 - q0w * q0y * s022 -
              qperpw * qperpy * s022 - qperpx * qperpz * s022 + s102 -
              q0y * q0y * s102 - q0z * q0z * s102 - qperpy * qperpy * s102 -
              qperpz * qperpz * s102 - q0w * q0z * s112 +
              qperpx * qperpy * s112 - qperpw * qperpz * s112 +
              q0w * q0y * s122 + qperpw * qperpy * s122 +
              qperpx * qperpz * s122 +
              q0x * (-(q0y * s012) - q0z * s022 + q0y * s112 + q0z * s122));

      c2[0] = DerivativeTerm(
          0.,
          -(qperpy * qperpy * s000) - qperpz * qperpz * s000 +
              qperpx * qperpy * s010 - qperpw * qperpz * s010 +
              qperpw * qperpy * s020 + qperpx * qperpz * s020 +
              q0y * q0y * (s000 - s100) + q0z * q0z * (s000 - s100) +
              qperpy * qperpy * s100 + qperpz * qperpz * s100 -
              qperpx * qperpy * s110 + qperpw * qperpz * s110 -
              qperpw * qperpy * s120 - qperpx * qperpz * s120 +
              2 * q0x * qperpy * s010 * theta -
              2 * q0w * qperpz * s010 * theta +
              2 * q0w * qperpy * s020 * theta +
              2 * q0x * qperpz * s020 * theta +
              q0y * (q0x * (-s010 + s110) + q0w * (-s020 + s120) +
                     2 * (-2 * qperpy * s000 + qperpx * s010 + qperpw * s020) *
                         theta) +
              q0z * (q0w * (s010 - s110) + q0x * (-s020 + s120) -
                     2 * (2 * qperpz * s000 + qperpw * s010 - qperpx * s020) *
                         theta),
          -(qperpy * qperpy * s001) - qperpz * qperpz * s001 +
              qperpx * qperpy * s011 - qperpw * qperpz * s011 +
              qperpw * qperpy * s021 + qperpx * qperpz * s021 +
              q0y * q0y * (s001 - s101) + q0z * q0z * (s001 - s101) +
              qperpy * qperpy * s101 + qperpz * qperpz * s101 -
              qperpx * qperpy * s111 + qperpw * qperpz * s111 -
              qperpw * qperpy * s121 - qperpx * qperpz * s121 +
              2 * q0x * qperpy * s011 * theta -
              2 * q0w * qperpz * s011 * theta +
              2 * q0w * qperpy * s021 * theta +
              2 * q0x * qperpz * s021 * theta +
              q0y * (q0x * (-s011 + s111) + q0w * (-s021 + s121) +
                     2 * (-2 * qperpy * s001 + qperpx * s011 + qperpw * s021) *
                         theta) +
              q0z * (q0w * (s011 - s111) + q0x * (-s021 + s121) -
                     2 * (2 * qperpz * s001 + qperpw * s011 - qperpx * s021) *
                         theta),
          -(qperpy * qperpy * s002) - qperpz * qperpz * s002 +
              qperpx * qperpy * s012 - qperpw * qperpz * s012 +
              qperpw * qperpy * s022 + qperpx * qperpz * s022 +
              q0y * q0y * (s002 - s102) + q0z * q0z * (s002 - s102) +
              qperpy * qperpy * s102 + qperpz * qperpz * s102 -
              qperpx * qperpy * s112 + qperpw * qperpz * s112 -
              qperpw * qperpy * s122 - qperpx * qperpz * s122 +
              2 * q0x * qperpy * s012 * theta -
              2 * q0w * qperpz * s012 * theta +
              2 * q0w * qperpy * s022 * theta +
              2 * q0x * qperpz * s022 * theta +
              q0y * (q0x * (-s012 + s112) + q0w * (-s022 + s122) +
                     2 * (-2 * qperpy * s002 + qperpx * s012 + qperpw * s022) *
                         theta) +
              q0z * (q0w * (s012 - s112) + q0x * (-s022 + s122) -
                     2 * (2 * qperpz * s002 + qperpw * s012 - qperpx * s022) *
                         theta));

      c3[0] = DerivativeTerm(
          0.,
          -2 *
              (q0x * qperpy * s010 - q0w * qperpz * s010 + q0w * qperpy * s020 +
               q0x * qperpz * s020 - q0x * qperpy * s110 + q0w * qperpz * s110 -
               q0w * qperpy * s120 - q0x * qperpz * s120 +
               q0y * (-2 * qperpy * s000 + qperpx * s010 + qperpw * s020 +
                      2 * qperpy * s100 - qperpx * s110 - qperpw * s120) +
               q0z * (-2 * qperpz * s000 - qperpw * s010 + qperpx * s020 +
                      2 * qperpz * s100 + qperpw * s110 - qperpx * s120)) *
              theta,
          -2 *
              (q0x * qperpy * s011 - q0w * qperpz * s011 + q0w * qperpy * s021 +
               q0x * qperpz * s021 - q0x * qperpy * s111 + q0w * qperpz * s111 -
               q0w * qperpy * s121 - q0x * qperpz * s121 +
               q0y * (-2 * qperpy * s001 + qperpx * s011 + qperpw * s021 +
                      2 * qperpy * s101 - qperpx * s111 - qperpw * s121) +
               q0z * (-2 * qperpz * s001 - qperpw * s011 + qperpx * s021 +
                      2 * qperpz * s101 + qperpw * s111 - qperpx * s121)) *
              theta,
          -2 *
              (q0x * qperpy * s012 - q0w * qperpz * s012 + q0w * qperpy * s022 +
               q0x * qperpz * s022 - q0x * qperpy * s112 + q0w * qperpz * s112 -
               q0w * qperpy * s122 - q0x * qperpz * s122 +
               q0y * (-2 * qperpy * s002 + qperpx * s012 + qperpw * s022 +
                      2 * qperpy * s102 - qperpx * s112 - qperpw * s122) +
               q0z * (-2 * qperpz * s002 - qperpw * s012 + qperpx * s022 +
                      2 * qperpz * s102 + qperpw * s112 - qperpx * s122)) *
              theta);

      c4[0] = DerivativeTerm(
          0.,
          -(q0x * qperpy * s010) + q0w * qperpz * s010 - q0w * qperpy * s020 -
              q0x * qperpz * s020 + q0x * qperpy * s110 - q0w * qperpz * s110 +
              q0w * qperpy * s120 + q0x * qperpz * s120 +
              2 * q0y * q0y * s000 * theta + 2 * q0z * q0z * s000 * theta -
              2 * qperpy * qperpy * s000 * theta -
              2 * qperpz * qperpz * s000 * theta +
              2 * qperpx * qperpy * s010 * theta -
              2 * qperpw * qperpz * s010 * theta +
              2 * qperpw * qperpy * s020 * theta +
              2 * qperpx * qperpz * s020 * theta +
              q0y *
                  (-(qperpx * s010) - qperpw * s020 +
                   2 * qperpy * (s000 - s100) + qperpx * s110 + qperpw * s120 -
                   2 * q0x * s010 * theta - 2 * q0w * s020 * theta) +
              q0z * (2 * qperpz * s000 + qperpw * s010 - qperpx * s020 -
                     2 * qperpz * s100 - qperpw * s110 + qperpx * s120 +
                     2 * q0w * s010 * theta - 2 * q0x * s020 * theta),
          -(q0x * qperpy * s011) + q0w * qperpz * s011 - q0w * qperpy * s021 -
              q0x * qperpz * s021 + q0x * qperpy * s111 - q0w * qperpz * s111 +
              q0w * qperpy * s121 + q0x * qperpz * s121 +
              2 * q0y * q0y * s001 * theta + 2 * q0z * q0z * s001 * theta -
              2 * qperpy * qperpy * s001 * theta -
              2 * qperpz * qperpz * s001 * theta +
              2 * qperpx * qperpy * s011 * theta -
              2 * qperpw * qperpz * s011 * theta +
              2 * qperpw * qperpy * s021 * theta +
              2 * qperpx * qperpz * s021 * theta +
              q0y *
                  (-(qperpx * s011) - qperpw * s021 +
                   2 * qperpy * (s001 - s101) + qperpx * s111 + qperpw * s121 -
                   2 * q0x * s011 * theta - 2 * q0w * s021 * theta) +
              q0z * (2 * qperpz * s001 + qperpw * s011 - qperpx * s021 -
                     2 * qperpz * s101 - qperpw * s111 + qperpx * s121 +
                     2 * q0w * s011 * theta - 2 * q0x * s021 * theta),
          -(q0x * qperpy * s012) + q0w * qperpz * s012 - q0w * qperpy * s022 -
              q0x * qperpz * s022 + q0x * qperpy * s112 - q0w * qperpz * s112 +
              q0w * qperpy * s122 + q0x * qperpz * s122 +
              2 * q0y * q0y * s002 * theta + 2 * q0z * q0z * s002 * theta -
              2 * qperpy * qperpy * s002 * theta -
              2 * qperpz * qperpz * s002 * theta +
              2 * qperpx * qperpy * s012 * theta -
              2 * qperpw * qperpz * s012 * theta +
              2 * qperpw * qperpy * s022 * theta +
              2 * qperpx * qperpz * s022 * theta +
              q0y *
                  (-(qperpx * s012) - qperpw * s022 +
                   2 * qperpy * (s002 - s102) + qperpx * s112 + qperpw * s122 -
                   2 * q0x * s012 * theta - 2 * q0w * s022 * theta) +
              q0z * (2 * qperpz * s002 + qperpw * s012 - qperpx * s022 -
                     2 * qperpz * s102 - qperpw * s112 + qperpx * s122 +
                     2 * q0w * s012 * theta - 2 * q0x * s022 * theta));

      c5[0] = DerivativeTerm(
          0.,
          2 *
              (qperpy * qperpy * s000 + qperpz * qperpz * s000 -
               qperpx * qperpy * s010 + qperpw * qperpz * s010 -
               qperpw * qperpy * s020 - qperpx * qperpz * s020 -
               qperpy * qperpy * s100 - qperpz * qperpz * s100 +
               q0y * q0y * (-s000 + s100) + q0z * q0z * (-s000 + s100) +
               qperpx * qperpy * s110 - qperpw * qperpz * s110 +
               q0y * (q0x * (s010 - s110) + q0w * (s020 - s120)) +
               qperpw * qperpy * s120 + qperpx * qperpz * s120 +
               q0z * (-(q0w * s010) + q0x * s020 + q0w * s110 - q0x * s120)) *
              theta,
          2 *
              (qperpy * qperpy * s001 + qperpz * qperpz * s001 -
               qperpx * qperpy * s011 + qperpw * qperpz * s011 -
               qperpw * qperpy * s021 - qperpx * qperpz * s021 -
               qperpy * qperpy * s101 - qperpz * qperpz * s101 +
               q0y * q0y * (-s001 + s101) + q0z * q0z * (-s001 + s101) +
               qperpx * qperpy * s111 - qperpw * qperpz * s111 +
               q0y * (q0x * (s011 - s111) + q0w * (s021 - s121)) +
               qperpw * qperpy * s121 + qperpx * qperpz * s121 +
               q0z * (-(q0w * s011) + q0x * s021 + q0w * s111 - q0x * s121)) *
              theta,
          2 *
              (qperpy * qperpy * s002 + qperpz * qperpz * s002 -
               qperpx * qperpy * s012 + qperpw * qperpz * s012 -
               qperpw * qperpy * s022 - qperpx * qperpz * s022 -
               qperpy * qperpy * s102 - qperpz * qperpz * s102 +
               q0y * q0y * (-s002 + s102) + q0z * q0z * (-s002 + s102) +
               qperpx * qperpy * s112 - qperpw * qperpz * s112 +
               q0y * (q0x * (s012 - s112) + q0w * (s022 - s122)) +
               qperpw * qperpy * s122 + qperpx * qperpz * s122 +
               q0z * (-(q0w * s012) + q0x * s022 + q0w * s112 - q0x * s122)) *
              theta);

      c1[1] = DerivativeTerm(
          -t0y + t1y,
          -(qperpx * qperpy * s000) - qperpw * qperpz * s000 - s010 +
              q0z * q0z * s010 + qperpx * qperpx * s010 +
              qperpz * qperpz * s010 - q0y * q0z * s020 +
              qperpw * qperpx * s020 - qperpy * qperpz * s020 +
              qperpx * qperpy * s100 + qperpw * qperpz * s100 +
              q0w * q0z * (-s000 + s100) + q0x * q0x * (s010 - s110) + s110 -
              q0z * q0z * s110 - qperpx * qperpx * s110 -
              qperpz * qperpz * s110 +
              q0x * (q0y * (-s000 + s100) + q0w * (s020 - s120)) +
              q0y * q0z * s120 - qperpw * qperpx * s120 +
              qperpy * qperpz * s120,
          -(qperpx * qperpy * s001) - qperpw * qperpz * s001 - s011 +
              q0z * q0z * s011 + qperpx * qperpx * s011 +
              qperpz * qperpz * s011 - q0y * q0z * s021 +
              qperpw * qperpx * s021 - qperpy * qperpz * s021 +
              qperpx * qperpy * s101 + qperpw * qperpz * s101 +
              q0w * q0z * (-s001 + s101) + q0x * q0x * (s011 - s111) + s111 -
              q0z * q0z * s111 - qperpx * qperpx * s111 -
              qperpz * qperpz * s111 +
              q0x * (q0y * (-s001 + s101) + q0w * (s021 - s121)) +
              q0y * q0z * s121 - qperpw * qperpx * s121 +
              qperpy * qperpz * s121,
          -(qperpx * qperpy * s002) - qperpw * qperpz * s002 - s012 +
              q0z * q0z * s012 + qperpx * qperpx * s012 +
              qperpz * qperpz * s012 - q0y * q0z * s022 +
              qperpw * qperpx * s022 - qperpy * qperpz * s022 +
              qperpx * qperpy * s102 + qperpw * qperpz * s102 +
              q0w * q0z * (-s002 + s102) + q0x * q0x * (s012 - s112) + s112 -
              q0z * q0z * s112 - qperpx * qperpx * s112 -
              qperpz * qperpz * s112 +
              q0x * (q0y * (-s002 + s102) + q0w * (s022 - s122)) +
              q0y * q0z * s122 - qperpw * qperpx * s122 +
              qperpy * qperpz * s122);

      c2[1] = DerivativeTerm(
          0.,
          qperpx * qperpy * s000 + qperpw * qperpz * s000 + q0z * q0z * s010 -
              qperpx * qperpx * s010 - qperpz * qperpz * s010 -
              q0y * q0z * s020 - qperpw * qperpx * s020 +
              qperpy * qperpz * s020 - qperpx * qperpy * s100 -
              qperpw * qperpz * s100 + q0x * q0x * (s010 - s110) -
              q0z * q0z * s110 + qperpx * qperpx * s110 +
              qperpz * qperpz * s110 + q0y * q0z * s120 +
              qperpw * qperpx * s120 - qperpy * qperpz * s120 +
              2 * q0z * qperpw * s000 * theta +
              2 * q0y * qperpx * s000 * theta -
              4 * q0z * qperpz * s010 * theta +
              2 * q0z * qperpy * s020 * theta +
              2 * q0y * qperpz * s020 * theta +
              q0x * (q0w * s020 + q0y * (-s000 + s100) - q0w * s120 +
                     2 * qperpy * s000 * theta - 4 * qperpx * s010 * theta -
                     2 * qperpw * s020 * theta) +
              q0w * (-(q0z * s000) + q0z * s100 + 2 * qperpz * s000 * theta -
                     2 * qperpx * s020 * theta),
          qperpx * qperpy * s001 + qperpw * qperpz * s001 + q0z * q0z * s011 -
              qperpx * qperpx * s011 - qperpz * qperpz * s011 -
              q0y * q0z * s021 - qperpw * qperpx * s021 +
              qperpy * qperpz * s021 - qperpx * qperpy * s101 -
              qperpw * qperpz * s101 + q0x * q0x * (s011 - s111) -
              q0z * q0z * s111 + qperpx * qperpx * s111 +
              qperpz * qperpz * s111 + q0y * q0z * s121 +
              qperpw * qperpx * s121 - qperpy * qperpz * s121 +
              2 * q0z * qperpw * s001 * theta +
              2 * q0y * qperpx * s001 * theta -
              4 * q0z * qperpz * s011 * theta +
              2 * q0z * qperpy * s021 * theta +
              2 * q0y * qperpz * s021 * theta +
              q0x * (q0w * s021 + q0y * (-s001 + s101) - q0w * s121 +
                     2 * qperpy * s001 * theta - 4 * qperpx * s011 * theta -
                     2 * qperpw * s021 * theta) +
              q0w * (-(q0z * s001) + q0z * s101 + 2 * qperpz * s001 * theta -
                     2 * qperpx * s021 * theta),
          qperpx * qperpy * s002 + qperpw * qperpz * s002 + q0z * q0z * s012 -
              qperpx * qperpx * s012 - qperpz * qperpz * s012 -
              q0y * q0z * s022 - qperpw * qperpx * s022 +
              qperpy * qperpz * s022 - qperpx * qperpy * s102 -
              qperpw * qperpz * s102 + q0x * q0x * (s012 - s112) -
              q0z * q0z * s112 + qperpx * qperpx * s112 +
              qperpz * qperpz * s112 + q0y * q0z * s122 +
              qperpw * qperpx * s122 - qperpy * qperpz * s122 +
              2 * q0z * qperpw * s002 * theta +
              2 * q0y * qperpx * s002 * theta -
              4 * q0z * qperpz * s012 * theta +
              2 * q0z * qperpy * s022 * theta +
              2 * q0y * qperpz * s022 * theta +
              q0x * (q0w * s022 + q0y * (-s002 + s102) - q0w * s122 +
                     2 * qperpy * s002 * theta - 4 * qperpx * s012 * theta -
                     2 * qperpw * s022 * theta) +
              q0w * (-(q0z * s002) + q0z * s102 + 2 * qperpz * s002 * theta -
                     2 * qperpx * s022 * theta));

      c3[1] =
          DerivativeTerm(0.,
                         2 *
                             (-(q0x * qperpy * s000) - q0w * qperpz * s000 +
                              2 * q0x * qperpx * s010 + q0x * qperpw * s020 +
                              q0w * qperpx * s020 + q0x * qperpy * s100 +
                              q0w * qperpz * s100 - 2 * q0x * qperpx * s110 -
                              q0x * qperpw * s120 - q0w * qperpx * s120 +
                              q0z * (2 * qperpz * s010 - qperpy * s020 +
                                     qperpw * (-s000 + s100) -
                                     2 * qperpz * s110 + qperpy * s120) +
                              q0y * (-(qperpx * s000) - qperpz * s020 +
                                     qperpx * s100 + qperpz * s120)) *
                             theta,
                         2 *
                             (-(q0x * qperpy * s001) - q0w * qperpz * s001 +
                              2 * q0x * qperpx * s011 + q0x * qperpw * s021 +
                              q0w * qperpx * s021 + q0x * qperpy * s101 +
                              q0w * qperpz * s101 - 2 * q0x * qperpx * s111 -
                              q0x * qperpw * s121 - q0w * qperpx * s121 +
                              q0z * (2 * qperpz * s011 - qperpy * s021 +
                                     qperpw * (-s001 + s101) -
                                     2 * qperpz * s111 + qperpy * s121) +
                              q0y * (-(qperpx * s001) - qperpz * s021 +
                                     qperpx * s101 + qperpz * s121)) *
                             theta,
                         2 *
                             (-(q0x * qperpy * s002) - q0w * qperpz * s002 +
                              2 * q0x * qperpx * s012 + q0x * qperpw * s022 +
                              q0w * qperpx * s022 + q0x * qperpy * s102 +
                              q0w * qperpz * s102 - 2 * q0x * qperpx * s112 -
                              q0x * qperpw * s122 - q0w * qperpx * s122 +
                              q0z * (2 * qperpz * s012 - qperpy * s022 +
                                     qperpw * (-s002 + s102) -
                                     2 * qperpz * s112 + qperpy * s122) +
                              q0y * (-(qperpx * s002) - qperpz * s022 +
                                     qperpx * s102 + qperpz * s122)) *
                             theta);

      c4[1] = DerivativeTerm(
          0.,
          -(q0x * qperpy * s000) - q0w * qperpz * s000 +
              2 * q0x * qperpx * s010 + q0x * qperpw * s020 +
              q0w * qperpx * s020 + q0x * qperpy * s100 + q0w * qperpz * s100 -
              2 * q0x * qperpx * s110 - q0x * qperpw * s120 -
              q0w * qperpx * s120 + 2 * qperpx * qperpy * s000 * theta +
              2 * qperpw * qperpz * s000 * theta +
              2 * q0x * q0x * s010 * theta + 2 * q0z * q0z * s010 * theta -
              2 * qperpx * qperpx * s010 * theta -
              2 * qperpz * qperpz * s010 * theta +
              2 * q0w * q0x * s020 * theta -
              2 * qperpw * qperpx * s020 * theta +
              2 * qperpy * qperpz * s020 * theta +
              q0y * (-(qperpx * s000) - qperpz * s020 + qperpx * s100 +
                     qperpz * s120 - 2 * q0x * s000 * theta) +
              q0z *
                  (2 * qperpz * s010 - qperpy * s020 + qperpw * (-s000 + s100) -
                   2 * qperpz * s110 + qperpy * s120 - 2 * q0w * s000 * theta -
                   2 * q0y * s020 * theta),
          -(q0x * qperpy * s001) - q0w * qperpz * s001 +
              2 * q0x * qperpx * s011 + q0x * qperpw * s021 +
              q0w * qperpx * s021 + q0x * qperpy * s101 + q0w * qperpz * s101 -
              2 * q0x * qperpx * s111 - q0x * qperpw * s121 -
              q0w * qperpx * s121 + 2 * qperpx * qperpy * s001 * theta +
              2 * qperpw * qperpz * s001 * theta +
              2 * q0x * q0x * s011 * theta + 2 * q0z * q0z * s011 * theta -
              2 * qperpx * qperpx * s011 * theta -
              2 * qperpz * qperpz * s011 * theta +
              2 * q0w * q0x * s021 * theta -
              2 * qperpw * qperpx * s021 * theta +
              2 * qperpy * qperpz * s021 * theta +
              q0y * (-(qperpx * s001) - qperpz * s021 + qperpx * s101 +
                     qperpz * s121 - 2 * q0x * s001 * theta) +
              q0z *
                  (2 * qperpz * s011 - qperpy * s021 + qperpw * (-s001 + s101) -
                   2 * qperpz * s111 + qperpy * s121 - 2 * q0w * s001 * theta -
                   2 * q0y * s021 * theta),
          -(q0x * qperpy * s002) - q0w * qperpz * s002 +
              2 * q0x * qperpx * s012 + q0x * qperpw * s022 +
              q0w * qperpx * s022 + q0x * qperpy * s102 + q0w * qperpz * s102 -
              2 * q0x * qperpx * s112 - q0x * qperpw * s122 -
              q0w * qperpx * s122 + 2 * qperpx * qperpy * s002 * theta +
              2 * qperpw * qperpz * s002 * theta +
              2 * q0x * q0x * s012 * theta + 2 * q0z * q0z * s012 * theta -
              2 * qperpx * qperpx * s012 * theta -
              2 * qperpz * qperpz * s012 * theta +
              2 * q0w * q0x * s022 * theta -
              2 * qperpw * qperpx * s022 * theta +
              2 * qperpy * qperpz * s022 * theta +
              q0y * (-(qperpx * s002) - qperpz * s022 + qperpx * s102 +
                     qperpz * s122 - 2 * q0x * s002 * theta) +
              q0z *
                  (2 * qperpz * s012 - qperpy * s022 + qperpw * (-s002 + s102) -
                   2 * qperpz * s112 + qperpy * s122 - 2 * q0w * s002 * theta -
                   2 * q0y * s022 * theta));

      c5[1] = DerivativeTerm(
          0.,
          -2 *
              (qperpx * qperpy * s000 + qperpw * qperpz * s000 +
               q0z * q0z * s010 - qperpx * qperpx * s010 -
               qperpz * qperpz * s010 - q0y * q0z * s020 -
               qperpw * qperpx * s020 + qperpy * qperpz * s020 -
               qperpx * qperpy * s100 - qperpw * qperpz * s100 +
               q0w * q0z * (-s000 + s100) + q0x * q0x * (s010 - s110) -
               q0z * q0z * s110 + qperpx * qperpx * s110 +
               qperpz * qperpz * s110 +
               q0x * (q0y * (-s000 + s100) + q0w * (s020 - s120)) +
               q0y * q0z * s120 + qperpw * qperpx * s120 -
               qperpy * qperpz * s120) *
              theta,
          -2 *
              (qperpx * qperpy * s001 + qperpw * qperpz * s001 +
               q0z * q0z * s011 - qperpx * qperpx * s011 -
               qperpz * qperpz * s011 - q0y * q0z * s021 -
               qperpw * qperpx * s021 + qperpy * qperpz * s021 -
               qperpx * qperpy * s101 - qperpw * qperpz * s101 +
               q0w * q0z * (-s001 + s101) + q0x * q0x * (s011 - s111) -
               q0z * q0z * s111 + qperpx * qperpx * s111 +
               qperpz * qperpz * s111 +
               q0x * (q0y * (-s001 + s101) + q0w * (s021 - s121)) +
               q0y * q0z * s121 + qperpw * qperpx * s121 -
               qperpy * qperpz * s121) *
              theta,
          -2 *
              (qperpx * qperpy * s002 + qperpw * qperpz * s002 +
               q0z * q0z * s012 - qperpx * qperpx * s012 -
               qperpz * qperpz * s012 - q0y * q0z * s022 -
               qperpw * qperpx * s022 + qperpy * qperpz * s022 -
               qperpx * qperpy * s102 - qperpw * qperpz * s102 +
               q0w * q0z * (-s002 + s102) + q0x * q0x * (s012 - s112) -
               q0z * q0z * s112 + qperpx * qperpx * s112 +
               qperpz * qperpz * s112 +
               q0x * (q0y * (-s002 + s102) + q0w * (s022 - s122)) +
               q0y * q0z * s122 + qperpw * qperpx * s122 -
               qperpy * qperpz * s122) *
              theta);

      c1[2] = DerivativeTerm(
          -t0z + t1z,
          (qperpw * qperpy * s000 - qperpx * qperpz * s000 - q0y * q0z * s010 -
           qperpw * qperpx * s010 - qperpy * qperpz * s010 - s020 +
           q0y * q0y * s020 + qperpx * qperpx * s020 + qperpy * qperpy * s020 -
           qperpw * qperpy * s100 + qperpx * qperpz * s100 +
           q0x * q0z * (-s000 + s100) + q0y * q0z * s110 +
           qperpw * qperpx * s110 + qperpy * qperpz * s110 +
           q0w * (q0y * (s000 - s100) + q0x * (-s010 + s110)) +
           q0x * q0x * (s020 - s120) + s120 - q0y * q0y * s120 -
           qperpx * qperpx * s120 - qperpy * qperpy * s120),
          (qperpw * qperpy * s001 - qperpx * qperpz * s001 - q0y * q0z * s011 -
           qperpw * qperpx * s011 - qperpy * qperpz * s011 - s021 +
           q0y * q0y * s021 + qperpx * qperpx * s021 + qperpy * qperpy * s021 -
           qperpw * qperpy * s101 + qperpx * qperpz * s101 +
           q0x * q0z * (-s001 + s101) + q0y * q0z * s111 +
           qperpw * qperpx * s111 + qperpy * qperpz * s111 +
           q0w * (q0y * (s001 - s101) + q0x * (-s011 + s111)) +
           q0x * q0x * (s021 - s121) + s121 - q0y * q0y * s121 -
           qperpx * qperpx * s121 - qperpy * qperpy * s121),
          (qperpw * qperpy * s002 - qperpx * qperpz * s002 - q0y * q0z * s012 -
           qperpw * qperpx * s012 - qperpy * qperpz * s012 - s022 +
           q0y * q0y * s022 + qperpx * qperpx * s022 + qperpy * qperpy * s022 -
           qperpw * qperpy * s102 + qperpx * qperpz * s102 +
           q0x * q0z * (-s002 + s102) + q0y * q0z * s112 +
           qperpw * qperpx * s112 + qperpy * qperpz * s112 +
           q0w * (q0y * (s002 - s102) + q0x * (-s012 + s112)) +
           q0x * q0x * (s022 - s122) + s122 - q0y * q0y * s122 -
           qperpx * qperpx * s122 - qperpy * qperpy * s122));

      c2[2] = DerivativeTerm(
          0.,
          (q0w * q0y * s000 - q0x * q0z * s000 - qperpw * qperpy * s000 +
           qperpx * qperpz * s000 - q0w * q0x * s010 - q0y * q0z * s010 +
           qperpw * qperpx * s010 + qperpy * qperpz * s010 + q0x * q0x * s020 +
           q0y * q0y * s020 - qperpx * qperpx * s020 - qperpy * qperpy * s020 -
           q0w * q0y * s100 + q0x * q0z * s100 + qperpw * qperpy * s100 -
           qperpx * qperpz * s100 + q0w * q0x * s110 + q0y * q0z * s110 -
           qperpw * qperpx * s110 - qperpy * qperpz * s110 - q0x * q0x * s120 -
           q0y * q0y * s120 + qperpx * qperpx * s120 + qperpy * qperpy * s120 -
           2 * q0y * qperpw * s000 * theta + 2 * q0z * qperpx * s000 * theta -
           2 * q0w * qperpy * s000 * theta + 2 * q0x * qperpz * s000 * theta +
           2 * q0x * qperpw * s010 * theta + 2 * q0w * qperpx * s010 * theta +
           2 * q0z * qperpy * s010 * theta + 2 * q0y * qperpz * s010 * theta -
           4 * q0x * qperpx * s020 * theta - 4 * q0y * qperpy * s020 * theta),
          (q0w * q0y * s001 - q0x * q0z * s001 - qperpw * qperpy * s001 +
           qperpx * qperpz * s001 - q0w * q0x * s011 - q0y * q0z * s011 +
           qperpw * qperpx * s011 + qperpy * qperpz * s011 + q0x * q0x * s021 +
           q0y * q0y * s021 - qperpx * qperpx * s021 - qperpy * qperpy * s021 -
           q0w * q0y * s101 + q0x * q0z * s101 + qperpw * qperpy * s101 -
           qperpx * qperpz * s101 + q0w * q0x * s111 + q0y * q0z * s111 -
           qperpw * qperpx * s111 - qperpy * qperpz * s111 - q0x * q0x * s121 -
           q0y * q0y * s121 + qperpx * qperpx * s121 + qperpy * qperpy * s121 -
           2 * q0y * qperpw * s001 * theta + 2 * q0z * qperpx * s001 * theta -
           2 * q0w * qperpy * s001 * theta + 2 * q0x * qperpz * s001 * theta +
           2 * q0x * qperpw * s011 * theta + 2 * q0w * qperpx * s011 * theta +
           2 * q0z * qperpy * s011 * theta + 2 * q0y * qperpz * s011 * theta -
           4 * q0x * qperpx * s021 * theta - 4 * q0y * qperpy * s021 * theta),
          (q0w * q0y * s002 - q0x * q0z * s002 - qperpw * qperpy * s002 +
           qperpx * qperpz * s002 - q0w * q0x * s012 - q0y * q0z * s012 +
           qperpw * qperpx * s012 + qperpy * qperpz * s012 + q0x * q0x * s022 +
           q0y * q0y * s022 - qperpx * qperpx * s022 - qperpy * qperpy * s022 -
           q0w * q0y * s102 + q0x * q0z * s102 + qperpw * qperpy * s102 -
           qperpx * qperpz * s102 + q0w * q0x * s112 + q0y * q0z * s112 -
           qperpw * qperpx * s112 - qperpy * qperpz * s112 - q0x * q0x * s122 -
           q0y * q0y * s122 + qperpx * qperpx * s122 + qperpy * qperpy * s122 -
           2 * q0y * qperpw * s002 * theta + 2 * q0z * qperpx * s002 * theta -
           2 * q0w * qperpy * s002 * theta + 2 * q0x * qperpz * s002 * theta +
           2 * q0x * qperpw * s012 * theta + 2 * q0w * qperpx * s012 * theta +
           2 * q0z * qperpy * s012 * theta + 2 * q0y * qperpz * s012 * theta -
           4 * q0x * qperpx * s022 * theta - 4 * q0y * qperpy * s022 * theta));

      c3[2] = DerivativeTerm(
          0.,
          -2 *
              (-(q0w * qperpy * s000) + q0x * qperpz * s000 +
               q0x * qperpw * s010 + q0w * qperpx * s010 -
               2 * q0x * qperpx * s020 + q0w * qperpy * s100 -
               q0x * qperpz * s100 - q0x * qperpw * s110 - q0w * qperpx * s110 +
               q0z * (qperpx * s000 + qperpy * s010 - qperpx * s100 -
                      qperpy * s110) +
               2 * q0x * qperpx * s120 +
               q0y * (qperpz * s010 - 2 * qperpy * s020 +
                      qperpw * (-s000 + s100) - qperpz * s110 +
                      2 * qperpy * s120)) *
              theta,
          -2 *
              (-(q0w * qperpy * s001) + q0x * qperpz * s001 +
               q0x * qperpw * s011 + q0w * qperpx * s011 -
               2 * q0x * qperpx * s021 + q0w * qperpy * s101 -
               q0x * qperpz * s101 - q0x * qperpw * s111 - q0w * qperpx * s111 +
               q0z * (qperpx * s001 + qperpy * s011 - qperpx * s101 -
                      qperpy * s111) +
               2 * q0x * qperpx * s121 +
               q0y * (qperpz * s011 - 2 * qperpy * s021 +
                      qperpw * (-s001 + s101) - qperpz * s111 +
                      2 * qperpy * s121)) *
              theta,
          -2 *
              (-(q0w * qperpy * s002) + q0x * qperpz * s002 +
               q0x * qperpw * s012 + q0w * qperpx * s012 -
               2 * q0x * qperpx * s022 + q0w * qperpy * s102 -
               q0x * qperpz * s102 - q0x * qperpw * s112 - q0w * qperpx * s112 +
               q0z * (qperpx * s002 + qperpy * s012 - qperpx * s102 -
                      qperpy * s112) +
               2 * q0x * qperpx * s122 +
               q0y * (qperpz * s012 - 2 * qperpy * s022 +
                      qperpw * (-s002 + s102) - qperpz * s112 +
                      2 * qperpy * s122)) *
              theta);

      c4[2] = DerivativeTerm(
          0.,
          q0w * qperpy * s000 - q0x * qperpz * s000 - q0x * qperpw * s010 -
              q0w * qperpx * s010 + 2 * q0x * qperpx * s020 -
              q0w * qperpy * s100 + q0x * qperpz * s100 + q0x * qperpw * s110 +
              q0w * qperpx * s110 - 2 * q0x * qperpx * s120 -
              2 * qperpw * qperpy * s000 * theta +
              2 * qperpx * qperpz * s000 * theta -
              2 * q0w * q0x * s010 * theta +
              2 * qperpw * qperpx * s010 * theta +
              2 * qperpy * qperpz * s010 * theta +
              2 * q0x * q0x * s020 * theta + 2 * q0y * q0y * s020 * theta -
              2 * qperpx * qperpx * s020 * theta -
              2 * qperpy * qperpy * s020 * theta +
              q0z * (-(qperpx * s000) - qperpy * s010 + qperpx * s100 +
                     qperpy * s110 - 2 * q0x * s000 * theta) +
              q0y *
                  (-(qperpz * s010) + 2 * qperpy * s020 +
                   qperpw * (s000 - s100) + qperpz * s110 - 2 * qperpy * s120 +
                   2 * q0w * s000 * theta - 2 * q0z * s010 * theta),
          q0w * qperpy * s001 - q0x * qperpz * s001 - q0x * qperpw * s011 -
              q0w * qperpx * s011 + 2 * q0x * qperpx * s021 -
              q0w * qperpy * s101 + q0x * qperpz * s101 + q0x * qperpw * s111 +
              q0w * qperpx * s111 - 2 * q0x * qperpx * s121 -
              2 * qperpw * qperpy * s001 * theta +
              2 * qperpx * qperpz * s001 * theta -
              2 * q0w * q0x * s011 * theta +
              2 * qperpw * qperpx * s011 * theta +
              2 * qperpy * qperpz * s011 * theta +
              2 * q0x * q0x * s021 * theta + 2 * q0y * q0y * s021 * theta -
              2 * qperpx * qperpx * s021 * theta -
              2 * qperpy * qperpy * s021 * theta +
              q0z * (-(qperpx * s001) - qperpy * s011 + qperpx * s101 +
                     qperpy * s111 - 2 * q0x * s001 * theta) +
              q0y *
                  (-(qperpz * s011) + 2 * qperpy * s021 +
                   qperpw * (s001 - s101) + qperpz * s111 - 2 * qperpy * s121 +
                   2 * q0w * s001 * theta - 2 * q0z * s011 * theta),
          q0w * qperpy * s002 - q0x * qperpz * s002 - q0x * qperpw * s012 -
              q0w * qperpx * s012 + 2 * q0x * qperpx * s022 -
              q0w * qperpy * s102 + q0x * qperpz * s102 + q0x * qperpw * s112 +
              q0w * qperpx * s112 - 2 * q0x * qperpx * s122 -
              2 * qperpw * qperpy * s002 * theta +
              2 * qperpx * qperpz * s002 * theta -
              2 * q0w * q0x * s012 * theta +
              2 * qperpw * qperpx * s012 * theta +
              2 * qperpy * qperpz * s012 * theta +
              2 * q0x * q0x * s022 * theta + 2 * q0y * q0y * s022 * theta -
              2 * qperpx * qperpx * s022 * theta -
              2 * qperpy * qperpy * s022 * theta +
              q0z * (-(qperpx * s002) - qperpy * s012 + qperpx * s102 +
                     qperpy * s112 - 2 * q0x * s002 * theta) +
              q0y *
                  (-(qperpz * s012) + 2 * qperpy * s022 +
                   qperpw * (s002 - s102) + qperpz * s112 - 2 * qperpy * s122 +
                   2 * q0w * s002 * theta - 2 * q0z * s012 * theta));

      c5[2] = DerivativeTerm(
          0.,
          2 *
              (qperpw * qperpy * s000 - qperpx * qperpz * s000 +
               q0y * q0z * s010 - qperpw * qperpx * s010 -
               qperpy * qperpz * s010 - q0y * q0y * s020 +
               qperpx * qperpx * s020 + qperpy * qperpy * s020 +
               q0x * q0z * (s000 - s100) - qperpw * qperpy * s100 +
               qperpx * qperpz * s100 +
               q0w * (q0y * (-s000 + s100) + q0x * (s010 - s110)) -
               q0y * q0z * s110 + qperpw * qperpx * s110 +
               qperpy * qperpz * s110 + q0y * q0y * s120 -
               qperpx * qperpx * s120 - qperpy * qperpy * s120 +
               q0x * q0x * (-s020 + s120)) *
              theta,
          2 *
              (qperpw * qperpy * s001 - qperpx * qperpz * s001 +
               q0y * q0z * s011 - qperpw * qperpx * s011 -
               qperpy * qperpz * s011 - q0y * q0y * s021 +
               qperpx * qperpx * s021 + qperpy * qperpy * s021 +
               q0x * q0z * (s001 - s101) - qperpw * qperpy * s101 +
               qperpx * qperpz * s101 +
               q0w * (q0y * (-s001 + s101) + q0x * (s011 - s111)) -
               q0y * q0z * s111 + qperpw * qperpx * s111 +
               qperpy * qperpz * s111 + q0y * q0y * s121 -
               qperpx * qperpx * s121 - qperpy * qperpy * s121 +
               q0x * q0x * (-s021 + s121)) *
              theta,
          2 *
              (qperpw * qperpy * s002 - qperpx * qperpz * s002 +
               q0y * q0z * s012 - qperpw * qperpx * s012 -
               qperpy * qperpz * s012 - q0y * q0y * s022 +
               qperpx * qperpx * s022 + qperpy * qperpy * s022 +
               q0x * q0z * (s002 - s102) - qperpw * qperpy * s102 +
               qperpx * qperpz * s102 +
               q0w * (q0y * (-s002 + s102) + q0x * (s012 - s112)) -
               q0y * q0z * s112 + qperpw * qperpx * s112 +
               qperpy * qperpz * s112 + q0y * q0y * s122 -
               qperpx * qperpx * s122 - qperpy * qperpy * s122 +
               q0x * q0x * (-s022 + s122)) *
              theta);
    }
  }

  ~AnimationTransform(){};

  void Decompose(const Transform *t, Vector3f *T, Quaternion *R, Matrix4f *S);

  void Interpolate(Float time, Transform *t) const;

  Bounds3f MotionBounds(const Bounds3f &b) const;

  Bounds3f BoundPointMotion(const Point3f &p) const;

  Point3f operator()(Float time, const Point3f &p);

  Vector3f operator()(Float time, const Vector3f &v);

  Ray operator()(const Ray &r);

  RayDifferential operator()(const RayDifferential &r);
};

};  // namespace eric