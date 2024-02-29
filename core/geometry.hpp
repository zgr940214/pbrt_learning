#pragma once
#include "math.hpp"

#include <cassert>
#include <cmath>
#include <algorithm>
#include <limits>

namespace eric {

double Radians(double degrees) {
    return degrees * M_PI / 180.0;
}

Float Radians(Float degrees) {
    return degrees * M_PI / 180.0;
}

Float Clamp(Float v, Float min, Float max) {
    v = v > min ? v : min;
    v = v < max ? v : max;
    return v;
}

/*
-----------------------------------------------------------------------------------------
** Vector:
----------------------------------------------------------------------------------------- 
*/
template <typename Ty>
class Normal3;

template <typename Ty>
class Vector2 {
    public:
        union {
            struct {
                Ty x = 0;
                Ty y = 0;
            };
            struct {
                Ty r;
                Ty s;
            };
            struct {
                Ty u;
                Ty v;
            };
        };

    public: 
        Vector2()=default;
        Vector2(Ty x, Ty y): x(x), y(y){assert(!has_nans())};
        template<typename U>
        explicit Vector2(const Vector3<U> &v):x(static_cast<Ty>(v.x)), y(static_cast<Ty>(p.y)){assert(!has_nans());};

        Ty dot(const Vector2<Ty> &v) const;
        Vector2<Ty>& abs();

        Ty length() const;
        Ty length_squared() const;

        Vector2<Ty>& normalize();
        Vector2<Ty> normalized() const;
        
        bool has_nans()const;

        Ty operator[](int index) const;
        Vector2<Ty> operator+(const Vector2<Ty> &v) const;
        Vector2<Ty> & operator+=(const Vector2<Ty> &v);
        Vector2<Ty> operator-(const Vector2<Ty> &v) const;
        Vector2<Ty> & operator-=(const Vector2<Ty> &v);
        Vector2<Ty> operator*(Ty v) const;
        Vector2<Ty>& operator*=(Ty v); 
        Vector2<Ty> operator/(Ty v) const;
        Vector2<Ty>& operator/=(Ty v);
        Vector2<Ty>& operator-(); // unary operator -, negated all val;
};

template <typename Ty>
class Vector3 {
    public:
        union {
            struct {
                Ty x;
                Ty y;
                Ty z;
            };
            struct {
                Ty r;
                Ty s;
                Ty t;
            };
        };

    public:
        Vector3():x(0), y(0), z(0){};
        Vector3(Ty x, Ty y, Ty z):x(x), y(y), z(z){assert(!has_nans())};
        explicit Vector3(const Normal3<Ty> &n);

        Ty dot(const Vector3<Ty> &v) const;
        Vector3<Ty> cross(const Vector3<Ty> &v) const;
        Vector3<Ty>& abs();

        Ty length() const;
        Ty length_squared() const;

        Vector3<Ty> normalized() const;
        Vector3<Ty>& normalize();
        
        bool has_nans() const;

        explicit operator Normal3<Ty>() const;

        Ty operator [](int index) const;
        Vector3<Ty> operator+(const Vector3<Ty> &v) const;
        Vector3<Ty> & operator+=(const Vector3<Ty> &v);
        Vector3<Ty> operator-(const Vector3<Ty> &v) const;
        Vector3<Ty> & operator-=(const Vector3<Ty> &v);
        Vector3<Ty> operator*(Ty v) const;
        Vector3<Ty>& operator*=(Ty v); 
        Vector3<Ty> operator/(Ty v) const;
        Vector3<Ty>& operator/=(Ty v);
        Vector3<Ty>& operator-(); // unary operator -, negated all val;
    
};

/*
 ** Non-Member inline tool Function for Vector calculation 
*/

template<typename Ty> inline Vector2<Ty> 
operator*(Ty s, const Vector2<Ty> &v) {
    return v * s;
};

template<typename Ty> inline Vector3<Ty> 
operator*(Ty s, const Vector3<Ty> &v) {
    return v * s;
};

template<typename Ty> inline Vector2<Ty> 
Abs(const Vector2<Ty> &v) {
    return Vector2<Ty>{std::abs(v.x), std::abs(v.y)};
};

template<typename Ty> inline Vector3<Ty> 
Abs(const Vector3<Ty> &v) {
    return Vector3<Ty>{std::abs(v.x), std::abs(v.y), std::abs(v.z)};
};

template<typename Ty> inline Ty
Dot(const Vector2<Ty> &v1, const Vector2<Ty> &v2) {
    return v1.dot(v2);
};


template<typename Ty> inline Ty
Dot(const Vector3<Ty> &v1, const Vector3<Ty> &v2) {
    return v1.dot(v2);
};

template<typename Ty> inline Vector3<Ty>
Cross(const Vector3<Ty> &v1, const Vector3<Ty> &v2) {
    return v1.cross(v2);
}; 

template<typename Ty> inline Ty
AbsDot(const Vector2<Ty>& v1, const Vector2<Ty> &v2) {
    return std::abs(Dot(v1, v2));
};

template<typename Ty> inline Ty
AbsDot(const Vector3<Ty>& v1, const Vector3<Ty> &v2) {
    return std::abs(Dot(v1, v2));
};

template<typename Ty> inline Vector2<Ty>
Normalize (const Vector2<Ty> &v) {
    return v.normalized();
};

template<typename Ty> inline Vector3<Ty>
Normalize (const Vector3<Ty> &v) {
    return v.normalized();
};

template<typename Ty> inline Ty
MinComponent(const Vector2<Ty> &v) {
    return v.x < v.y ? v.x : v.y;
};

template<typename Ty> inline Ty
MinComponent(const Vector3<Ty> &v) {
    return (v.x < v.y) ? (v.x < v.z ? v.x : v.z) : (v.y < v.z ? v.y : v.z);
};

template<typename Ty> inline Ty
MaxComponent(const Vector2<Ty>) {
    return v.x > v.y ? v.x : v.y;
};

template<typename Ty> inline Ty
MaxComponent(const Vector3<Ty> &v) {
    return (v.x > v.y) ? (v.x > v.z ? v.x : v.z) :(v.y > v.z ? v.y : v.z);
};

template<typename Ty> inline int
MinComponentIndex(const Vector2<Ty> &v) {
    return (v.x < v.y) ? 0 : 1;
};

template<typename Ty> inline int
MinComponentIndex(const Vector3<Ty> &v) {
    return (v.x < v.y) ? (v.x < v.z ? 0 : 2) :(v.y < v.z ? 1 : 2);
};

template<typename Ty> inline int
MaxComponentIndex(const Vector2<Ty> &v) {
    return (v.x > v.y) ? 0 : 1;
};

template<typename Ty> inline int
MaxComponentIndex(const Vector3<Ty> &v) {
    return (v.x > v.y) ? (v.x > v.z ? 0 : 2) :(v.y > v.z ? 1 : 2);
};

template<typename Ty> inline Vector2<Ty>
Min (const Vector2<Ty> &v1, const Vector2<Ty> &v2) {
    return Vector2<Ty> {
        std::min(v1.x, v2.x),
        std::min(v1.y, v2.y)};
};

template<typename Ty> inline Vector2<Ty>
Min (const Vector3<Ty> &v1, const Vector3<Ty> &v2) {
    return Vector3<Ty> {
        std::min(v1.x, v2.x),
        std::min(v1.y, v2.y),
        std::min(v1.z, v2.z)};
};

template<typename Ty> inline Vector2<Ty>
Max (const Vector2<Ty> &v1, const Vector2<Ty> &v2) {
    return Vector2<Ty> {
        std::max(v1.x, v2.x),
        std::max(v1.y, v2.y)};
};

template<typename Ty> inline Vector2<Ty>
Min (const Vector3<Ty> &v1, const Vector3<Ty> &v2) {
    return Vector3<Ty> {
        std::max(v1.x, v2.x),
        std::max(v1.y, v2.y),
        std::max(v1.z, v2.z)};
};

template <typename Ty> inline void 
CoordinateSystem (const Vector3<Ty> &v1, Vector3<Ty> *v2, Vector3<Ty> *v3) {
    if (std::abs(v1.x) > std::abs(v1.y)) {
        *v2 = Vector2<Ty> {-v1.z, 0, v1.x};
    } else {
        *v2 = Vector2<Ty> {0, v1.z, -v1.y};
    }

    *v3 = Cross(v1, *v2);
};

template <typename Ty> inline Vector2<Ty>
Faceforward(const Vector2<Ty> &v1, const Vector2<Ty> &v2) {
return (Dot(v1, v2) < 0.f) ? -v1 : v1;
};

template <typename Ty> inline Vector3<Ty>
Faceforward(const Vector3<Ty> &v1, const Vector3<Ty> &v2) {
return (Dot(v1, v2) < 0.f) ? -v1 : v1;
};

using Vector2f = Vector2<Float>;
using Vector2i = Vector2<int>;
using Vector3f = Vector3<Float>;
using Vector3i = Vector3<int>;

/*
-----------------------------------------------------------------------------------------
** Point: represent a position in coordinate system
----------------------------------------------------------------------------------------- 
*/

template<typename Ty>
class Point2 {
    public:
        union {
            struct {
                Ty x;
                Ty y;
            };
            struct {
                Ty s;
                Ty r;
            };
            struct {
                Ty u;
                Ty v;
            }
        };

    public:
        Point2():x(0), y(0){};
        Point2(Ty x, Ty y):x(x), y(y) {assert(!has_nans());};
        template<typename U>
        explicit Point2(const Point3<U> &p):x(static_cast<Ty>(p.x)), y(static_cast<Ty>(p.x)){assert(!has_nans())};
        
        bool has_nans() const;

        Ty length() const;
        Ty length_squared() const;


        explicit operator Vector2<Ty>();
        explicit operator Vector3<Ty>();

        Point2<Ty> operator+(const Point2<Ty> &v) const;
        Point2<Ty> & operator+=(const Point2<Ty> &v);
        Point2<Ty> operator-(const Vector2<Ty> &v) const;
        Vector2<Ty> operator-(const Point2<Ty> &v) const;
        Point2<Ty> & operator-=(const Point2<Ty> &v);
        Point2<Ty> operator*(Ty v) const;
        Point2<Ty>& operator*=(Ty v); 
        Point2<Ty> operator/(Ty v) const;
        Point2<Ty>& operator/=(Ty v);
        Point2<Ty>& operator-(); // unary operator -, negated all val;
};

template<typename Ty>
class Point3 {
    public:
        union {
            struct {
                Ty x;
                Ty y;
                Ty z;
            };
            struct {
                Ty s;
                Ty r;
                Ty t;
            };
        };

    public:
        Point3():x(0), y(0), z(0){};
        Point3(Ty x, Ty y, Ty z):x(x), y(y), z(z){assert(!has_nans());};

        bool has_nans() const;

        Ty length() const;
        Ty length_squared() const;


        explicit operator Vector3<Ty>();
        Point3<Ty> operator+(const Vector3<Ty> &v) const;
        Point3<Ty>& operator+=(const Point3<Ty> &p);
        Point3<Ty>& operator+=(const Vector3<Ty> &v);
        Point3<Ty> operator+(const Point3<Ty> &p) const;
        Point3<Ty> operator-(const Vector3<Ty> &v) const;
        Vector3<Ty> operator-(const Point3<Ty> &v) const;
        Point3<Ty>& operator-=(const Point3<Ty> &v);
        Point3<Ty> operator*(Ty v) const;
        Point3<Ty>& operator*=(Ty v); 
        Point3<Ty> operator/(Ty v) const;
        Point3<Ty>& operator/=(Ty v);
        Point3<Ty>& operator-(); // unary operator -, negated all val;
        Ty operator[](int index) const;
};


/*
 ** Non-Member inline Tool Function 
*/

template<typename Ty> inline Point2<Ty>
operator*(Float s, const Point2<Ty>& p) {
    return p * s;
}

template<typename Ty> inline Point3<Ty>
operator*(Float s, const Point3<Ty>& p) {
    return p * s;
}

template<typename Ty> inline Ty
Distance(const Point2<Ty> &p1, const Point2<Ty> &p2) {
    return std::sqrt(DistanceSquared(p1, p2));
};

template<typename Ty> inline Ty
DistanceSquared(const Point2<Ty> &p1, const Point2<Ty> &p2) {
    return (p1 - p2).length_squared();
};

template<typename Ty> inline Ty
Distance(const Point3<Ty> &p1, const Point3<Ty> &p2) {
    return std::sqrt(DistanceSquared(p1, p2));
};

template<typename Ty> inline Ty
DistanceSquared(const Point3<Ty> &p1, const Point3<Ty> &p2) {
    return (p1 - p2).length_squared();
};

template<typename Ty> inline Point2<Ty>
Lerp(Float t, const Point2<Ty> &p1, const Point2<Ty> &p2) {
    assert(t <= 1 && t >= 0);
    return (1 - t) * p1 + t * p2;
};

template<typename Ty> inline Point3<Ty>
Lerp(Float t, const Point3<Ty> &p1, const Point3<Ty> &p2) {
    assert(t <= 1 && t >= 0);
    return (1 - t) * p1 + t * p2;
};

template <typename Ty> inline Ty 
Lerp(Float t, const Ty v1, const Ty v2) {
    assert(t <= 1 && t >= 0);
    return (1 - t) * v1 + t * v2;
}

template <typename Ty> Point3<Ty>
Min(const Point3<Ty> &p1, const Point3<Ty> &p2) {
return Point3<Ty>(std::min(p1.x, p2.x), std::min(p1.y, p2.y),
std::min(p1.z, p2.z));
};

template <typename Ty> Point3<Ty>
Max(const Point3<Ty> &p1, const Point3<Ty> &p2) {
return Point3<Ty>(std::max(p1.x, p2.x), std::max(p1.y, p2.y),
std::max(p1.z, p2.z));
};

template <typename Ty> Point3<Ty> Floor(const Point3<Ty> &p) {
return Point3<T>(std::floor(p.x), std::floor(p.y), std::floor(p.z));
};

template <typename Ty> Point3<Ty> Ceil(const Point3<Ty> &p) {
return Point3<Ty>(std::ceil(p.x), std::ceil(p.y), std::ceil(p.z));
};

template <typename Ty> Point3<Ty> Abs(const Point3<Ty> &p) {
return Point3<Ty>(std::abs(p.x), std::abs(p.y), std::abs(p.z));
};

template <typename Ty> Point3<Ty> Permute(const Point3<Ty> &p, int x, int y, int z) {
    return Point3<Ty>(p[x], p[y], p[z]);
};

using Point2f = Point2<Float>;
using Point2i = Point2<int>;
using Point3f = Point3<Float>;
using Point3i = Point3<int>;

/*
-----------------------------------------------------------------------------------------
** Normal: perpendicular to a surface, similar to vector3
----------------------------------------------------------------------------------------- 
*/
template <typename Ty> 
class Normal3 {
    public:
        union {
            struct {
                Ty x;
                Ty y;
                Ty z;
            };
        };

    public:
        Normal3() = default;
        Normal3(Ty x, Ty y, Ty z):x(x), y(y), z(z){assert(!has_nans())};
        explicit Normal3(const Vector3<Ty> &v):x(v.x), y(v.y), z(v.z) {assert(!has_nans())};

        bool has_nans() const;

        Ty dot(const Vector3<Ty> &v) const;
        Vector3<Ty> cross(const Vector3<Ty> &v) const;

        Normal3<Ty>& scale(Ty s);
        Normal3<Ty> scaled(Ty s) const;

        Normal3<Ty>& normalize();
        Normal3<Ty> normalized() const;

        Ty length() const;
        Ty length_squared() const;


        explicit operator Vector3<Ty>() {
            return Vector3<Ty>(x, y, z);
        };

        Vector3<Ty> operator+(const Vector3<Ty> &v) const {
            Vector3<Ty> m(x, y, z);
            m += v;
            return m;
        };

        Normal3<Ty> & operator+=(const Vector3<Ty> &v) {
            x += v.x;
            y += v.y;
            z += v.z;
            return *this;
        };

        Vector3<Ty> operator-(const Vector3<Ty> &v) const {
            Vector3<Ty> m(x, y, z);
            m -= v;
            return m;
        };

        Normal3<Ty> & operator-=(const Vector3<Ty> &v) {
            x -= v.x;
            y -= v.y;
            z -= v.z;
            return *this;
        };

        Normal3<Ty> operator*(Ty v) const {
            Normal3<Ty> m(*this);
            m *= v;
            return m;
        };

        Normal3<Ty>& operator*=(Ty v) {
            x *= v;
            y *= v;
            z *= v;
            return *this;
        }; 

        Normal3<Ty> operator/(Ty v) const {
            Normal3<Ty> m(*this);
            m /= v;
            return m;
        };

        Normal3<Ty>& operator/=(Ty v) {
            x /= v;
            y /= v;
            z /= v;
            return *this;
        };

        Normal3<Ty>& operator-() {// unary operator -, negated all val;
            x = -x;
            y = -y;
            z = -z;
            return *this;
        }; 

        bool operator==(const Normal3<Ty>& n) const {
            Ty dx = x - n.x, dy = y - n.y, dz = z - n.z;
            #define NOT_ZERO(x) (x > 0.001f || x < -0.001f)
                return !NOT_ZERO(dx) && !NOT_ZERO(dy) && !NOT_ZERO(dz);
            #undef NOT_ZERO 
        };

        bool operator!=(const Normal3<Ty>& n) const {
            Ty dx = x - n.x, dy = y - n.y, dz = z - n.z;
            #define NOT_ZERO(x) (x > 0.001f || x < -0.001f)
                return NOT_ZERO(dx) || NOT_ZERO(dy) || NOT_ZERO(dz);
            #undef NOT_ZERO 
        };
};

/*
 ** Non-Member inline Tool Function 
*/

template <typename Ty> inline Normal3<Ty>
operator* (Ty v, const Normal3<Ty> &n) {
    return Normal3<Ty>(n.x * v, n.y * v, n.z * v);
};

template <typename Ty> inline Normal3<Ty> 
Abs(const Normal3<Ty> &v) {
    return Normal3<Ty>{std::abs(v.x), std::abs(v.y), std::abs(v.z)};
};

template <typename Ty> inline Ty
Dot(const Normal3<Ty> &n, const Vector3<Ty> &v) {
    return n.dot(v);
};

template <typename Ty> inline Ty
Dot(const Vector3<Ty> &v, const Normal3<Ty> &n) {
    return v.dot(n);
};

template <typename Ty> inline Vector3<Ty>
Cross(const Normal3<Ty> &n, const Vector3<Ty> &v) {
    return n.cross(v);
}; 

template <typename Ty> inline Vector3<Ty>
Cross(const Vector3<Ty> &v, const Normal3<Ty> &n) {
    return v.cross(n);
}; 

template <typename Ty> inline Ty
AbsDot(const Normal3<Ty>& n, const Vector3<Ty> &v) {
    return std::abs(Dot(n, v));
};

template <typename Ty> inline Ty
AbsDot(const Vector3<Ty> &v, const Normal3<Ty>& n) {
    return std::abs(Dot(v, n));
};

template <typename Ty> inline Normal3<Ty>
Normalize (const Normal3<Ty> &n) {
    return n.normalized();
};

template <typename Ty> inline Ty
MinComponent(const Normal3<Ty> &n) {
    return (n.x < n.y) ? (n.x < n.z ? n.x : n.z) : (n.y < n.z ? n.y : n.z);
};

template <typename Ty> inline Ty
MaxComponent(const Normal3<Ty> &n) {
    return (n.x > n.y) ? (n.x > n.z ? n.x : n.z) :(n.y > n.z ? n.y : n.z);
};

template <typename Ty> inline int
MinComponentIndex(const Normal3<Ty> &n) {
    return (n.x < n.y) ? (n.x < n.z ? 0 : 2) :(n.y < n.z ? 1 : 2);
};

template <typename Ty> inline int
MaxComponentIndex(const Normal3<Ty> &n) {
    return (n.x > n.y) ? (n.x > n.z ? 0 : 2) :(n.y > n.z ? 1 : 2);
};

template <typename T> inline Normal3<T>
Faceforward(const Normal3<T> &n, const Vector3<T> &v) {
return (Dot(n, v) < 0.f) ? -n : n;
}


using Normal3f = Normal3<Float>;


}//eric