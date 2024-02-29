#include "bounds.hpp"
#include "ray.hpp"
#include "intervals.hpp"

namespace eric {
/*
-----------------------------------------------------------------------------------------
** Bounds class
----------------------------------------------------------------------------------------- 
*/  
// Bounds2
template <typename Ty>
const Point2<Ty> Bounds2<Ty>::operator[](int i) const {
    assert(i <= 3 && i >= 0);
    Ty x = (i & 1) == 0 ? pMin.x : pMax.x;
    Ty y = (i & 2) == 0 ? pMin.y : pMax.y;
    return Point2<Ty>(x, y);
};

template <typename Ty>
Point2<Ty> Bounds2<Ty>::operator[](int i) {
    assert(i <= 3 && i >= 0);
    Ty x = (i & 1) == 0 ? pMin.x : pMax.x;
    Ty y = (i & 2) == 0 ? pMin.y : pMax.y;
    return Point2<Ty>(x, y);
};

template <typename Ty>
Point2<Ty> Bounds2<Ty>::_corner(int i) const {
    return (*this)[i];
};

template <typename Ty>
Bounds2<Ty> Bounds2<Ty>::_union(const Point2<Ty> &p) const {
    return Bounds2<Ty>(
            Point2<Ty>(
                std::min(pMin.x, p.x),
                std::min(pMin.y, p.y)
            ),
            Point2<Ty>(
                std::max(pMax.x, p.x),
                std::max(pMax.y, p.y)
            )
        );
};

template <typename T>
Vector2<T> Bounds2<T>::_diagonal() const{
    return pMax - pMin;
};

//Bounds3
template <typename Ty>
const Point3<Ty> Bounds3<Ty>::operator[](int i) const {
    assert(i <= 7 && i >= 0);
    return Point3<Ty>(
        (i & 1) == 0 ? pMin.x : pMax.x,
        (i & 2) == 0 ? pMin.y : pMax.y,
        (i & 4) == 0 ? pMin.z : pMax.z
    );
};

template <typename Ty>
Point3<Ty> Bounds3<Ty>::operator[](int i) {
    assert(i <= 7 && i >= 0);
    return Point3<Ty>(
        (i & 1) == 0 ? pMin.x : pMax.x,
        (i & 2) == 0 ? pMin.y : pMax.y,
        (i & 4) == 0 ? pMin.z : pMax.z
    );
};

template <typename Ty>
Point3<Ty> Bounds3<Ty>::_corner(int i) const {
    return (*this)[i];
};

template <typename Ty>
Bounds3<Ty> Bounds3<Ty>::_union(const Point3<Ty> &p) const {
    return Bounds3<Ty>(
            Point3<Ty>(
                std::min(pMin.x, p.x),
                std::min(pMin.y, p.y),
                std::min(pMin.z, p.z)
            ),
            Point3<Ty>(
                std::max(pMax.x, p.x),
                std::max(pMax.y, p.y),
                std::max(pMax.z, p.z)
            )
        );
};


template <typename T>
Vector3<T> Bounds3<T>::_diagonal() const {
    return pMax - pMin;
};

template <typename T>
T Bounds3<T>::_surface_area() const {
    Vector3<T> diagonal = _diagonal();
    return 2 * (diagonal.x * diagonal.y + diagonal.y * diagonal.z + diagonal.z * diagonal.x);
};

template <typename T>
T Bounds3<T>::_volume() const {
    Vector3<T> diagonal = _diagonal();
    return diagonal.x * diagonal.y * diagonal.z;
};

template <typename T>
int Bounds3<T>::_maximum_extent() const {
    Vector3<T> diagonal = _diagonal();
    return (diagonal.x > diagonal.y ? 
        (diagonal.x > diagonal.z ? 0 : 2) : (diagonal.y > diagonal.z ? 1 : 2));
};

template <typename T>
Point3<T> Bounds3<T>::_lerp(Point3<T> &t) const {
    return Point3<T>(
        ::Lerp(t.x, pMin.x, pMax.x),
        ::Lerp(t.y, pMin.y, pMax.y),
        ::Lerp(t.z, pMin.z, pMax.z)
    ); 
};

template <typename T>
Vector3<T> Bounds3<T>::_offset(Point3<T> &p) const {
    Vector3<T> o = p - pMin;
    if (pMax.x > pMin.x) o.x /= pMax.x - pMin.x;
    if (pMax.y > pMin.y) o.y /= pMax.y - pMin.y;
    if (pMax.z > pMin.z) o.z /= pMax.z - pMin.z;

    return o;
};

template <typename T>
void Bounds3<T>::_bounding_sphere(Point3<T> *center, T *radius) const {
    *center = static_cast<Point3<T>>(pMin + (pMax - pMin));
    *radius = Inside(*center, *this) ? Distance(pMax, *center) : 0;
}; 

template <typename T>
bool Bounds3<T>::_intersect_p(const Ray& r, Float *hitt0, Float *hitt1) const {
    assert(hitt0 && hitt1);
    Interval it(0, std::numeric_limits<T>::max());
    for (int i = 0; i < 3; i++) {
        Float t1, t2;
        if (r.d[i] == 0) {
            t1 = t2 = std::numeric_limits<T>::max();
        } else {
            t1 = (pMin[i] - r.o[i]) / r.d[i];
            t2 = (pMax[i] - r.o[i]) / r.d[i];
        }
        if (t1 > t2) {
            Float tmp = t2;
            t2 = t1;
            t1 = tmp;
        }
        it = it._intersect(Interval(t1, t2)); 
    }
    if (!it._degenerated()) {
        *hitt0 = it.low;
        *hitt1 = it.high;
        return true;
    } else {
        return false;
    }
};

template <typename T>
bool Bounds3<T>::_intersect_p(const Ray& r, const Vector3<T>& invDir, const int dirIsNeg[3]) const {
    Float t1, t2;
    const Bounds3<T>& b = *this;
    Interval it(std::numeric_limits<T>::min(), std::numeric_limits<T>::max());

    for (int i = 0; i < 3; i++) {
        t1 = (b[dirIsNeg[i] << i][i] - r.o[i]) * invDir[i];
        t2 = (b[(1 - dirIsNeg[i]) << i][i] - r.o[i]) * invDir[i];
        it = it._intersect(Interval(t1, t2));
    }
    if (it._degenerated()) {
        return false;
    } else {
        return true;
    }
};

};