#pragma once
#include "geometry.hpp"

namespace eric {
/*
-----------------------------------------------------------------------------------------
** Bounds: for collision test , in pbr we use AABB axis aligned bounding box
** for AABB we have two options to represent it . 1: one of its vertice along with 3 length
** 2: two opposite vertices , here we use option 2
----------------------------------------------------------------------------------------- 
*/  
template <typename Ty>
class Bounds2 {
    public:
        Point2<Ty> pMin;
        Point2<Ty> pMax;
    public:
        Bounds2() {
            Ty min = std::numeric_limits<Ty>::lowest();
            Ty max = std::numeric_limits<Ty>::max();
            pMin = Point2<Ty>(min, min);
            pMax = Point2<Ty>(max, max);
        };
        
        Bounds2(const Point2<Ty> &p):pMin(p), pMax(p){};
        
        /// @brief p1 p2 represent 2 opposite corners of the bounding box
        Bounds2(const Point2<Ty> &p1, const Point2<Ty> &p2) { 
            pMin = Point2<Ty>(std::min(p1.x, p2.x), std::min(p1.y, p2.y));
            pMax = Point2<Ty>(std::max(p1.x, p2.x), std::max(p1.y, p2.y));
        };

        /// @brief return the one of four corners of the bounding box based on index 
        /// @param i first bit select the x , second select the y
        /// @return corner point
        Point2<Ty> _corner(int i) const;
        /// @brief create a new bounds2 with a point2 , if the points is outside the bounds volume, then return the bounds, else set the new bound according to point
        /// @param p new bound point
        /// @return new Bounds2
        Bounds2<Ty> _union(const Point2<Ty> &p) const;

        Vector2<Ty> _diagonal() const;

        const Point2<Ty> operator[](int i) const;
        Point2<Ty> operator[](int i);

};

class Ray;
template <typename Ty>
class Bounds3 {
    public:
        Point3<Ty> pMin;
        Point3<Ty> pMax;
    public:
        Bounds3() {
            Ty min = std::numeric_limits<Ty>::lowest();
            Ty max = std::numeric_limits<Ty>::max();
            pMin = Point3<Ty>(min, min, min);
            pMax = Point3<Ty>(max, max, max);
        };
        
        Bounds3(const Point3<Ty> &p):pMin(p), pMax(p){};
        
        /// @brief p1 p2 represent 2 opposite corners of the bounding box
        Bounds3(const Point3<Ty> &p1, const Point3<Ty> &p2) { 
            pMin = Point3<Ty>(std::min(p1.x, p2.x), std::min(p1.y, p2.y), std::min(p1.z, p2.z));
            pMax = Point3<Ty>(std::max(p1.x, p2.x), std::max(p1.y, p2.y), std::max(p1.z, p2.z));
        };

        /// @brief return the one of four corners of the bounding box based on index 
        /// @param i first bit select the x , second select the y, third select the z
        /// @return corner point
        Point3<Ty> _corner(int i) const;
        /// @brief create a new bounds3 with a point3 , if the points is outside the bounds volume, then return the bounds, else set the new bound according to point
        /// @param p new bound point
        /// @return new Bounds3
        Bounds3<Ty> _union(const Point3<Ty> &p) const;

        Vector3<Ty> _diagonal() const;

        Ty _surface_area() const;

        Ty _volume() const;

        int _maximum_extent() const;

        Point3<Ty> _lerp(Point3<Ty> &t) const;

        Vector3<Ty> _offset(Point3<Ty> &p) const;

        void _bounding_sphere(Point3<Ty> *center, Ty *radius) const;

        bool _intersect_p(const Ray& r, Float *hitt0, Float *hitt1) const;

        bool _intersect_p(const Ray& r, const Vector3<Ty>& invDir, const int dirIsNeg[3]) const;

        const Point3<Ty> operator[](int i) const;

        Point3<Ty> operator[](int i);
};

/*
 ** Non-Member inline Tool Function 
*/
template <typename T> inline Bounds3 <T>
Union(const Bounds3<T> &b, const Point3<T> &p) {
    return b._union(p);
}

template <typename T> inline Bounds3 <T>
Union(const Bounds3<T> &b1, const Bounds3<T> &b2) {
    return Bounds3<T>(
        Point3<T>(
            std::min(b1.pMin.x, b2.pMin.x),
            std::min(b1.pMin.y, b2.pMin.y),
            std::min(b1.pMin.z, b2.pMin.z)
        ),
        Point3<T>(
            std::max(b1.pMax.x, b2.pMax.x),
            std::max(b1.pMax.y, b2.pMax.y),
            std::max(b1.pMax.z, b2.pMax.z)
        )
    ); 
}

template <typename T> inline Bounds3 <T> 
Intersect(const Bounds3<T> &b1, const Bounds3<T> &b2) {
    return Bounds3<T>(
        Point3<T>(
            std::max(b1.pMin.x, b2.pMin.x),
            std::max(b1.pMin.y, b2.pMin.y),
            std::max(b1.pMin.z, b2.pMin.z)
        ),
        Point3<T>(
            std::min(b1.pMax.x, b2.pMax.x),
            std::min(b1.pMax.y, b2.pMax.y),
            std::min(b1.pMax.z, b2.pMax.z)
        )
    ); 
};

template <typename T> inline
bool Overlaps(const Bounds3<T> &b1, const Bounds3<T> &b2) {
    bool x = (b1.pMax.x > b2.pMin.x) && (b2.pMax.x > b1.pMin.x);
    bool y = (b1.pMax.y > b2.pMin.y) && (b2.pMax.y > b1.pMin.y);
    bool z = (b1.pMax.z > b2.pMin.z) && (b2.pMax.z > b1.pMin.z);
    return (x && y && z);
};

template <typename T> inline
bool Inside(const Point3<T> &p, const Bounds3<T> &b) {
    bool x = (p.x >= b.pMin.x) && (p.x <= b.pMax.x);
    bool y = (p.y >= b.pMin.y) && (p.y <= b.pMax.y);
    bool z = (p.z >= b.pMin.z) && (p.z <= b.pMax.z);
    return (x && y && z); 
};

///@brief   The InsideExclusive() variant of Inside() doesn’t consider points on the upper 
///       boundary to be inside the bounds. It is mostly useful with integer-typed bounds 
template <typename T> inline
bool InsideExclusive(const Point3<T> &p, const Bounds3<T> &b) {
    bool x = (p.x >= b.pMin.x) && (p.x < b.pMax.x);
    bool y = (p.y >= b.pMin.y) && (p.y < b.pMax.y);
    bool z = (p.z >= b.pMin.z) && (p.z < b.pMax.z);
    return (x && y && z); 
};

template <typename T, typename U> inline
Bounds3<T> Expand(const Bounds3<T> &b, U delta) {
    return Bounds3<T>(
        (b.pMin - Vector3<T>(delta, delta, delta)),
        (b.pMax + Vector3<T>(delta, delta, delta)));
};

template <typename T> inline
Vector3<T> Diagonal(const Bounds3<T> &b) {
    return b.pMax - b.pMin;
};

template <typename T> inline
T SurfaceArea(const Bounds3<T> &b) {
    Vector3<T> diagonal = Diagonal(b);
    return 2 * (diagonal.x * diagonal.y + diagonal.y * diagonal.z + diagonal.z * diagonal.x);
};

template <typename T> inline
T Volume(const Bounds3<T> &b) {
    Vector3<T> diagonal = Diagonal(b);
    return diagonal.x * diagonal.y * diagonal.z;
};

using Bounds2f = Bounds2<Float>;
using Bounds2i = Bounds2<int>;
using Bounds3f = Bounds3<Float>;
using Bounds3i = Bounds3<int>;

};
