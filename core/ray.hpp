#pragma once
#include "geometry.hpp"

namespace eric {
/*
-----------------------------------------------------------------------------------------
** Ray: represent by a point along with a direction
** RayDifferential: contains additional two auxiliary rays which offsets the origin by one sample in the x and y
----------------------------------------------------------------------------------------- 
*/  
constexpr Float infinity =  std::numeric_limits<Float>::infinity();

class Medium;
class Ray {
    public:
        Point3f o;
        Vector3f d;
        const Medium* medium;
        mutable Float tMax;
        Float time;
    public:
        Ray():tMax(infinity), time(0.f), medium(nullptr){};
        explicit Ray(const Point3f &o, const Vector3f &d, 
            Float t = infinity, Float time = 0.f, const Medium *m = nullptr):
            o(o), d(d), tMax(t), time(t), medium(m){};

        Point3f operator()(Float t) const;
};

class RayDifferential : public Ray {
    public:
        bool hasDifferentials;
        Point3f rxOrigin, ryOrigin;
        Vector3f rxDirection, ryDirection;
        
    public:
        RayDifferential(){hasDifferentials = false;};
        RayDifferential(const Point3f &o, const Vector3f &d, 
            Float t = infinity, Float time = 0.f,const Medium *m = nullptr):
            Ray(o, d, t, time, m){hasDifferentials = false;};
        RayDifferential(const Ray &ray):Ray(ray){hasDifferentials = false;};

        void ScaleDifferential(Float s) {
            rxOrigin = o + (rxOrigin - o) * s;
            ryOrigin = o + (ryOrigin - o) * s;
            rxDirection = d + (rxDirection - d) * s;
            ryDirection = d + (ryDirection - d) * s;
        }
};
};