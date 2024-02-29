#pragma once
#include "core/shape.hpp"

namespace eric {
class Sphere : public Shape {
    private:
        const Float radius;
        const Float zMin, zMax;
        const Float thetaMin, thetaMax;
        const Float phiMax;
    public:
        Sphere(const Transform *ObjectToWorld, const Transform *WorldToObject , bool reverseOrientation,
            Float radius, Float zMin, Float zMax, Float phiMax) : 
                Shape(ObjectToWorld, WorldToObject, reverseOrientation), radius(radius),
                zMin(Clamp(std::min(zMin, zMax), -radius, radius)),
                zMax(Clamp(std::max(zMin, zMax), -radius, radius)),  
                thetaMin(acos(Clamp(zMin / radius, -1, 1))), 
                thetaMax(acos(Clamp(zMax / radius, -1, 1))), 
                phiMax(Radians(Clamp(phiMax, 0, 360))) {};
        
        Bounds3f ObjectBound() const override;

        bool Intersect(
            const Ray &ray, Float *tHit, 
            SurfaceInteraction *isect, 
            bool testAlphaTexture = true) const override;

        bool IntersectP(const Ray &ray,
            bool testAlphaTexture = true) const override;

        Float SurfaceArea() const;
};

};