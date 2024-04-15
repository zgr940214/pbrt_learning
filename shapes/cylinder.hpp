#pragma once
#include "core/math.hpp"
#include "core/shape.hpp"

namespace eric {
class Cylinder : public Shape {
    private:
        const Float radius, zMax, zMin, phiMax;
    public:
        Cylinder(const Transform* ObjectToWorld, const Transform* WorldToObject, 
            bool reverseOrientation, Float radius, Float zMin, Float zMax, Float phiMax)
            : Shape(ObjectToWorld, WorldToObject, reverseOrientation), radius(radius),
            zMin(zMin), zMax(zMax), phiMax(Radians(Clamp(phiMax, 0.f, 360.f))) {};
        
        Bounds3f ObjectBound() const override;

        bool Intersect(const Ray &ray, Float* tHit,  SurfaceInteraction *isect, bool testAlphaTexture = true) const override;

        bool IntersectP(const Ray &ray, bool testAlphaTexture = true) const override;

        Float SurfaceArea() const override;
}; 
};