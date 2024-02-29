#include "core/math.hpp"
#include "core/shape.hpp"

namespace eric {
class Disk : public Shape {
    private:
        const Float height, radius, innerRadius, phiMax;
    public:
        Disk(const Transform *ObjectToWorld, const Transform *WorldToObject, bool reverseOrientation,
            Float height, Float raidus, Float innerRadius, Float phiMax): 
            Shape(ObjectToWorld, WorldToObject, reverseOrientation), height(height),
            radius(radius), innerRadius(innerRadius), phiMax(Radians(Clamp(phiMax, 0.f, 360.f))) {};
        
        Bounds3f ObjectBound() const override;

        bool Intersect (const Ray &ray, Float *tHit, SurfaceInteraction *isect, bool testAlphaTexture = true) const override;

        bool IntersectP (const Ray &ray, bool testAlphaTexture = true) const override;

        Float SurfaceArea() const override;
};

}