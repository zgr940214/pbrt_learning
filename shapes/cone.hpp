#include "core/math.hpp"
#include "core/shape.hpp"

namespace eric {
class Cone : public Shape {
    private:
        Float radius, height; 
        Float phiMax; 
    public:
        Cone(const Transform *ObjectToWorld, const Transform *WorldToObject, 
            bool reverseOrientation, Float raidus, Float height, Float phiMax) : Shape(ObjectToWorld, WorldToObject, reverseOrientation),
            radius(radius), height(height), phiMax(Radians(Clamp(phiMax, 0.f, 360.f))) {};
        
        Bounds3f ObjectBound() const override;

        bool Intersect(const Ray &ray, Float *tHit, SurfaceInteraction *isect, bool testAlphaTexture = true) const override;

        bool IntersectP(const Ray &ray, bool testAlphaTexture = true) const override;

        Float SurfaceArea() const override;
};
}