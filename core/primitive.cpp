#include "primitive.hpp"

namespace eric {
    bool GeometricPrimitive::Intersect(const Ray &ray, SurfaceInteraction *isect) const {
        Float tHit;
        if (!shape->Intersect(ray, &tHit, isect)){
            return false;
        }
        ray.tMax = tHit;
        // initialize mediumInterface
        return true;
    };

    bool GeometricPrimitive::IntersectP(const Ray &ray) const {
        return shape->IntersectP(ray);    
    }

    void GeometricPrimitive::ComputeScatteringFunctions(
        SurfaceInteraction *isect, MemoryArena &arena, TransportMode mode,
        bool allowMultipleLobes) const {
            if (material)
            material->ComputeScatteringFunctions(isect, arena, mode,
            allowMultipleLobes);
    }

    bool TransformedPrimitive::Intersect(const Ray &r, SurfaceInteraction *in) const {
        Transform p2w;
        PrimitiveToWorld.Interpolate(r.time, &p2w);
        Ray pRay = Inverse(p2w)(r);

        if (! primitive->Intersect(pRay, in)) {
            return false;
        }
        r.tMax = pRay.tMax; // update tmax for original ray
        if (!IsIdentity(p2w)) {
            *in = p2w(*in);
        }
    };

    bool TransformedPrimitive::IntersectP(const Ray &r) const {
        Transform p2w;
        PrimitiveToWorld.Interpolate(r.time, &p2w);
        Ray pRay = Inverse(p2w)(r);

        return primitive->IntersectP(pRay);
    };

    Bounds3f Aggregate::WorldBound() const {
        
    };

    bool Aggregate::Intersect(const Ray& ray, SurfaceInteraction* isect) const {

    };

    bool Aggregate::IntersectP(const Ray& ray) const {

    };
}