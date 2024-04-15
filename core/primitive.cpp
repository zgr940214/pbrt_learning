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
}