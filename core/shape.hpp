#pragma once
#include "transform.hpp"

namespace eric {
class Shape {
    public:
        const Transform *ObjectToWorld, *WorldToObject;
        const bool reverseOrientation;
        const bool transformSwapsHandedness;
    public:
        Shape::Shape(const Transform *ObjectToWorld,
            const Transform *WorldToObject, bool reverseOrientation)
            : ObjectToWorld(ObjectToWorld), WorldToObject(WorldToObject),
            reverseOrientation(reverseOrientation),
            transformSwapsHandedness(ObjectToWorld->SwapsHandedness()) {
        }

        virtual Bounds3f ObjectBound() const = 0; // return bounds in ojbect space
        
        virtual Bounds3f WorldBound() const {
            return (*ObjectToWorld)(ObjectBound());
        }

        virtual bool Intersect(const Ray &ray, Float *tHit, SurfaceInteraction *isect, bool testAlphaTexture = true) const = 0;

        virtual bool IntersectP(const Ray &ray,
            bool testAlphaTexture = true) const {
            Float tHit = ray.tMax;
            SurfaceInteraction isect;
            return Intersect(ray, &tHit, &isect, testAlphaTexture);
        }
        
        virtual Float SurfaceArea() const = 0;


};
};