#pragma once
#include "shape.hpp"
#include "debug/debug.h"
#include <memory>

namespace eric {
class AreaLight;
class Material;
class MemoryArena;
class TransportMode;
class BSDF;
class BSSRDF;

class Primitive {
    public:
        BSDF *bsdf = nullptr;
        BSSRDF *bssrf = nullptr;
    public:
        virtual Bounds3f WorldBound() const = 0;
        virtual bool Intersect(const Ray& ray, SurfaceInteraction* isect) const = 0;
        virtual bool IntersectP(const Ray& ray) const = 0;
        virtual const AreaLight *GetAreaLight() const = 0;
        virtual const Material *GetMaterial() const = 0;
        virtual void ComputeScatteringFunction(SurfaceInteraction *isect, 
            MemoryArena &arena, TransportMode mode, bool allowMultipleLobes) const = 0;
};      

class GeometricPrimitive : public Primitive {
    private:
        std::shared_ptr<Shape> shape;
        std::shared_ptr<Material> material;
        std::shared_ptr<AreaLight> areaLight;
        MediumInterface mediumInterface;
    public:
        GeometricPrimitive(std::shared_ptr<Shape> shape, 
            std::shared_ptr<Material> material,
            std::shared_ptr<AreaLight> areaLight,
            MediumInterface mediumInterface): 
            shape(shape), material(material), areaLight(areaLight), mediumInterface(mediumInterface){};
    
        bool Intersect(const Ray& ray, SurfaceInteraction *isect) const override;
        
        bool IntersectP(const Ray& ray) const override;

        void GeometricPrimitive::ComputeScatteringFunctions(
            SurfaceInteraction *isect, MemoryArena &arena, TransportMode mode,
            bool allowMultipleLobes) const override;
};

class TransformedPrimitive : public Primitive {
    private:
        std::shared_ptr<Primitive> primitive;
        const AnimatedTransform PrimitiveToWorld;
    public:
        TransformedPrimitive(std::shared_ptr<Primitive> &primitive, 
            const AnimatedTransform &PrimitiveToWorld): 
                primitive(primitive), PrimitiveToWorld(PrimitiveToWorld) {};

        bool Intersect(const Ray &r, SurfaceInteraction *in) const;

        bool IntersectP(const Ray &r) const;

        const AreaLight *GetAreaLight() const { return nullptr; }
        
        const Material *GetMaterial() const { return nullptr; }

        void ComputeScatteringFunctions(SurfaceInteraction *isect, MemoryArena &arena,
            TransportMode mode, bool allowMultipleLobes) const {
            Severe("TransformedPrimitive::ComputeScatteringFunctions() shouldn't be called");
        }

        Bounds3f WorldBound() const { 
            return PrimitiveToWorld.MotionBounds(primitive->WorldBound());
        }
};

class Aggregate : public Primitive {
    private:
    
    public:
        Bounds3f WorldBound() const override;
        bool Intersect(const Ray& ray, SurfaceInteraction* isect) const override ;
        bool IntersectP(const Ray& ray) const override;

        // never should be called
        const AreaLight *GetAreaLight() const override {
            Severe("TransformedPrimitive::ComputeScatteringFunctions() shouldn't be called");
        };
        const Material *GetMaterial() const override {
            Severe("TransformedPrimitive::ComputeScatteringFunctions() shouldn't be called");
        };
        void ComputeScatteringFunction(SurfaceInteraction *isect, 
            MemoryArena &arena, TransportMode mode, bool allowMultipleLobes) const override {
                Severe("TransformedPrimitive::ComputeScatteringFunctions() shouldn't be called");
            };
}

}