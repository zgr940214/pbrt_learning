#pragma once
#include "shape.hpp"
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
    
    public:
    
}

}