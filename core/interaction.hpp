#pragma once
#include "geometry.hpp"
#include "shape.hpp"
#include "primitive.hpp"

namespace eric {

class MediumInterface;

struct Interaction {
    Point3f p;
    Vector3f pError;
    /*
            For interactions that lie along a ray (either from a ray–shape intersection or from a ray
        passing through participating media), the negative ray direction is stored in wo, which
        corresponds to ωo, the notation we use for the outgoing direction when computing
        lighting at points. For other types of interaction points where the notion of an outgoing
        direction doesn’t apply (e.g., those found by randomly sampling points on the surface of
        shapes), wo has the value (0, 0, 0).
    */
    Vector3f wo;
    Normal3f n;
    MediumInterface mediumInterface;
    Float time;

    Interaction(){};
    Interaction(const Point3f &p, const Normal3f &n, const Vector3f &pError, 
        const Vector3f &wo, Float time, const MediumInterface &mediumInterface ):
        p(p), n(n), pError(pError), wo(wo), time(time), mediumInterface(mediumInterface){};

    bool IsSurfaceInteraction() const{
        return n != (Normal3f());
    }

};

class Shape;
class SurfaceInteraction : public Interaction {
    public:
        Point2f uv;
        Vector3f dpdu, dpdv;
        Normal3f dndu, dndv;
        const Shape *shape = nullptr;
        const Primitive * primitive = nullptr;
        struct {
            Normal3f n;
            Vector3f dpdu, dpdv;
            Normal3f dndu, dndv;
        } shading;

    public:
        SurfaceInteraction() {};
        SurfaceInteraction(const Point3f &p,
            const Vector3f &pError, const Point2f &uv, const Vector3f &wo,
            const Vector3f &dpdu, const Vector3f &dpdv,
            const Normal3f &dndu, const Normal3f &dndv,
            Float time, const Shape *shape)
            : Interaction(p, Normal3f(Normalize(Cross(dpdu, dpdv))), pError, wo,
            time, nullptr),
            uv(uv), dpdu(dpdu), dpdv(dpdv), dndu(dndu), dndv(dndv),
            shape(shape) {
                //Initialize shading geometry from true geometry
                shading.n = n;
                shading.dpdu = dpdu;
                shading.dpdv = dpdv;
                shading.dndu = dndu;
                shading.dndv = dndv;
                //Adjust normal based on orientation and handedness
                if (shape && (shape->reverseOrientation ^
                    shape->transformSwapsHandedness)) {
                    n *= -1;
                    shading.n *= -1;
                }
            };
        void SurfaceInteraction::SetShadingGeometry(const Vector3f &dpdus,
            const Vector3f &dpdvs, const Normal3f &dndus,
            const Normal3f &dndvs, bool orientationIsAuthoritative) {
            //Compute shading.n for SurfaceInteraction
                shading.n = Normalize((Normal3f)Cross(dpdus, dpdvs));
                if (shape && (shape->reverseOrientation ^
                    shape->transformSwapsHandedness))
                    shading.n = -shading.n;
                if (orientationIsAuthoritative)
                    n = Faceforward(n, shading.n);
                else
                    shading.n = Faceforward(shading.n, n);
                //Initialize shading partial derivative values
                shading.dpdu = dpdus;
                shading.dpdv = dpdvs;
                shading.dndu = dndus;
                shading.dndv = dndvs;
            };
    
};
}