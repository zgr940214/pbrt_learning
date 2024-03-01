#include "triangle_mesh.hpp"

namespace eric {
Bounds3f Triangle::ObjectBound() const {
    Bounds3f b(mesh->P[*v], mesh->P[*(v + 1)]);
    return (*WorldToObject)(Union(b, mesh->P[*(v + 2)]));
};

bool Triangle::Intersect(
    const Ray &ray, Float *tHit, 
    SurfaceInteraction *isect, 
    bool testAlphaTexture = true) const {
        Point3f p1 = mesh->P[*v];
        Point3f p2 = mesh->P[*(v + 1)];
        Point3f p3 = mesh->P[*(v + 2)];
        Normal3f norm = (Normal3f)Cross(p2 - p1, p3 - p2);
        Float A = norm.x;
        Float B = norm.y;
        Float C = norm.z;
        EFloat dx(ray.d.x), dy(ray.d.y), dz(ray.d.z);
        EFloat ox(ray.o.x), oy(ray.o.y), oz(ray.o.z); 
        EFloat D = -1.f *(EFloat(A * p1.x * p1.x) + EFloat(B * p1.y * p1.y) + EFloat(C * p1.z * p1.z));

        EFloat a = dx * dx * EFloat(A) + dy * dy * EFloat(B) + dz * dz * EFloat(C);
        EFloat b = 2.f * (dx * ox + dy * oy + dz * oz);
        EFloat c = ox * ox * EFloat(A) + oy * oy * EFloat(B) + oz * oz * EFloat(C) + D;

        EFloat t0, t1;
        if (!Quadratic(a, b, c, &t0, &t1)) {
            return false;
        }

        if (t0.LowerBound() < 0.f || t0.UpperBound() > ray.tMax) {
            return false;
        }
        /*
        **  Möller–Trumbore u = (P - V1) dot (E2 cross E1) / N dot E1;
        **                  v = (P - V1) dot (E1 cross E2) / N dot E2;
        */
        Point3f p = ray(t0);
        Vector3f E2 = (p2 - p1), E1(p3 - p1);
        Float u0 = Dot(p - p1, Cross(E2, E1)) / Dot(Vector3f(norm), E1);
        Float v0 = Dot(p - p1, Cross(E1, E2)) / Dot(Vector3f(norm), E2);
        
        if (u0 > 1.f || u0 < 0.f || 
            v0 > 1.f || v0 < 0.f ||
            u0 + v0 > 1.f) {
                return false;
            }
        
        Vector3f dpdu = Vector3f(
            
        );
    };

bool Triangle::IntersectP(const Ray &ray,
    bool testAlphaTexture = true) const {

    };

Float Triangle::SurfaceArea() const {

};

}