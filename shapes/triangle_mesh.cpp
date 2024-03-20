#include "triangle_mesh.hpp"

namespace eric {
Bounds3f Triangle::ObjectBound() const {
    Bounds3f b(mesh->P[*v], mesh->P[*(v + 1)]);
    return (*WorldToObject)(Union(b, mesh->P[*(v + 2)]));
};

Bounds3f Triangle::WorldBound() const {
    return Union(Bounds3f(mesh->P[*v], mesh->P[*(v + 1)]), mesh->P[*(v + 2)]);
}

bool Triangle::Intersect(
    const Ray &ray, Float *tHit, 
    SurfaceInteraction *isect, 
    bool testAlphaTexture = true) const {
        Point3f p1 = mesh->P[*v] - Vector3f(ray.o.x, ray.o.y, ray.o.z);
        Point3f p2 = mesh->P[*(v + 1)] - Vector3f(ray.o.x, ray.o.y, ray.o.z);
        Point3f p3 = mesh->P[*(v + 2)] - Vector3f(ray.o.x, ray.o.y, ray.o.z);
        Point2f uv1 = mesh->uv[*v];
        Point2f uv2 = mesh->uv[*(v + 1)];
        Point2f uv3 = mesh->uv[*(v + 2)];
        int kz = MaxComponentIndex(ray.d);
        int kx = kz + 1;
        if (kx > 2){
            kx == 0;
        }
        int ky = kx + 1;
        if (ky > 2) {
            ky == 0;
        }
        int x, y, z;
        x = ray.d[kx];
        y = ray.d[ky];
        z = ray.d[kz];

        p1.x = p1[kx] - x / z * p1[kz];
        p1.y = p1[ky] - y / z * p1[kz];
        p1.z = p1[kz] * 1 / z; 

        p2.x = p2[kx] - x / z * p2[kz];
        p2.y = p2[ky] - y / z * p2[kz];
        p2.z = p2[kz] * 1 / z; 

        p3.x = p3[kx] - x / z * p3[kz];
        p3.y = p3[ky] - y / z * p3[kz];
        p3.z = p3[kz] * 1 / z; 

        Vector2f p1p2(p2.x - p1.x, p2.y - p1.y);
        Vector2f p2p3(p3.x - p2.x, p3.y - p2.y);
        Vector2f p3p1(p1.x - p3.x, p1.y - p3.y);

        Float e0 = p1.x * p2.y - p1.y * p2.x;
        Float e1 = p3.x * p2.y - p3.y * p2.x;
        Float e2 = p1.x * p3.y - p1.y * p3.x;

        if (e0 < 0 || e1 < 0 || e2 < 0 && e0 > 0 || e1 > 0 || e2 > 0) {
            return false;
        }
        Float det = e0 + e1 + e2;
        if (det == 0) {
            return false;
        }
        
        Float tScaled = e0 * p1.z + e1 * p2.z + e2 * p3.z;
        if (det > 0 && (tScaled < 0 || tScaled > ray.tMax * det)) {
            return false;
        }
        if (det < 0 && (tScaled > 0 || tScaled < ray.tMax * det)) {
            return false;
        }
        
        Float invDet = 1 / det;
        Float b0 = e0 * invDet;
        Float b1 = e1 * invDet;
        Float b2 = e2 * invDet;
        Float t = tScaled * invDet;
        
        Point2f uvHit = b0 * uv1 + b1 * uv2 + b2 * uv3;
        Point3f pHit = ray.o + ray.d * t;

        Point3f _p1 = mesh->P[*v] - Vector3f(ray.o.x, ray.o.y, ray.o.z);
        Point3f _p2 = mesh->P[*(v + 1)] - Vector3f(ray.o.x, ray.o.y, ray.o.z);
        Point3f _p3 = mesh->P[*(v + 2)] - Vector3f(ray.o.x, ray.o.y, ray.o.z);

        Point2f delta_uv21 = Point2f(uv2.x - uv1.x, uv2.y - uv1.y);
        Point2f delta_uv32 = Point2f(uv3.x - uv2.x, uv3.y - uv2.y);

        Normal3f N = (Normal3f)Cross(Vector3f(_p3 - _p2), Vector3f(_p2 - _p1));
        /*
        *  dpdu * du1 + dpdv * dv1 = dp1;
        *  dpdu * du2 + dpdv * dv2 = dp2;
        */
        Vector3f dpdu, dpdv;
        if ((delta_uv21.y * delta_uv32.x - delta_uv32.y * delta_uv21.x) == 0.f) {
            dpdu = Vector3f(
                ((_p3.x - _p2.x) * delta_uv21.y - (_p2.x - _p1.x) * delta_uv32.y) / (delta_uv21.y * delta_uv32.x - delta_uv32.y * delta_uv21.x),
                ((_p3.y - _p2.y) * delta_uv21.y - (_p2.y - _p1.y) * delta_uv32.y) / (delta_uv21.y * delta_uv32.x - delta_uv32.y * delta_uv21.x),
                ((_p3.z - _p2.z) * delta_uv21.y - (_p2.z - _p1.z) * delta_uv32.y) / (delta_uv21.y * delta_uv32.x - delta_uv32.y * delta_uv21.x)
            );  

            dpdv = Vector3f(
                ((_p3.x - _p2.x) * delta_uv21.x - (_p2.x - _p1.x) * delta_uv32.x) / (delta_uv21.x * delta_uv32.y - delta_uv32.x * delta_uv21.y),
                ((_p3.x - _p2.x) * delta_uv21.x - (_p2.x - _p1.x) * delta_uv32.x) / (delta_uv21.x * delta_uv32.y - delta_uv32.x * delta_uv21.y),
                ((_p3.x - _p2.x) * delta_uv21.x - (_p2.x - _p1.x) * delta_uv32.x) / (delta_uv21.x * delta_uv32.y - delta_uv32.x * delta_uv21.y)
            );
        } else {
            CoordinateSystem(Vector3f(N), &dpdu, &dpdv);
        }

        if (testAlphaTexture && mesh->alphaMask) {
            SurfaceInteraction isectLocal =
                SurfaceInteraction(pHit, Vector3f(0, 0, 0), uvHit, -ray.d, 
                    dpdu, dpdv, Normal3f(0, 0, 0), Normal3f(0, 0, 0),
                    ray.time, this);
            if (mesh->alphaMask->Evaluate(isectLocal) == 0) {
                return false;
            }
        }

        *isect = SurfaceInteraction(pHit, Vector3f(0, 0, 0), uvHit, -ray.d, 
                    dpdu, dpdv, Normal3f(0, 0, 0), Normal3f(0, 0, 0),
                    ray.time, this);
        //Override surface normal in isect for triangle
        isect->shading.n = Normal3f(Cross(Vector3f(_p2 - _p1), Vector3f(_p3 - _p1)));
        //Ensure correct orientation of the geometric normal
        if (mesh->N) {
            isect->n = Faceforward(isect->n, Vector3f(isect->shading.n));
        } else if (reverseOrientation ^ transformSwapsHandedness) {
            isect->n = isect->shading.n = -isect->n;
        }
        //Initialize Triangle shading geometry
        if (mesh->N || mesh->S) {
            Normal3f ns;
            Vector3f ss, ts;
            if (mesh->N) {
                ns = Normal3f(
                    b0 * mesh->N[*v] + 
                    b1 * mesh->N[*(v + 1)] + 
                    b2 * mesh->N[*(v + 2)]
                );
                if (sqrt(Dot(Vector3f(ns), Vector3f(ns))) < 1e-7) {
                    ns = isect->n;
                }
            } else {
                ns = isect->n;
            }
            ns = Normalize(ns);

            if (mesh->S) {
                ss = Vector3f(
                    b0 * mesh->S[*v] +
                    b1 * mesh->S[*(v + 1)] +
                    b2 * mesh->S[*(v + 2)]
                );
                if (sqrt(Dot(ss, ss)) < 1e-7) {
                    ss = dpdu;
                }
            } else {
                ss = dpdu;
            }
            ss = Normalize(ss);

            ts = Cross(ns, ss);
            if (ts.length() > 0.f) {
                Normalize(ts);
                ss = Cross(ts, ns);
            } else {
                CoordinateSystem(Vector3f(ns), &ts, &ss);
            }

            Normal3f dndu, dndv;
            //TODO: normal partial derivative
        }
    };

bool Triangle::IntersectP(const Ray &ray,
    bool testAlphaTexture = true) const {
        Point3f p1 = mesh->P[*v] - Vector3f(ray.o.x, ray.o.y, ray.o.z);
        Point3f p2 = mesh->P[*(v + 1)] - Vector3f(ray.o.x, ray.o.y, ray.o.z);
        Point3f p3 = mesh->P[*(v + 2)] - Vector3f(ray.o.x, ray.o.y, ray.o.z);
        Point2f uv1 = mesh->uv[*v];
        Point2f uv2 = mesh->uv[*(v + 1)];
        Point2f uv3 = mesh->uv[*(v + 2)];
        int kz = MaxComponentIndex(ray.d);
        int kx = kz + 1;
        if (kx > 2){
            kx == 0;
        }
        int ky = kx + 1;
        if (ky > 2) {
            ky == 0;
        }
        int x, y, z;
        x = ray.d[kx];
        y = ray.d[ky];
        z = ray.d[kz];

        p1.x = p1[kx] - x / z * p1[kz];
        p1.y = p1[ky] - y / z * p1[kz];
        p1.z = p1[kz] * 1 / z; 

        p2.x = p2[kx] - x / z * p2[kz];
        p2.y = p2[ky] - y / z * p2[kz];
        p2.z = p2[kz] * 1 / z; 

        p3.x = p3[kx] - x / z * p3[kz];
        p3.y = p3[ky] - y / z * p3[kz];
        p3.z = p3[kz] * 1 / z; 

        Vector2f p1p2(p2.x - p1.x, p2.y - p1.y);
        Vector2f p2p3(p3.x - p2.x, p3.y - p2.y);
        Vector2f p3p1(p1.x - p3.x, p1.y - p3.y);

        Float e0 = p1.x * p2.y - p1.y * p2.x;
        Float e1 = p3.x * p2.y - p3.y * p2.x;
        Float e2 = p1.x * p3.y - p1.y * p3.x;

        if (e0 < 0 || e1 < 0 || e2 < 0 && e0 > 0 || e1 > 0 || e2 > 0) {
            return false;
        }
        Float det = e0 + e1 + e2;
        if (det == 0) {
            return false;
        }
        
        Float tScaled = e0 * p1.z + e1 * p2.z + e2 * p3.z;
        if (det > 0 && (tScaled < 0 || tScaled > ray.tMax * det)) {
            return false;
        }
        if (det < 0 && (tScaled > 0 || tScaled < ray.tMax * det)) {
            return false;
        }
        return true;
    };

Float Triangle::SurfaceArea() const {
    Point3f p1 = mesh->P[*v];
    Point3f p2 = mesh->P[*(v + 1)];
    Point3f p3 = mesh->P[*(v + 2)];
    Vector3f E1(p2 - p1), E2(p3 - p1);
    
    return 0.5f * Cross(E1.normalized(), E2.normalized()).length() * E1.length() * E2.length(); 
};

}