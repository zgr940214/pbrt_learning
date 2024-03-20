#include "core/math.hpp"
#include "core/shape.hpp"
#include <memory>
#include <vector>

namespace eric {

template <typename T>
class Texture;

struct TriangleMesh {
    const Transform *w2o;
    const Transform *o2w;
    uint32_t nVertices;
    uint32_t nTriangles;
    std::vector<int> vertexIndices;
    std::unique_ptr<Point3f[]> P;
    std::unique_ptr<Vector3f[]> S; //tangent
    std::unique_ptr<Vector3f[]> N;
    std::unique_ptr<Point2f[]> uv;
    std::shared_ptr<Texture<Float>> alphaMask;

    TriangleMesh(const Transform *_w2o, const Transform *_o2w, uint32_t _nVertices, uint32_t _nTriangles, 
        const Point3f *_p, const Vector3f *_s, const Vector3f *_n, const Point2f *_uv, 
        const std::shared_ptr<Texture<Float>> _alphaMask):
        w2o(_w2o), o2w(_o2w), nVertices(_nVertices), nTriangles(_nTriangles) {
            P.reset(new Point3f[nVertices]);
            for (int i = 0; i < nVertices; i++) {
                P[i] = (*o2w)(_p[i]);
            }
            S.reset(new Vector3f[nTriangles]);
            N.reset(new Vector3f[nTriangles]);
            Transform tiO2W = Inverse(Transpose(*o2w));
            for (int i = 0; i < nTriangles; i++) {
                S[i] = (*o2w)(_s[i]);
                N[i] = (tiO2W)(_n[i]);
            }

        }
};

class Triangle : public Shape {
    private:
        const int *v;
        std::shared_ptr<TriangleMesh> mesh;
        uint32_t triangleIndex;
    public:
        Triangle(const Transform *ObjectToWorld, const Transform *WorldToObject, 
            bool reverseOrientation, const std::shared_ptr<TriangleMesh> &mesh, uint32_t triangleIndex):
            Shape(ObjectToWorld, WorldToObject, reverseOrientation), mesh(mesh), triangleIndex(triangleIndex)
            {
                v = &(mesh->vertexIndices[triangleIndex * 3]);
            };

        Bounds3f ObjectBound() const override;

        Bounds3f WorldBound() const override;

        bool Intersect(
            const Ray &ray, Float *tHit, 
            SurfaceInteraction *isect, 
            bool testAlphaTexture = true) const override;

        bool IntersectP(const Ray &ray,
            bool testAlphaTexture = true) const override;

        Float SurfaceArea() const;
};

std::vector<std::shared_ptr<Shape>> CreateTriangleMesh(
    const Transform *w2o, const Transform *o2w, bool reverseOrientation, 
    uint32_t nVertices, uint32_t nTriangles, const Point3f *p, 
    const Vector3f *s, const Vector3f *n, const Point2f *uv,
    const std::shared_ptr<Texture<Float>> alphaMask) {
        std::vector<std::shared_ptr<Shape>> triangles;
        std::shared_ptr<TriangleMesh> mesh = 
            std::make_shared<TriangleMesh>(w2o, o2w, nVertices, 
                nTriangles, p, s, n, uv, alphaMask);
        
        for (int i = 0; i < nTriangles; i++) {
            std::shared_ptr<Triangle> t = 
                std::make_shared<Triangle>(o2w, w2o, reverseOrientation, mesh, i);
            triangles.push_back(t);
        }
        return triangles;
    };




}