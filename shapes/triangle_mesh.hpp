#include "core/math.hpp"
#include "core/transform.hpp"

namespace eric {
struct TriangleMesh {
    const Transform *w2o;
    const Transform *o2w;
    uint32_t nVertices;
    uint32_t nTriangles;
    uint32_t *vertexIndices;
    Point3f *P;
    Vector3f *S;
    Vector3f *N;
    Point2f *uv;
    bool alphaMask;
}
}