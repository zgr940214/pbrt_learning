#include "shape.hpp"
#include <memory>

namespace eric {

enum class CurveType {
    FLAT,
    CYLINDER,
    RIBBON,
};

struct CurveCommon {
    const CurveType type;
    const Point3f cpOjb[4];
    const Float width[2];
    Normal3f n[2];
    Float normalAngle, invSinNormalAngle;
};

class Curve : public Shape {
    private:
        std::shared_ptr<CurveCommon> common;
        Float uMin, uMax;
    public:
        Curve(const Transform *object2World, const Transform *world2Ojbect, bool reverseOrientation,
            std::shared_ptr<CurveCommon> com, Float min, Float max): 
                Shape(object2World, world2Ojbect, reverseOrientation), 
                common(com), uMin(min), uMax(max) {};
        ~Curve() {};


};
}