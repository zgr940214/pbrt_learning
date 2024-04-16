#pragma once
#include "core/bounds.hpp"
#include "core/primitive.hpp"
#include <vector>

namespace eric {

class BVHAccel : Aggregate{
    public:
        enum class SplitMethod {SAH, HLBVH, Middle, EqualCounts};
        struct BVHPrimitiveInfo {
            BVHPrimitiveInfo(size_t primitiveNumber, 
                Bounds3f bound): primitiveNumber(primitiveNumber), bounds(bounds),
                    centroid(0.5f * bound.pMax + 0.5f * bound.pMin) {};

            size_t primitiveNumber;
            Bounds3f bounds;
            Point3f centroid;
        };

        struct BVHBuildNode {
            void InitLeaf(int first, int n, const Bounds3f b) {
                firstPrimOffset = first;
                nPrimitives = n;
                bounds = b;
                children[0] = children[1] = nullptr;
            };
            
            void InitInterior(int first, int n, BVHBuildNode *c0, BVHBuildNode *c1){
                firstPrimOffset = first;
                nPrimitives = n;
                children[0] = c0;
                children[1] = c1;
                bounds = Union(c0->bounds, c1->bounds);
            };

            Bounds3f bounds;
            BVHBuildNode *children[2];
            int splitAxis, firstPrimOffset, nPrimitives;
        };

    private:
        const int MaxPrimsInNode;
        const SplitMethod splitMethod;
        std::vector<std::shared_ptr<Primitive>> primitives;
    public:
        BVHAccel(const std::vector<std::shared_ptr<Primitive>> prims, 
            const int MaxPrims, const SplitMethod method):
                MaxPrimsInNode(std::min(255, MaxPrims)), splitMethod(method), primitives(prims) {
                    std::vector<BVHPrimitiveInfo> primitiveInfo;
                    for (int i = 0; i < prims.size(); i++) {

                    }
                    // construct bvh tree from primitives
                };

        BVHBuildNode* HLBVHBuild();

        BVHBuildNode* RecursiveBuild();
        

};

}