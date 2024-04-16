#pragma once
#include <stdint.h>
#include <vector>

namespace eric {

constexpr int PageSize = 4 << 10;
constexpr int AlignSize = 64;

#define Align(ptr, alignment) \
    (uintptr_t)(((uintptr_t)(ptr) + alignment) & ~(alignment - 1))

#define AlignPtr(ptr) Align((ptr), (AlignSize))

struct MemoryPage {
    uint32_t refCount;
    
};

class MemoryArena {
    private:
        std::vector<uint8_t*> usedPages;
        std::vector<uint8_t*> freePages;
        uint8_t *current_page;
        uint32_t nUsed;
        uint32_t nFree;
        uint32_t offset;

    public:
        MemoryArena():nUsed(0), nFree(0), offset(0){};

        uint8_t* Alloc(uint32_t size) {
            uintptr_t end = (uintptr_t)current_page + PageSize;
            uintptr_t p = AlignPtr((uintptr_t)current_page + offset);
            
            if (end - p < size - 1) {
                usedPages.push_back()
            }
        };

};
}