#include <stdio.h>
#include <intrin.h> 
namespace eric {
#define Severe(msg) \
    {\
        fprintf(stderr, "%s\n", (msg));\
        __debugbreak();\
    }\

}