#include "util.h"

Uint32 rand32(Uint64 *rng)
{
    *rng = *rng*0x3243f6a8885a308d + 1;
    return *rng >> 32;
}