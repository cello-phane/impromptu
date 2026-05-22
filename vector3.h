#ifndef VECTOR3_H
#define VECTOR3_H

#include "matrix4.h"

struct Vector3 {
    float x, y, z, w;
}; 

// We move Vector3 instances by value.

struct Vector3 Vector3_create_point(float x, float y, float z);
struct Vector3 Vector3_create_direction(float x, float y, float z);

// Operations.
struct Vector3 Vector3_add(struct Vector3 v1, struct Vector3 v2);
struct Vector3 Vector3_sub(struct Vector3 v1, struct Vector3 v2);
struct Vector3 Vector3_smul(struct Vector3 v, float s);
struct Vector3 Vector3_sdiv(struct Vector3 v, float s);
struct Vector3 Vector3_cross(struct Vector3 v1, struct Vector3 v2);
float          Vector3_dot(struct Vector3 v1, struct Vector3 v2);
float          Vector3_norm(struct Vector3 v);
float          Vector3_norm_squared(struct Vector3 v);
struct Vector3 Vector3_normalize(struct Vector3 v);

// Utility.
void           Vector3_print(const char *str, struct Vector3 v);

#endif