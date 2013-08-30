#ifndef CLIPPLANE_H_GUARD
#define CLIPPLANE_H_GUARD
#include <linealg.h>
void clip_triangle(Vector4f plane, unsigned int flags);
int classifyTriangle(const Vector4f& v1, const Vector4f& v2, const Vector4f& v3);
#endif
