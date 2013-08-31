#ifndef RASTERIZER_H_GUARD
#define RASTERIZER_H_GUARD
#include <linealg.h>

void DrawTriangle();

/* Computes coefficients for the vertex-scalars to interpolate.
   z, w, texture coordinates, normals, etc are multiplied with this matrix and later used
   in DrawTriangle */
Matrix3f ComputeCoeffMatrix(const Vector4f& v1, const Vector4f& v2, const Vector4f& v3);

//extern int zmin;
//extern int zmax;

void SR_Render(unsigned int flags);
#endif

