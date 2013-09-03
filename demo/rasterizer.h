#ifndef RASTERIZER_H_GUARD
#define RASTERIZER_H_GUARD
#include <linealg.h>

void DrawTriangles(unsigned int flags);

/* Computes coefficients for the vertex-scalars to interpolate.
   z, w, texture coordinates, normals, etc are multiplied with this matrix and later used
   in DrawTriangle */
Matrix3f ComputeCoeffMatrix(const VectorPOD4f& v1, const VectorPOD4f& v2, const VectorPOD4f& v3);

//extern int zmin;
//extern int zmax;

void SR_Render(unsigned int flags);
#endif

