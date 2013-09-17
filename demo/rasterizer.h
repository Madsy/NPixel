#ifndef RASTERIZER_H_GUARD
#define RASTERIZER_H_GUARD
#include <linealg.h>
#include "framebuffer.h"
//Tile size base. Must be POT
const int Q = 3;
//Actual Tile size
const int q = (1<<Q);

//Fixedpoint base for coefficients
const int coeff_precision_base = 15;
//Fixedpoint base for NDC coordinates (we convert screen-space coords to NDCs later)
const int ndc_precision_base = 26;
//Need 16-bits precision for z-buffer
const int depth_precision_base = 16;
//The actual scalars for the bases above (precision = 1 << base)
const float f_coeff_precision = (float)(1 << coeff_precision_base);
const float f_ndc_precision = (float)(1 << ndc_precision_base);
const float f_depth_precision = (float)(1 << depth_precision_base);
const int i_coeff_precision = 1 << coeff_precision_base;
const int i_ndc_precision = 1 << ndc_precision_base;
const int i_depth_precision = 1 << depth_precision_base;
//0.5 in coeff fixedpoint to subtract from u and v
const int i_uv_half = i_coeff_precision / 2;
//Base to use when converting from NDC coordinates to coefficients
const int base_diff = (ndc_precision_base - coeff_precision_base);
//Base to use when converting from NDC coordinates to depth
const int base_diff_z = (ndc_precision_base - depth_precision_base);


struct Tile {
    int FDX12, FDX23, FDX31;
    int FDY12, FDY23, FDY31;
    int CY1, CY2, CY3;
    int bw0, bw1, bw2, bw3; //w corner values
    int bz0, bz1, bz2, bz3; //w corner values
    int bu0, bu1, bu2, bu3; //w corner values
    int bv0, bv1, bv2, bv3; //w corner values
    int zMin, zMax; //for early z-culling
    bool operator<(const Tile& t) const {
        return zMin < t.zMin;
    }
};

struct TileSet {
    TileSet() : count(0) {
        tiles.resize(500);
    }
    int count;
    std::vector<Tile> tiles;
};

/* Convert from BGRA32 to the HW framebuffer pixelformat
   OBS! This throws away alpha! */
inline unsigned int ConvertToInternalColorFormat(unsigned int color)
{
    unsigned int b = (color & 0x000000FF);
    unsigned int g = (color & 0x0000FF00) >>  8;
    unsigned int r = (color & 0x00FF0000) >> 16;
    unsigned int a = (color & 0xFF000000) >> 24;

    /* Must shift down into new 1.0 range */
    b >>= wc_bLoss;
    g >>= wc_gLoss;
    r >>= wc_rLoss;
    a >>= wc_aLoss;

    return (r << wc_rShift) | (g << wc_gShift) | (b << wc_bShift);
}

//rasterizer_new.cpp
void DrawTriangles(unsigned int flags);
//rasterizer_deferred.cpp
void DrawTrianglesDeferred(unsigned int flags);

/* Computes coefficients for the vertex-scalars to interpolate.
   z, w, texture coordinates, normals, etc are multiplied with this matrix and later used
   in DrawTriangle */
Matrix3f ComputeCoeffMatrix(const VectorPOD4f& v1, const VectorPOD4f& v2, const VectorPOD4f& v3);

//extern int zmin;
//extern int zmax;

void SR_Render(unsigned int flags);
#endif

