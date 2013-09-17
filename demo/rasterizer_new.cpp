#include <utility>
#include <algorithm>
#include <vector>
#include <cstdio>
#include <SDL/SDL.h>
#include <linealg.h>
#include <fixedpoint.h>
#include "rasterizer.h"
#include "clipplane.h"
#include "vertexdata.h"
#include "framebuffer.h"
#include "texture.h"
#include "myassert.h"
#if 0
void DrawTriangles(unsigned int flags)
{
    using std::min;
    using std::max;

    /* pointer for 15, 16, 24 and 32-bit bpp */
    void* colorbuffer = wc_colorbuffer->Lock();
    unsigned short* depthbuffer = (unsigned short*)wc_depthbuffer->Ptr();

    const unsigned int* tbuf = &wc_texture0->texels[0];
    const unsigned int width = wc_colorbuffer->w;
    const unsigned int height = wc_colorbuffer->h;
    const unsigned int iTw = wc_texture0->width;
    const unsigned int iTh = wc_texture0->height;

    for(int i=0; i<wc_vertices->size(); i+=3) {
        VectorPOD4f& v1 = (*wc_vertices)[i+0];
        VectorPOD4f& v2 = (*wc_vertices)[i+2];
        VectorPOD4f& v3 = (*wc_vertices)[i+1];

        VectorPOD4f& tc1 = (*wc_tcoords0)[i+0];
        VectorPOD4f& tc2 = (*wc_tcoords0)[i+1];
        VectorPOD4f& tc3 = (*wc_tcoords0)[i+2];

        // 28.4 fixed-point coordinates
        const int Y1 = (int)(16.0f * v1.y);
        const int Y2 = (int)(16.0f * v2.y);
        const int Y3 = (int)(16.0f * v3.y);
        const int X1 = (int)(16.0f * v1.x);
        const int X2 = (int)(16.0f * v2.x);
        const int X3 = (int)(16.0f * v3.x);

        // Deltas
        const int DX12 = X1 - X2;
        const int DX23 = X2 - X3;
        const int DX31 = X3 - X1;
        const int DY12 = Y1 - Y2;
        const int DY23 = Y2 - Y3;
        const int DY31 = Y3 - Y1;

        // Fixed-point deltas
        const int FDX12 = DX12 << 4;
        const int FDX23 = DX23 << 4;
        const int FDX31 = DX31 << 4;
        const int FDY12 = DY12 << 4;
        const int FDY23 = DY23 << 4;
        const int FDY31 = DY31 << 4;

        // Bounding rectangle
        int minx = (min(X1, min(X2, X3)) + 0xF) >> 4;
        int maxx = (max(X1, max(X2, X3)) + 0xF) >> 4;
        int miny = (min(Y1, min(Y2, Y3)) + 0xF) >> 4;
        int maxy = (max(Y1, max(Y2, Y3)) + 0xF) >> 4;

        // Start in corner of a 8x8 block alligned to block size
        minx &= ~(q - 1);
        miny &= ~(q - 1);
        maxx = (maxx + (q - 1)) & ~(q - 1);
        maxy = (maxy + (q - 1)) & ~(q - 1);

        // Half-edge constants
        int C1 = DY12 * X1 - DX12 * Y1;
        int C2 = DY23 * X2 - DX23 * Y2;
        int C3 = DY31 * X3 - DX31 * Y3;

        // Correct for fill convention
        if(DY12 < 0 || (DY12 == 0 && DX12 > 0)) C1++;
        if(DY23 < 0 || (DY23 == 0 && DX23 > 0)) C2++;
        if(DY31 < 0 || (DY31 == 0 && DX31 > 0)) C3++;

        // Loop through blocks
        for(int y = miny; y < maxy; y += q) {
            for(int x = minx; x < maxx; x += q) {

                // Corners of block
                int x0 = x;
                int x1 = (x + q - 1);
                int y0 = y;
                int y1 = (y + q - 1);

                // test block against x and y frustum planes
                bool px0min = x0 >= 0;
                bool px0max = x0 < width;
                bool py0min = y0 >= 0;
                bool py0max = y0 < height;

                bool px1min = x1 >= 0;
                bool px1max = x1 < width;
                bool py1min = y1 >= 0;
                bool py1max = y1 < height;

                int pflags0 = (px0min << 3) | (px1min << 2) | (px0max << 1) | px1max;
                int pflags1 = (py0min << 3) | (py1min << 2) | (py0max << 1) | py1max;

                int pflags = (pflags0<<4) | pflags1;

                if(pflags != 0xFF)
                    continue;

                x0 <<= 4;
                x1 <<= 4;
                y0 <<= 4;
                y1 <<= 4;

                // Evaluate half-space functions
                bool a00 = C1 + DX12 * y0 - DY12 * x0 > 0;
                bool a10 = C1 + DX12 * y0 - DY12 * x1 > 0;
                bool a01 = C1 + DX12 * y1 - DY12 * x0 > 0;
                bool a11 = C1 + DX12 * y1 - DY12 * x1 > 0;

                int a = (a00 << 0) | (a10 << 1) | (a01 << 2) | (a11 << 3);

                bool b00 = C2 + DX23 * y0 - DY23 * x0 > 0;
                bool b10 = C2 + DX23 * y0 - DY23 * x1 > 0;
                bool b01 = C2 + DX23 * y1 - DY23 * x0 > 0;
                bool b11 = C2 + DX23 * y1 - DY23 * x1 > 0;

                int b = (b00 << 0) | (b10 << 1) | (b01 << 2) | (b11 << 3);

                bool c00 = C3 + DX31 * y0 - DY31 * x0 > 0;
                bool c10 = C3 + DX31 * y0 - DY31 * x1 > 0;
                bool c01 = C3 + DX31 * y1 - DY31 * x0 > 0;
                bool c11 = C3 + DX31 * y1 - DY31 * x1 > 0;

                int c = (c00 << 0) | (c10 << 1) | (c01 << 2) | (c11 << 3);

                // Skip block when outside an edge
                if(a == 0x0 || b == 0x0 || c == 0x0) continue;

                // Coefficients for the equation s/w = Ax + By + C
                // Where x and y are in NDC space
                const int Az = v1.z * f_depth_precision;
                const int Bz = v3.z * f_depth_precision;
                const int Cz = v2.z * f_depth_precision;
                const int Aw = v1.w  * f_coeff_precision;
                const int Bw = v3.w  * f_coeff_precision;
                const int Cw = v2.w  * f_coeff_precision;
                const int Au = tc1.x * f_coeff_precision;
                const int Bu = tc2.x * f_coeff_precision;
                const int Cu = tc3.x * f_coeff_precision;
                const int Av = tc1.y * f_coeff_precision;
                const int Bv = tc2.y * f_coeff_precision;
                const int Cv = tc3.y * f_coeff_precision;

                const float fHalfWidthInv = 2.0f / (float)width;
                const float fHalfHeightInv = 2.0f / (float)height;

                int NDC_x_step = fHalfWidthInv * f_ndc_precision; //1 subtracted later
                int NDC_y_step = fHalfHeightInv * f_ndc_precision; //1 subtracted later
                int NDC_x0 = x * NDC_x_step;  //min x
                int NDC_y0 = y * NDC_y_step;  //min y
                int NDC_x1 = (x + q - 1) * NDC_x_step; //max x
                int NDC_y1 = (y + q - 1) * NDC_y_step; //max y

                int bwx0 = (NDC_x0 - i_ndc_precision) >> base_diff;
                int bwx1 = (NDC_x1 - i_ndc_precision) >> base_diff;
                int bwy0 = (NDC_y0 - i_ndc_precision) >> base_diff;
                int bwy1 = (NDC_y1 - i_ndc_precision) >> base_diff;

                int bzx0 = (NDC_x0 - i_ndc_precision) >> base_diff_z;
                int bzx1 = (NDC_x1 - i_ndc_precision) >> base_diff_z;
                int bzy0 = (NDC_y0 - i_ndc_precision) >> base_diff_z;
                int bzy1 = (NDC_y1 - i_ndc_precision) >> base_diff_z;

                int bw0 = (((long long)Aw*bwx0 + (long long)Bw*bwy0)>>coeff_precision_base) + Cw; //top left
                int bw1 = (((long long)Aw*bwx0 + (long long)Bw*bwy1)>>coeff_precision_base) + Cw; //bottom left
                int bw2 = (((long long)Aw*bwx1 + (long long)Bw*bwy0)>>coeff_precision_base) + Cw; //top right
                int bw3 = (((long long)Aw*bwx1 + (long long)Bw*bwy1)>>coeff_precision_base) + Cw; //bottom right

                int bz0 = (((long long)Az*bzx0 + (long long)Bz*bzy0)>>depth_precision_base) + Cz; //top left
                int bz1 = (((long long)Az*bzx0 + (long long)Bz*bzy1)>>depth_precision_base) + Cz; //bottom left
                int bz2 = (((long long)Az*bzx1 + (long long)Bz*bzy0)>>depth_precision_base) + Cz; //top right
                int bz3 = (((long long)Az*bzx1 + (long long)Bz*bzy1)>>depth_precision_base) + Cz; //bottom right

                //Compute u for the corners of the tile
                int bu0 = (((long long)Au*bwx0 + (long long)Bu*bwy0)>>coeff_precision_base) + Cu; //top left
                int bu1 = (((long long)Au*bwx0 + (long long)Bu*bwy1)>>coeff_precision_base) + Cu; //bottom left
                int bu2 = (((long long)Au*bwx1 + (long long)Bu*bwy0)>>coeff_precision_base) + Cu; //top right
                int bu3 = (((long long)Au*bwx1 + (long long)Bu*bwy1)>>coeff_precision_base) + Cu; //bottom right

                //Compute v for the corners of the tile
                int bv0 = (((long long)Av*bwx0 + (long long)Bv*bwy0)>>coeff_precision_base) + Cv; //top left
                int bv1 = (((long long)Av*bwx0 + (long long)Bv*bwy1)>>coeff_precision_base) + Cv; //bottom left
                int bv2 = (((long long)Av*bwx1 + (long long)Bv*bwy0)>>coeff_precision_base) + Cv; //top right
                int bv3 = (((long long)Av*bwx1 + (long long)Bv*bwy1)>>coeff_precision_base) + Cv; //bottom right


                // Delta *and* slope. Since we know that a tile width
                // and height equals a constant power of two (which is constant q), we can use
                // it for slope as-is without dividing by deltaY or deltaX
                // bwSlope0 = ((bw1 - bw0) << Q) >> Q; is pointless.
                // Instead, pretend that bwSlope is delta and slope, and is
                // in Q fixedpoint. Accumulate it and shift down with Q later
                int bwSlopeY0 = bw1 - bw0; //vertical slope left
                int bwSlopeY1 = bw3 - bw2; //vertical slope right
                int bzSlopeY0 = bz1 - bz0; //vertical slope left
                int bzSlopeY1 = bz3 - bz2; //vertical slope right
                int buSlopeY0 = bu1 - bu0; //vertical slope left
                int buSlopeY1 = bu3 - bu2; //vertical slope right
                int bvSlopeY0 = bv1 - bv0; //vertical slope left
                int bvSlopeY1 = bv3 - bv2; //vertical slope right

                //Since bwSlopeY0 is delta and slope, we're really using
                //a fixedpoint base of Q, so correct the start value
                //bwSlopeYAccum represents the vertical interpolated w values
                int bwSlopeYAccum0 = bw0 << Q; //start vertical w left
                int bwSlopeYAccum1 = bw2 << Q; //start vertical w right
                int bzSlopeYAccum0 = bz0 << Q; //start vertical z left
                int bzSlopeYAccum1 = bz2 << Q; //start vertical z right
                int buSlopeYAccum0 = bu0 << Q; //start vertical u left
                int buSlopeYAccum1 = bu2 << Q; //start vertical u right
                int bvSlopeYAccum0 = bv0 << Q; //start vertical v left
                int bvSlopeYAccum1 = bv2 << Q; //start vertical v right

                // Accept whole block when totally covered
                if(a == 0xF && b == 0xF && c == 0xF) {
                    unsigned int col = y*width;
                    for(int iy = y; iy < y + q; iy++) {
                        //Slope and delta for horisontal interpolation for w
                        //Same rule applies here as for bwSlopeYAccum0 and bwSlopeYAccum1
                        int bwSlopeX0 = bwSlopeYAccum1 - bwSlopeYAccum0;
                        int bzSlopeX0 = bzSlopeYAccum1 - bzSlopeYAccum0;
                        int buSlopeX0 = buSlopeYAccum1 - buSlopeYAccum0;
                        int bvSlopeX0 = bvSlopeYAccum1 - bvSlopeYAccum0;
                        //Since we pretend that delta and slope is the same
                        //bwSlopeXAccum0 is in Q+Q fixedpoint, so correct start value
                        int bwSlopeXAccum0 = bwSlopeYAccum0 << Q; //start horisontal w left
                        int bzSlopeXAccum0 = bzSlopeYAccum0 << Q; //start horisontal z left
                        int buSlopeXAccum0 = buSlopeYAccum0 << Q; //start horisontal u left
                        int bvSlopeXAccum0 = bvSlopeYAccum0 << Q; //start horisontal v left

                        for(int ix = x; ix < x + q; ix++) {
                            unsigned short z = bzSlopeXAccum0 >> (Q*2);
                            if(z < depthbuffer[ix + col]) {
                                depthbuffer[ix + col] = z;
                                int uw = buSlopeXAccum0>>(Q*2);
                                int vw = bvSlopeXAccum0>>(Q*2);
                                int w = bwSlopeXAccum0>>(Q*2);
                                //If we interpolate w instead of 1/w ..
                                if(w) w = ((long long)1<<(coeff_precision_base * 2)) / w;
                                unsigned int u = ((long long)uw*w*(iTw - 1)) >> (coeff_precision_base);
                                unsigned int v = ((long long)vw*w*(iTh - 1)) >> (coeff_precision_base);
                                u >>= coeff_precision_base;
                                v >>= coeff_precision_base;
                                u = clamp(u, 0u, iTw-1u);
                                v = clamp(v, 0u, iTh-1u);
                                //colorbuffer[ix + col] = tbuf[u + v*iTw];
                                switch(wc_bpp) {
                                case 15:
                                case 16: {
                                    unsigned int color32 = tbuf[u + v*iTw];
                                    unsigned short color16 = ConvertToInternalColorFormat(color32);
                                    ((unsigned short*)colorbuffer)[ix + col] = color16;
                                    break;
                                }
                                case 24: //just ignores alpha. 1-byte padded anyway
                                case 32: {
                                    unsigned int color32 = tbuf[u + v*iTw];
                                    ((unsigned int*)colorbuffer)[ix + col] = color32;
                                    break;
                                }
                                }
                            }
                            bwSlopeXAccum0 += bwSlopeX0;
                            bzSlopeXAccum0 += bzSlopeX0;
                            buSlopeXAccum0 += buSlopeX0;
                            bvSlopeXAccum0 += bvSlopeX0;
                        }
                        //vertically interpolated values
                        //bSlopeXAccum0 is the result of horisontally
                        //interpolating between these two points
                        bwSlopeYAccum0 += bwSlopeY0;
                        bwSlopeYAccum1 += bwSlopeY1;
                        bzSlopeYAccum0 += bzSlopeY0;
                        bzSlopeYAccum1 += bzSlopeY1;
                        buSlopeYAccum0 += buSlopeY0;
                        buSlopeYAccum1 += buSlopeY1;
                        bvSlopeYAccum0 += bvSlopeY0;
                        bvSlopeYAccum1 += bvSlopeY1;
                        col += width;
                    }
                } else { // Partially covered
                    int CY1 = C1 + DX12 * y0 - DY12 * x0;
                    int CY2 = C2 + DX23 * y0 - DY23 * x0;
                    int CY3 = C3 + DX31 * y0 - DY31 * x0;
                    unsigned int col = y*width;
                    for(int iy = y; iy < y + q; iy++) {
                        int CX1 = CY1;
                        int CX2 = CY2;
                        int CX3 = CY3;
                        int bwSlopeX0 = bwSlopeYAccum1 - bwSlopeYAccum0;
                        int bzSlopeX0 = bzSlopeYAccum1 - bzSlopeYAccum0;
                        int buSlopeX0 = buSlopeYAccum1 - buSlopeYAccum0;
                        int bvSlopeX0 = bvSlopeYAccum1 - bvSlopeYAccum0;
                        int bwSlopeXAccum0 = bwSlopeYAccum0 << Q;
                        int bzSlopeXAccum0 = bzSlopeYAccum0 << Q;
                        int buSlopeXAccum0 = buSlopeYAccum0 << Q;
                        int bvSlopeXAccum0 = bvSlopeYAccum0 << Q;

                        for(int ix = x; ix < x + q; ix++) {
                            if(CX1 > 0 && CX2 > 0 && CX3 > 0) {
                                unsigned short z = bzSlopeXAccum0 >> (Q*2);
                                if(z < depthbuffer[ix + col]) {
                                    depthbuffer[ix + col] = z;
                                    int uw = buSlopeXAccum0>>(Q*2);
                                    int vw = bvSlopeXAccum0>>(Q*2);
                                    int w = bwSlopeXAccum0>>(Q*2);
                                    //If we interpolate w instead of 1/w ..
                                    if(w) w = ((long long)1<<(coeff_precision_base * 2)) / w;
                                    unsigned int u = ((long long)uw*w*(iTw - 1)) >> (coeff_precision_base);
                                    unsigned int v = ((long long)vw*w*(iTh - 1)) >> (coeff_precision_base);
                                    u >>= coeff_precision_base;
                                    v >>= coeff_precision_base;
                                    u = clamp(u, 0u, iTw-1u);
                                    v = clamp(v, 0u, iTh-1u);
                                    //colorbuffer[ix + col] = tbuf[u + v*iTw];
                                    switch(wc_bpp) {
                                    case 15:
                                    case 16: {
                                        unsigned int color32 = tbuf[u + v*iTw];
                                        unsigned short color16 = ConvertToInternalColorFormat(color32);
                                        ((unsigned short*)colorbuffer)[ix + col] = color16;
                                        break;
                                    }
                                    case 24: //just ignores alpha. 1-byte padded anyway
                                    case 32: {
                                        unsigned int color32 = tbuf[u + v*iTw];
                                        ((unsigned int*)colorbuffer)[ix + col] = color32;
                                        break;
                                    }
                                    }
                                }
                            }
                            bwSlopeXAccum0 += bwSlopeX0;
                            bzSlopeXAccum0 += bzSlopeX0;
                            buSlopeXAccum0 += buSlopeX0;
                            bvSlopeXAccum0 += bvSlopeX0;
                            CX1 -= FDY12;
                            CX2 -= FDY23;
                            CX3 -= FDY31;
                        }
                        bwSlopeYAccum0 += bwSlopeY0;
                        bwSlopeYAccum1 += bwSlopeY1;
                        bzSlopeYAccum0 += bzSlopeY0;
                        bzSlopeYAccum1 += bzSlopeY1;
                        buSlopeYAccum0 += buSlopeY0;
                        buSlopeYAccum1 += buSlopeY1;
                        bvSlopeYAccum0 += bvSlopeY0;
                        bvSlopeYAccum1 += bvSlopeY1;
                        CY1 += FDX12;
                        CY2 += FDX23;
                        CY3 += FDX31;
                        col += width;
                    }
                }
            }
        }
    }
    wc_colorbuffer->Unlock();
}
#endif

struct TileInfo {
    int bwSlopeY0;
    int bwSlopeY1;
    int bzSlopeY0;
    int bzSlopeY1;
    int buSlopeY0;
    int buSlopeY1;
    int bvSlopeY0;
    int bvSlopeY1;
    int bwSlopeYAccum0;
    int bwSlopeYAccum1;
    int bzSlopeYAccum0;
    int bzSlopeYAccum1;
    int buSlopeYAccum0;
    int buSlopeYAccum1;
    int bvSlopeYAccum0;
    int bvSlopeYAccum1;
    /* used for partial coverage */
    int CY1; // = C1 + DX12 * y0 - DY12 * x0;
    int CY2; // = C2 + DX23 * y0 - DY23 * x0;
    int CY3; // = C3 + DX31 * y0 - DY31 * x0;
    int FDY12;
    int FDY23;
    int FDY31;
    int FDX12;
    int FDX23;
    int FDX31;
    /* misc */
    unsigned int x;
    unsigned int y;
};


/*
    void* buf_color;
    void* buf_depth;
    void* buf_texture;
    unsigned int texWidth;
    unsigned int texHeight;
*/
template<bool filled>
static inline void DrawTileFilled888(TileInfo info)
{
    /* pointer for 15, 16, 24 and 32-bit bpp */
    unsigned int* buf_color = (unsigned int*)wc_colorbuffer->Lock();
    unsigned short* buf_depth = (unsigned short*)wc_depthbuffer->Ptr();
    const unsigned int* buf_texture = &wc_texture0->texels[0];
    const unsigned short texWidth = wc_texture0->width;
    const unsigned short texHeight = wc_texture0->height;

    int col = info.y * wc_width;

    for(int iy = info.y; iy < info.y + q; iy++) {
        int bwSlopeX0 = info.bwSlopeYAccum1 - info.bwSlopeYAccum0;
        int bzSlopeX0 = info.bzSlopeYAccum1 - info.bzSlopeYAccum0;
        int buSlopeX0 = info.buSlopeYAccum1 - info.buSlopeYAccum0;
        int bvSlopeX0 = info.bvSlopeYAccum1 - info.bvSlopeYAccum0;
        int bwSlopeXAccum0 = info.bwSlopeYAccum0 << Q;
        int bzSlopeXAccum0 = info.bzSlopeYAccum0 << Q;
        int buSlopeXAccum0 = info.buSlopeYAccum0 << Q;
        int bvSlopeXAccum0 = info.bvSlopeYAccum0 << Q;
        int CX1_0, CX2_0, CX3_0;
        if(!filled) {
            CX1_0 = info.CY1;
            CX2_0 = info.CY2;
            CX3_0 = info.CY3;
        }

        // Compiler hint. x and y are always aligned to tile size
        info.x &= ~(q - 1);
        info.y &= ~(q - 1);
        
        for(int ix = info.x; ix < info.x + q; ix += 2) {
            unsigned short z0 = bzSlopeXAccum0 >> (Q*2);
            unsigned short z1 = (bzSlopeXAccum0 + bzSlopeX0) >> (Q*2);
            unsigned char z0Cmp = z0 < buf_depth[ix + col];
            unsigned char z1Cmp = z1 < buf_depth[1 + ix + col];
            if(!filled){
                int CX1_1 = CX1_0 - info.FDY12;
                int CX2_1 = CX2_0 - info.FDY23;
                int CX3_1 = CX3_0 - info.FDY31;
                // Both depth-test and interior test must pass
                z0Cmp &= (CX1_0 > 0 && CX2_0 > 0 && CX3_0 > 0);
                z1Cmp &= (CX1_1 > 0 && CX2_1 > 0 && CX3_1 > 0);
            }
            unsigned char z2Flag = (z1Cmp << 1) | z0Cmp;
            // Two pixels, four cases:
            // 0 : Both failed z-test or interior test. Only update accumulators.
            // 1 : Pixel 0 passed tests. Single write
            // 2 : Pixel 1 passed test. Single write
            // 3 : Both pixels passed. Write both with a single write.
            switch(z2Flag) {
            case 0: {     
                bwSlopeXAccum0 += (bwSlopeX0 << 1);
                bzSlopeXAccum0 += (bzSlopeX0 << 1);
                buSlopeXAccum0 += (buSlopeX0 << 1);
                bvSlopeXAccum0 += (bvSlopeX0 << 1);
                // Skip two pixels ahead
                if(!filled){
                    CX1_0 -= (info.FDY12 << 1);
                    CX2_0 -= (info.FDY23 << 1);
                    CX3_0 -= (info.FDY31 << 1);
                }
            }
            case 1: {
                int uw = buSlopeXAccum0>>(Q*2);
                int vw = bvSlopeXAccum0>>(Q*2);
                int w = bwSlopeXAccum0>>(Q*2);
                //If we interpolate w instead of 1/w ..
                if(w) w = ((long long)1<<(coeff_precision_base * 2)) / w;
                unsigned short u = ((long long)uw*w*(texWidth - 1)) >> (coeff_precision_base * 2);
                unsigned short v = ((long long)vw*w*(texHeight - 1)) >> (coeff_precision_base * 2);
                u = clamp((unsigned int)u, 0u, texWidth-1u);
                v = clamp((unsigned int)v, 0u, texHeight-1u);
                unsigned int color32 = buf_texture[u + v*texWidth];
                buf_depth[ix + col] = z0;
                buf_color[ix + col] = color32;

                bwSlopeXAccum0 += (bwSlopeX0 << 1);
                bzSlopeXAccum0 += (bzSlopeX0 << 1);
                buSlopeXAccum0 += (buSlopeX0 << 1);
                bvSlopeXAccum0 += (bvSlopeX0 << 1);
                // Skip two pixels ahead
                if(!filled){
                    CX1_0 -= (info.FDY12 << 1);
                    CX2_0 -= (info.FDY23 << 1);
                    CX3_0 -= (info.FDY31 << 1);
                }
                break;
            }
            case 2: {
                int uw = (buSlopeXAccum0 + buSlopeX0)>>(Q*2);
                int vw = (bvSlopeXAccum0 + bvSlopeX0)>>(Q*2);
                int w  = (bwSlopeXAccum0 + bwSlopeX0)>>(Q*2);
                //If we interpolate w instead of 1/w ..
                if(w) w = ((long long)1<<(coeff_precision_base * 2)) / w;
                unsigned short u = ((long long)uw*w*(texWidth - 1)) >> (coeff_precision_base);
                unsigned short v = ((long long)vw*w*(texHeight - 1)) >> (coeff_precision_base);
                u >>= coeff_precision_base;
                v >>= coeff_precision_base;
                u = clamp((unsigned int)u, 0u, texWidth-1u);
                v = clamp((unsigned int)v, 0u, texHeight-1u);
                unsigned int color32 = buf_texture[u + v*texWidth];
                buf_depth[1 + ix + col] = z1;
                buf_color[1 + ix + col] = color32;

                bwSlopeXAccum0 += (bwSlopeX0 << 1);
                bzSlopeXAccum0 += (bzSlopeX0 << 1);
                buSlopeXAccum0 += (buSlopeX0 << 1);
                bvSlopeXAccum0 += (bvSlopeX0 << 1);
                // Skip two pixels ahead
                if(!filled){
                    CX1_0 -= (info.FDY12 << 1);
                    CX2_0 -= (info.FDY23 << 1);
                    CX3_0 -= (info.FDY31 << 1);
                }
                break;
            }
            case 3: {
                int uw0 =  buSlopeXAccum0>>(Q*2);
                int uw1 = (buSlopeXAccum0 + buSlopeX0)>>(Q*2);
                int vw0 =  bvSlopeXAccum0>>(Q*2);
                int vw1 = (bvSlopeXAccum0 + bvSlopeX0)>>(Q*2);
                int w0  =  bwSlopeXAccum0>>(Q*2);
                int w1  = (bwSlopeXAccum0 + bwSlopeX0)>>(Q*2);
                //If we interpolate w instead of 1/w ..
                if(w0) w0 = ((long long)1<<(coeff_precision_base * 2)) / w0;
                if(w1) w1 = ((long long)1<<(coeff_precision_base * 2)) / w1;
                unsigned short u0 = ((long long)uw0*w0*(texWidth  - 1)) >> (coeff_precision_base * 2);
                unsigned short u1 = ((long long)uw1*w1*(texWidth  - 1)) >> (coeff_precision_base * 2);
                unsigned short v0 = ((long long)vw0*w0*(texHeight - 1)) >> (coeff_precision_base * 2);
                unsigned short v1 = ((long long)vw1*w1*(texHeight - 1)) >> (coeff_precision_base * 2);
                u0 = clamp((unsigned int)u0, 0u, texWidth-1u);
                u1 = clamp((unsigned int)u1, 0u, texWidth-1u);
                v0 = clamp((unsigned int)v0, 0u, texHeight-1u);
                v1 = clamp((unsigned int)v1, 0u, texHeight-1u);
                unsigned int color32_0 = buf_texture[u0 + v0*texWidth];
                unsigned int color32_1 = buf_texture[u1 + v1*texWidth];
                ((unsigned int*)buf_depth)[(ix + col)>>1] = ((unsigned int)z1 << 16) | z0;
                ((unsigned long long*)buf_color)[(ix + col)>>1] = ((unsigned long long)color32_1 << 32) | color32_0;

                bwSlopeXAccum0 += (bwSlopeX0 << 1);
                bzSlopeXAccum0 += (bzSlopeX0 << 1);
                buSlopeXAccum0 += (buSlopeX0 << 1);
                bvSlopeXAccum0 += (bvSlopeX0 << 1);
                // Skip two pixels ahead
                if(!filled){
                    CX1_0 -= (info.FDY12 << 1);
                    CX2_0 -= (info.FDY23 << 1);
                    CX3_0 -= (info.FDY31 << 1);
                }
                break;
            }
            } //end switch
        } //end ix
        info.bwSlopeYAccum0 += info.bwSlopeY0;
        info.bwSlopeYAccum1 += info.bwSlopeY1;
        info.bzSlopeYAccum0 += info.bzSlopeY0;
        info.bzSlopeYAccum1 += info.bzSlopeY1;
        info.buSlopeYAccum0 += info.buSlopeY0;
        info.buSlopeYAccum1 += info.buSlopeY1;
        info.bvSlopeYAccum0 += info.bvSlopeY0;
        info.bvSlopeYAccum1 += info.bvSlopeY1;
        col += wc_width;
        if(!filled) {
            info.CY1 += info.FDX12;
            info.CY2 += info.FDX23;
            info.CY3 += info.FDX31;
        }
    }
}

template<bool filled>
static inline void DrawTileFilled565(TileInfo info)
{
    /* pointer for 15, 16, 24 and 32-bit bpp */
    unsigned short* buf_color = (unsigned short*)wc_colorbuffer->Lock();
    unsigned short* buf_depth = (unsigned short*)wc_depthbuffer->Ptr();
    const unsigned int* buf_texture = &wc_texture0->texels[0];
    const unsigned short texWidth = wc_texture0->width;
    const unsigned short texHeight = wc_texture0->height;

    int col = info.y * wc_width;

    for(int iy = info.y; iy < info.y + q; iy++) {
        int bwSlopeX0 = info.bwSlopeYAccum1 - info.bwSlopeYAccum0;
        int bzSlopeX0 = info.bzSlopeYAccum1 - info.bzSlopeYAccum0;
        int buSlopeX0 = info.buSlopeYAccum1 - info.buSlopeYAccum0;
        int bvSlopeX0 = info.bvSlopeYAccum1 - info.bvSlopeYAccum0;
        int bwSlopeXAccum0 = info.bwSlopeYAccum0 << Q;
        int bzSlopeXAccum0 = info.bzSlopeYAccum0 << Q;
        int buSlopeXAccum0 = info.buSlopeYAccum0 << Q;
        int bvSlopeXAccum0 = info.bvSlopeYAccum0 << Q;
        int CX1_0, CX2_0, CX3_0;
        if(!filled) {
            CX1_0 = info.CY1;
            CX2_0 = info.CY2;
            CX3_0 = info.CY3;
        }

        // Compiler hint. x and y are always aligned to tile size
        //info.x &= ~(q - 1);
        //info.y &= ~(q - 1);
        
        for(int ix = info.x; ix < info.x + q; ix += 2) {
            unsigned short z0 = bzSlopeXAccum0 >> (Q*2);
            unsigned short z1 = (bzSlopeXAccum0 + bzSlopeX0) >> (Q*2);
            unsigned char z0Cmp = z0 < buf_depth[ix + col];
            unsigned char z1Cmp = z1 < buf_depth[1 + ix + col];
            if(!filled){
                int CX1_1 = CX1_0 - info.FDY12;
                int CX2_1 = CX2_0 - info.FDY23;
                int CX3_1 = CX3_0 - info.FDY31;
                // Both depth-test and interior test must pass
                z0Cmp &= (unsigned char)(CX1_0 > 0 && CX2_0 > 0 && CX3_0 > 0);
                z1Cmp &= (unsigned char)(CX1_1 > 0 && CX2_1 > 0 && CX3_1 > 0);
            }
            unsigned char z2Flag = (z1Cmp << 1) | z0Cmp;
            // Two pixels, four cases:
            // 0 : Both failed z-test or interior test. Only update accumulators.
            // 1 : Pixel 0 passed tests. Single write
            // 2 : Pixel 1 passed test. Single write
            // 3 : Both pixels passed. Write both with a single write.
            switch(z2Flag) {
            case 0: {     
                bwSlopeXAccum0 += (bwSlopeX0 << 1);
                bzSlopeXAccum0 += (bzSlopeX0 << 1);
                buSlopeXAccum0 += (buSlopeX0 << 1);
                bvSlopeXAccum0 += (bvSlopeX0 << 1);
                // Skip two pixels ahead
                if(!filled){
                    CX1_0 -= (info.FDY12 * 2);
                    CX2_0 -= (info.FDY23 * 2);
                    CX3_0 -= (info.FDY31 * 2);
                }
            }
            case 1: {
                /*
                int uw = buSlopeXAccum0>>(Q*2);
                int vw = bvSlopeXAccum0>>(Q*2);
                int w = bwSlopeXAccum0>>(Q*2);
                //If we interpolate w instead of 1/w ..
                if(w) w = ((long long)1<<(coeff_precision_base * 2)) / w;
                unsigned short u = ((long long)uw*w*(texWidth - 1)) >> (coeff_precision_base * 2);
                unsigned short v = ((long long)vw*w*(texHeight - 1)) >> (coeff_precision_base * 2);
                u = clamp((unsigned int)u, 0u, texWidth-1u);
                v = clamp((unsigned int)v, 0u, texHeight-1u);
                unsigned int color32 = buf_texture[u + v*texWidth];
                unsigned short color16 = ConvertToInternalColorFormat(color32);
                buf_depth[ix + col] = z0;
                buf_color[ix + col] = color16;
                */
                buf_depth[ix + col] = z0;
                buf_color[ix + col] = ConvertToInternalColorFormat(0x00FF0000);
                
                bwSlopeXAccum0 += (bwSlopeX0 << 1);
                bzSlopeXAccum0 += (bzSlopeX0 << 1);
                buSlopeXAccum0 += (buSlopeX0 << 1);
                bvSlopeXAccum0 += (bvSlopeX0 << 1);
                // Skip two pixels ahead
                if(!filled){
                    CX1_0 -= (info.FDY12 * 2);
                    CX2_0 -= (info.FDY23 * 2);
                    CX3_0 -= (info.FDY31 * 2);
                }
                break;
            }
            case 2: {
                /*
                int uw = (buSlopeXAccum0 + buSlopeX0)>>(Q*2);
                int vw = (bvSlopeXAccum0 + bvSlopeX0)>>(Q*2);
                int w  = (bwSlopeXAccum0 + bwSlopeX0)>>(Q*2);
                //If we interpolate w instead of 1/w ..
                if(w) w = ((long long)1<<(coeff_precision_base * 2)) / w;
                unsigned short u = ((long long)uw*w*(texWidth - 1)) >> (coeff_precision_base * 2);
                unsigned short v = ((long long)vw*w*(texHeight - 1)) >> (coeff_precision_base * 2);
                u = clamp((unsigned int)u, 0u, texWidth-1u);
                v = clamp((unsigned int)v, 0u, texHeight-1u);
                unsigned int color32 = buf_texture[u + v*texWidth];
                unsigned short color16 = ConvertToInternalColorFormat(color32);
                buf_depth[1 + ix + col] = z1;
                buf_color[1 + ix + col] = color16;
                */
                buf_depth[1+ ix + col] = z1;
                buf_color[1+ ix + col] = ConvertToInternalColorFormat(0x0000FF00);
                
                bwSlopeXAccum0 += (bwSlopeX0 << 1);
                bzSlopeXAccum0 += (bzSlopeX0 << 1);
                buSlopeXAccum0 += (buSlopeX0 << 1);
                bvSlopeXAccum0 += (bvSlopeX0 << 1);
                // Skip two pixels ahead
                if(!filled){
                    CX1_0 -= (info.FDY12 * 2);
                    CX2_0 -= (info.FDY23 * 2);
                    CX3_0 -= (info.FDY31 * 2);
                }
                break;
            }
            case 3: {
                /*
                int uw0 =  buSlopeXAccum0>>(Q*2);
                int uw1 = (buSlopeXAccum0 + buSlopeX0)>>(Q*2);
                int vw0 =  bvSlopeXAccum0>>(Q*2);
                int vw1 = (bvSlopeXAccum0 + bvSlopeX0)>>(Q*2);
                int w0  =  bwSlopeXAccum0>>(Q*2);
                int w1  = (bwSlopeXAccum0 + bwSlopeX0)>>(Q*2);
                //If we interpolate w instead of 1/w ..
                if(w0) w0 = ((long long)1<<(coeff_precision_base * 2)) / w0;
                if(w1) w1 = ((long long)1<<(coeff_precision_base * 2)) / w1;
                unsigned short u0 = ((long long)uw0*w0*(texWidth  - 1)) >> (coeff_precision_base * 2);
                unsigned short u1 = ((long long)uw1*w1*(texWidth  - 1)) >> (coeff_precision_base * 2);
                unsigned short v0 = ((long long)vw0*w0*(texHeight - 1)) >> (coeff_precision_base * 2);
                unsigned short v1 = ((long long)vw1*w1*(texHeight - 1)) >> (coeff_precision_base * 2);
                u0 = clamp((unsigned int)u0, 0u, texWidth-1u);
                u1 = clamp((unsigned int)u1, 0u, texWidth-1u);
                v0 = clamp((unsigned int)v0, 0u, texHeight-1u);
                v1 = clamp((unsigned int)v1, 0u, texHeight-1u);
                unsigned int color32_0 = buf_texture[u0 + v0*texWidth];
                unsigned int color32_1 = buf_texture[u1 + v1*texWidth];
                unsigned short color16_0 = ConvertToInternalColorFormat(color32_0);
                unsigned short color16_1 = ConvertToInternalColorFormat(color32_1);
                ((unsigned int*)buf_depth)[(ix + col)>>1] = ((unsigned int)z1 << 16) | z0;
                //((unsigned int*)buf_color)[(ix + col)>>1] = ((unsigned int)color16_1 << 16) | color16_0;
                buf_color[    ix + col] = color16_0;
                buf_color[1 + ix + col] = color16_1;
                */
                buf_depth[   ix + col] = z0;
                buf_depth[1+ ix + col] = z1;
                buf_color[   ix + col] = ConvertToInternalColorFormat(0x000000FF);
                buf_color[1+ ix + col] = ConvertToInternalColorFormat(0x000000FF);
                
                bwSlopeXAccum0 += (bwSlopeX0 << 1);
                bzSlopeXAccum0 += (bzSlopeX0 << 1);
                buSlopeXAccum0 += (buSlopeX0 << 1);
                bvSlopeXAccum0 += (bvSlopeX0 << 1);
                // Skip two pixels ahead
                if(!filled){
                    CX1_0 -= (info.FDY12 * 2);
                    CX2_0 -= (info.FDY23 * 2);
                    CX3_0 -= (info.FDY31 * 2);
                }
                break;
            }
            } //end switch
        } //end ix
        info.bwSlopeYAccum0 += info.bwSlopeY0;
        info.bwSlopeYAccum1 += info.bwSlopeY1;
        info.bzSlopeYAccum0 += info.bzSlopeY0;
        info.bzSlopeYAccum1 += info.bzSlopeY1;
        info.buSlopeYAccum0 += info.buSlopeY0;
        info.buSlopeYAccum1 += info.buSlopeY1;
        info.bvSlopeYAccum0 += info.bvSlopeY0;
        info.bvSlopeYAccum1 += info.bvSlopeY1;
        col += wc_width;
        if(!filled) {
            info.CY1 += info.FDX12;
            info.CY2 += info.FDX23;
            info.CY3 += info.FDX31;
        }
    }
}

//unsigned short color16 = ConvertToInternalColorFormat(color32);

void DrawTriangles(unsigned int flags)
{
    using std::min;
    using std::max;

    /* pointer for 15, 16, 24 and 32-bit bpp */
    void* colorbuffer = wc_colorbuffer->Lock();
    unsigned short* depthbuffer = (unsigned short*)wc_depthbuffer->Ptr();

    const unsigned int* tbuf = &wc_texture0->texels[0];
    const unsigned int width = wc_colorbuffer->w;
    const unsigned int height = wc_colorbuffer->h;
    const unsigned int iTw = wc_texture0->width;
    const unsigned int iTh = wc_texture0->height;

    for(int i=0; i<wc_vertices->size(); i+=3) {
        VectorPOD4f& v1 = (*wc_vertices)[i+0];
        VectorPOD4f& v2 = (*wc_vertices)[i+2];
        VectorPOD4f& v3 = (*wc_vertices)[i+1];

        VectorPOD4f& tc1 = (*wc_tcoords0)[i+0];
        VectorPOD4f& tc2 = (*wc_tcoords0)[i+1];
        VectorPOD4f& tc3 = (*wc_tcoords0)[i+2];

        // 28.4 fixed-point coordinates
        const int Y1 = (int)(16.0f * v1.y);
        const int Y2 = (int)(16.0f * v2.y);
        const int Y3 = (int)(16.0f * v3.y);
        const int X1 = (int)(16.0f * v1.x);
        const int X2 = (int)(16.0f * v2.x);
        const int X3 = (int)(16.0f * v3.x);

        // Deltas
        const int DX12 = X1 - X2;
        const int DX23 = X2 - X3;
        const int DX31 = X3 - X1;
        const int DY12 = Y1 - Y2;
        const int DY23 = Y2 - Y3;
        const int DY31 = Y3 - Y1;

        // Fixed-point deltas
        const int FDX12 = DX12 << 4;
        const int FDX23 = DX23 << 4;
        const int FDX31 = DX31 << 4;
        const int FDY12 = DY12 << 4;
        const int FDY23 = DY23 << 4;
        const int FDY31 = DY31 << 4;

        // Bounding rectangle
        int minx = (min(X1, min(X2, X3)) + 0xF) >> 4;
        int maxx = (max(X1, max(X2, X3)) + 0xF) >> 4;
        int miny = (min(Y1, min(Y2, Y3)) + 0xF) >> 4;
        int maxy = (max(Y1, max(Y2, Y3)) + 0xF) >> 4;

        // Start in corner of a 8x8 block alligned to block size
        minx &= ~(q - 1);
        miny &= ~(q - 1);
        maxx = (maxx + (q - 1)) & ~(q - 1);
        maxy = (maxy + (q - 1)) & ~(q - 1);

        // Half-edge constants
        int C1 = DY12 * X1 - DX12 * Y1;
        int C2 = DY23 * X2 - DX23 * Y2;
        int C3 = DY31 * X3 - DX31 * Y3;

        // Correct for fill convention
        if(DY12 < 0 || (DY12 == 0 && DX12 > 0)) C1++;
        if(DY23 < 0 || (DY23 == 0 && DX23 > 0)) C2++;
        if(DY31 < 0 || (DY31 == 0 && DX31 > 0)) C3++;

        // Loop through blocks
        for(int y = miny; y < maxy; y += q) {
            for(int x = minx; x < maxx; x += q) {
                // Corners of block
                int x0 = x;
                int x1 = (x + q - 1);
                int y0 = y;
                int y1 = (y + q - 1);
                // test block against x and y frustum planes
                bool px0min = x0 >= 0;
                bool px0max = x0 < width;
                bool py0min = y0 >= 0;
                bool py0max = y0 < height;

                bool px1min = x1 >= 0;
                bool px1max = x1 < width;
                bool py1min = y1 >= 0;
                bool py1max = y1 < height;

                int pflags0 = (px0min << 3) | (px1min << 2) | (px0max << 1) | px1max;
                int pflags1 = (py0min << 3) | (py1min << 2) | (py0max << 1) | py1max;

                int pflags = (pflags0<<4) | pflags1;

                if(pflags != 0xFF)
                    continue;

                x0 <<= 4;
                x1 <<= 4;
                y0 <<= 4;
                y1 <<= 4;

                // Evaluate half-space functions
                bool a00 = C1 + DX12 * y0 - DY12 * x0 > 0;
                bool a10 = C1 + DX12 * y0 - DY12 * x1 > 0;
                bool a01 = C1 + DX12 * y1 - DY12 * x0 > 0;
                bool a11 = C1 + DX12 * y1 - DY12 * x1 > 0;

                int a = (a00 << 0) | (a10 << 1) | (a01 << 2) | (a11 << 3);

                bool b00 = C2 + DX23 * y0 - DY23 * x0 > 0;
                bool b10 = C2 + DX23 * y0 - DY23 * x1 > 0;
                bool b01 = C2 + DX23 * y1 - DY23 * x0 > 0;
                bool b11 = C2 + DX23 * y1 - DY23 * x1 > 0;

                int b = (b00 << 0) | (b10 << 1) | (b01 << 2) | (b11 << 3);

                bool c00 = C3 + DX31 * y0 - DY31 * x0 > 0;
                bool c10 = C3 + DX31 * y0 - DY31 * x1 > 0;
                bool c01 = C3 + DX31 * y1 - DY31 * x0 > 0;
                bool c11 = C3 + DX31 * y1 - DY31 * x1 > 0;

                int c = (c00 << 0) | (c10 << 1) | (c01 << 2) | (c11 << 3);

                // Skip block when outside an edge
                if(a == 0x0 || b == 0x0 || c == 0x0) continue;

                // Coefficients for the equation s/w = Ax + By + C
                // Where x and y are in NDC space
                const int Az = v1.z * f_depth_precision;
                const int Bz = v3.z * f_depth_precision;
                const int Cz = v2.z * f_depth_precision;
                const int Aw = v1.w  * f_coeff_precision;
                const int Bw = v3.w  * f_coeff_precision;
                const int Cw = v2.w  * f_coeff_precision;
                const int Au = tc1.x * f_coeff_precision;
                const int Bu = tc2.x * f_coeff_precision;
                const int Cu = tc3.x * f_coeff_precision;
                const int Av = tc1.y * f_coeff_precision;
                const int Bv = tc2.y * f_coeff_precision;
                const int Cv = tc3.y * f_coeff_precision;

                const float fHalfWidthInv = 2.0f / (float)width;
                const float fHalfHeightInv = 2.0f / (float)height;

                int NDC_x_step = fHalfWidthInv * f_ndc_precision; //1 subtracted later
                int NDC_y_step = fHalfHeightInv * f_ndc_precision; //1 subtracted later
                int NDC_x0 = x * NDC_x_step;  //min x
                int NDC_y0 = y * NDC_y_step;  //min y
                int NDC_x1 = (x + q - 1) * NDC_x_step; //max x
                int NDC_y1 = (y + q - 1) * NDC_y_step; //max y

                int bwx0 = (NDC_x0 - i_ndc_precision) >> base_diff;
                int bwx1 = (NDC_x1 - i_ndc_precision) >> base_diff;
                int bwy0 = (NDC_y0 - i_ndc_precision) >> base_diff;
                int bwy1 = (NDC_y1 - i_ndc_precision) >> base_diff;

                int bzx0 = (NDC_x0 - i_ndc_precision) >> base_diff_z;
                int bzx1 = (NDC_x1 - i_ndc_precision) >> base_diff_z;
                int bzy0 = (NDC_y0 - i_ndc_precision) >> base_diff_z;
                int bzy1 = (NDC_y1 - i_ndc_precision) >> base_diff_z;

                int bw0 = (((long long)Aw*bwx0 + (long long)Bw*bwy0)>>coeff_precision_base) + Cw; //top left
                int bw1 = (((long long)Aw*bwx0 + (long long)Bw*bwy1)>>coeff_precision_base) + Cw; //bottom left
                int bw2 = (((long long)Aw*bwx1 + (long long)Bw*bwy0)>>coeff_precision_base) + Cw; //top right
                int bw3 = (((long long)Aw*bwx1 + (long long)Bw*bwy1)>>coeff_precision_base) + Cw; //bottom right

                int bz0 = (((long long)Az*bzx0 + (long long)Bz*bzy0)>>depth_precision_base) + Cz; //top left
                int bz1 = (((long long)Az*bzx0 + (long long)Bz*bzy1)>>depth_precision_base) + Cz; //bottom left
                int bz2 = (((long long)Az*bzx1 + (long long)Bz*bzy0)>>depth_precision_base) + Cz; //top right
                int bz3 = (((long long)Az*bzx1 + (long long)Bz*bzy1)>>depth_precision_base) + Cz; //bottom right

                //Compute u for the corners of the tile
                int bu0 = (((long long)Au*bwx0 + (long long)Bu*bwy0)>>coeff_precision_base) + Cu; //top left
                int bu1 = (((long long)Au*bwx0 + (long long)Bu*bwy1)>>coeff_precision_base) + Cu; //bottom left
                int bu2 = (((long long)Au*bwx1 + (long long)Bu*bwy0)>>coeff_precision_base) + Cu; //top right
                int bu3 = (((long long)Au*bwx1 + (long long)Bu*bwy1)>>coeff_precision_base) + Cu; //bottom right

                //Compute v for the corners of the tile
                int bv0 = (((long long)Av*bwx0 + (long long)Bv*bwy0)>>coeff_precision_base) + Cv; //top left
                int bv1 = (((long long)Av*bwx0 + (long long)Bv*bwy1)>>coeff_precision_base) + Cv; //bottom left
                int bv2 = (((long long)Av*bwx1 + (long long)Bv*bwy0)>>coeff_precision_base) + Cv; //top right
                int bv3 = (((long long)Av*bwx1 + (long long)Bv*bwy1)>>coeff_precision_base) + Cv; //bottom right

                TileInfo info;
                // Delta *and* slope. Since we know that a tile width
                // and height equals a constant power of two (which is constant q), we can use
                // it for slope as-is without dividing by deltaY or deltaX
                // bwSlope0 = ((bw1 - bw0) << Q) >> Q; is pointless.
                // Instead, pretend that bwSlope is delta and slope, and is
                // in Q fixedpoint. Accumulate it and shift down with Q later
                info.bwSlopeY0 = bw1 - bw0; //vertical slope left
                info.bwSlopeY1 = bw3 - bw2; //vertical slope right
                info.bzSlopeY0 = bz1 - bz0; //vertical slope left
                info.bzSlopeY1 = bz3 - bz2; //vertical slope right
                info.buSlopeY0 = bu1 - bu0; //vertical slope left
                info.buSlopeY1 = bu3 - bu2; //vertical slope right
                info.bvSlopeY0 = bv1 - bv0; //vertical slope left
                info.bvSlopeY1 = bv3 - bv2; //vertical slope right

                //Since bwSlopeY0 is delta and slope, we're really using
                //a fixedpoint base of Q, so correct the start value
                //bwSlopeYAccum represents the vertical interpolated w values
                info.bwSlopeYAccum0 = bw0 << Q; //start vertical w left
                info.bwSlopeYAccum1 = bw2 << Q; //start vertical w right
                info.bzSlopeYAccum0 = bz0 << Q; //start vertical z left
                info.bzSlopeYAccum1 = bz2 << Q; //start vertical z right
                info.buSlopeYAccum0 = bu0 << Q; //start vertical u left
                info.buSlopeYAccum1 = bu2 << Q; //start vertical u right
                info.bvSlopeYAccum0 = bv0 << Q; //start vertical v left
                info.bvSlopeYAccum1 = bv2 << Q; //start vertical v right

                /* Slopes for coverage */
                info.FDY12 = FDY12;
                info.FDY23 = FDY23;
                info.FDY31 = FDY31;
                info.FDX12 = FDX12;
                info.FDX23 = FDX23;
                info.FDX31 = FDX31;

                info.CY1 = C1 + DX12 * y0 - DY12 * x0;
                info.CY2 = C2 + DX23 * y0 - DY23 * x0;
                info.CY3 = C3 + DX31 * y0 - DY31 * x0;

                info.x = x;
                info.y = y;

                // Accept whole block when totally covered
                if(wc_bpp == 24 || wc_bpp == 32) {
                    if(a == 0xF && b == 0xF && c == 0xF) {
                        DrawTileFilled888<true>(info);
                    } else { // Partially covered
                        DrawTileFilled888<false>(info);
                    }
                } else { // wc_bpp == 15 || wc_bpp == 16
                    if(a == 0xF && b == 0xF && c == 0xF) {
                        DrawTileFilled565<true>(info);
                    } else { // Partially covered
                        DrawTileFilled565<false>(info);
                    }
                }
            }
        }
    }
    wc_colorbuffer->Unlock();
}

bool ComputeCoeffMatrix(const VectorPOD4f& v1, const VectorPOD4f& v2, const VectorPOD4f& v3, MatrixPOD3f& m)
{
    //fesetexceptflag
    /*
    m = Matrix3f(Vector3f(v1.x, v1.y, v1.w),
    		   Vector3f(v2.x, v2.y, v2.w),
    		   Vector3f(v3.x, v3.y, v3.w));
    */
    m[0] = v1.x;
    m[1] = v1.y;
    m[2] = v1.w;
    m[3] = v2.x;
    m[4] = v2.y;
    m[5] = v2.w;
    m[6] = v3.x;
    m[7] = v3.y;
    m[8] = v3.w;

    const float eps = (1.0f / 128.0f);
    float det =
        m[0]*(m[4]*m[8] - m[5]*m[7]) +
        m[1]*(m[5]*m[6] - m[3]*m[8]) +
        m[2]*(m[3]*m[7] - m[4]*m[6]);

    //Degenerate or really small triangles
    //have really tiny determinants or equal
    //zero. We don't render these.

    if(std::abs(det) < 0.125f) {
        return false;
    }

    //Triangles with negative determinants
    //are backfaces, which don't need to be rendered
    if(det < 0.0f) {
        return false;
    }
    // Matrix cofactors
    float c00 = +(m[4]*m[8] - m[5]*m[7]); //4,5,7,8
    float c10 = -(m[3]*m[8] - m[5]*m[6]); //3,5,6,8
    float c20 = +(m[3]*m[7] - m[4]*m[6]); //3,4,6,7

    float c01 = -(m[1]*m[8] - m[2]*m[7]); //1,2,7,8
    float c11 = +(m[0]*m[8] - m[2]*m[6]); //0,2,6,8
    float c21 = -(m[0]*m[7] - m[1]*m[6]); //0,1,6,7

    float c02 = +(m[1]*m[5] - m[2]*m[4]); //1,2,4,5
    float c12 = -(m[0]*m[5] - m[2]*m[3]); //0,2,3,5
    float c22 = +(m[0]*m[4] - m[1]*m[3]); //0,1,3,4

    //Make an adjoint matrix (reuse m)
    m[0] = c00;
    m[1] = c01;
    m[2] = c02;
    m[3] = c10;
    m[4] = c11;
    m[5] = c12;
    m[6] = c20;
    m[7] = c21;
    m[8] = c22;

    det = 1.0f/det;
    m[0] *= det;
    m[1] *= det;
    m[2] *= det;
    m[3] *= det;
    m[4] *= det;
    m[5] *= det;
    m[6] *= det;
    m[7] *= det;
    m[8] *= det;

    return true;
}

inline static void SR_InterpTransform(float& f1, float& f2, float& f3, const MatrixPOD3f& m)
{
    /* Vector made up by one scalar from each vertex */
    VectorPOD3f v = {f1, f2, f3};

    /* Multiply by coefficient matrix */
    v = Mat3Vec3Mul(m, v);
    /* Assign the transformed values back */

    f1 = v.x;
    f2 = v.y;
    f3 = v.z;
}

inline static void SR_InterpTransform(VectorPOD4f& v1, VectorPOD4f& v2, VectorPOD4f& v3, const MatrixPOD3f& m)
{
    SR_InterpTransform(v1.x, v2.x, v3.x, m);
    SR_InterpTransform(v1.y, v2.y, v3.y, m);
    SR_InterpTransform(v1.z, v2.z, v3.z, m);
    SR_InterpTransform(v1.w, v2.w, v3.w, m);
}

void SR_Render(unsigned int flags)
{
    clip_triangle(Vector4f( 1.0f, 0.0f,  0.0f, 1.0f), flags);
    clip_triangle(Vector4f(-1.0f, 0.0f,  0.0f, 1.0f), flags);
    clip_triangle(Vector4f(0.0f,  1.0f,  0.0f, 1.0f), flags);
    clip_triangle(Vector4f(0.0f, -1.0f,  0.0f, 1.0f), flags);
    clip_triangle(Vector4f(0.0f,  0.0f,  1.0f, 1.0f), flags);
    clip_triangle(Vector4f(0.0f,  0.0f, -1.0f, 1.0f), flags);

    size_t oldSize = wc_vertices->size();
    wc_vertices->reserve(oldSize * 2);
    //Do the projection matrix multiply in main() instead, so we can make
    //a big batch of triangles instead of many few.
    /*
    Matrix4f modelviewProjection = wc_projection * wc_modelview;
    for(int i = 0; i < oldSize; i+=3){
      (*wc_vertices)[i+0] = modelviewProjection * (*wc_vertices)[i+0];
      (*wc_vertices)[i+1] = modelviewProjection * (*wc_vertices)[i+1];
      (*wc_vertices)[i+2] = modelviewProjection * (*wc_vertices)[i+2];
    }
    */

    for(int i = 0; i < oldSize; i+=3) {
        /* Compute [a,b,c] coefficients */
        MatrixPOD3f m;
        bool b = ComputeCoeffMatrix((*wc_vertices)[i+0], (*wc_vertices)[i+1], (*wc_vertices)[i+2], m);
        //Skip degenerate triangles, small triangles and backfaces
        if(!b) continue;

        //Project() :
        //Compute screen space coordinates for x and y
        //Normalize z into [0.0f, 1.0f> half-range, Q0.16 fixedpoint
        (*wc_vertices)[i+0] = project((*wc_vertices)[i+0], wc_colorbuffer->w, wc_colorbuffer->h);
        (*wc_vertices)[i+1] = project((*wc_vertices)[i+1], wc_colorbuffer->w, wc_colorbuffer->h);
        (*wc_vertices)[i+2] = project((*wc_vertices)[i+2], wc_colorbuffer->w, wc_colorbuffer->h);

        //Must interpolate z linearly in screenspace!
        //To get the coefficients required for an affine interpolation, simply multiply z with w
        (*wc_vertices)[i+0].z *= (*wc_vertices)[i+0].w;
        (*wc_vertices)[i+1].z *= (*wc_vertices)[i+1].w;
        (*wc_vertices)[i+2].z *= (*wc_vertices)[i+2].w;
        SR_InterpTransform((*wc_vertices)[i+0].z, (*wc_vertices)[i+1].z, (*wc_vertices)[i+2].z, m);

        // To get "1.0f / w", multiply the 3D Vector [1,1,1] with the coefficient matrix.
        // We don't need w (at the triangle vertices) anymore after this point.
        (*wc_vertices)[i+0].w = (*wc_vertices)[i+1].w = (*wc_vertices)[i+2].w = 1.0f;
        SR_InterpTransform((*wc_vertices)[i+0].w, (*wc_vertices)[i+1].w, (*wc_vertices)[i+2].w, m);

        /* Compute the coefficients for the rest of the per-vertex data */
        if(flags & SR_TEXCOORD0)
            SR_InterpTransform((*wc_tcoords0)[i+0], (*wc_tcoords0)[i+1], (*wc_tcoords0)[i+2], m);
        if(flags & SR_TEXCOORD1)
            SR_InterpTransform((*wc_tcoords1)[i+0], (*wc_tcoords1)[i+1], (*wc_tcoords1)[i+2], m);
        if(flags & SR_LIGHTING)
            SR_InterpTransform((*wc_normals)[i+0], (*wc_normals)[i+1], (*wc_normals)[i+2], m);
        if(flags & SR_COLOR)
            SR_InterpTransform((*wc_colors)[i+0], (*wc_colors)[i+1], (*wc_colors)[i+2], m);

        wc_vertices->push_back((*wc_vertices)[i+0]);
        wc_vertices->push_back((*wc_vertices)[i+1]);
        wc_vertices->push_back((*wc_vertices)[i+2]);
        wc_tcoords0->push_back((*wc_tcoords0)[i+0]);
        wc_tcoords0->push_back((*wc_tcoords0)[i+1]);
        wc_tcoords0->push_back((*wc_tcoords0)[i+2]);
    }
    //reuse these arrays but delete the previous data copy
    wc_vertices->erase(wc_vertices->begin(), wc_vertices->begin() + oldSize);
    wc_tcoords0->erase(wc_tcoords0->begin(), wc_tcoords0->begin() + oldSize);

    switch(flags) {
    case SR_TEXCOORD0:
    case SR_TEXCOORD1:
    case SR_LIGHTING:
    case SR_COLOR:
    case (SR_TEXCOORD0 | SR_TEXCOORD1):
    case (SR_TEXCOORD0 | SR_LIGHTING):
    case (SR_TEXCOORD0 | SR_COLOR):
    case (SR_TEXCOORD1 | SR_LIGHTING):
    case (SR_TEXCOORD1 | SR_COLOR):
    case (SR_LIGHTING | SR_COLOR):
    case (SR_TEXCOORD0 | SR_TEXCOORD1 | SR_LIGHTING):
    case (SR_TEXCOORD0 | SR_TEXCOORD1 | SR_COLOR):
    case (SR_TEXCOORD0 | SR_LIGHTING | SR_COLOR):
    case (SR_TEXCOORD1 | SR_LIGHTING | SR_COLOR):
    case (SR_TEXCOORD0 | SR_TEXCOORD1 | SR_LIGHTING | SR_COLOR):
    case 0: {

    }
    }
    //switch not done yet, so just call the test rasterizer
    DrawTriangles(flags);
    //DrawTrianglesDeferred(flags);
}

