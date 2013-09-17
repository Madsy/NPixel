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

static std::vector<TileSet> wc_tileListFilled; //completely filled tiles
static std::vector<TileSet> wc_tileList; //Partially filled

void BlitTilesFilled()
{
    unsigned int* colorbuffer = wc_colorbuffer->Lock();
    unsigned short* depthbuffer = wc_depthbuffer->Ptr();
    const unsigned int width = wc_colorbuffer->w;
    const unsigned int height = wc_colorbuffer->h;
    const unsigned int numTilesX = (width >> Q);
    const unsigned int numTilesY = (height >> Q);
    const int iTw = wc_texture0->width;
    const int iTh = wc_texture0->height;
    const unsigned int* tbuf = &wc_texture0->texels[0];
    for(int y = 0; y < height; y+= q) {
        for(int x = 0; x < width; x+= q) {
            const int tileX = x >> Q;
            const int tileY = y >> Q;
            TileSet& tileSet = wc_tileListFilled[tileX + tileY*numTilesX];
            if(tileSet.tiles.empty()) continue;
#if 0
            //The cover optimization doesn't always work,
            //and very few tiles have 100% cover, so it was
            //a lot of effort for little gain. Disabled for now
            //It's maybe worth to use in very specific scenes, with a high resolution
            //NOTE: The sorting probably spends more time than we save, but the idea
            //also works without sorting, just sub-optimal.
            int zMin0 = 0;
            int zMax0 = 65535;
            std::sort(tileSet.tiles.begin(), tileSet.tiles.begin() + tileSet.count);
#endif
            for(int i = 0; i < tileSet.count; ++i) {
                Tile& t = tileSet.tiles[i];
                int zMin1 = t.zMin;
                int zMax1 = t.zMax;
                bool skipZTest = false;
#if 0
                // z < buffer (65535)
                if(zMin1 > zMax0) {
                    // Occluded anyway, so skip
                    continue;
                } else if(zMax1 > zMax0 && zMin1 < zMax0) {
                    // Intersecting
                    //zMax0 = zMax1;
                    zMin0 = std::min(zMin0, zMin1);
                    skipZTest = false;
                } else if(zMax1 < zMin0) {
                    zMax0 = zMax1;
                    zMin0 = zMin1;
                    // totally at the front, so no need to z test
                    skipZTest = true;
                }
#endif
                //Gradients for y interpolation
                const int bwSlopeY0 = t.bw1 - t.bw0;
                const int bwSlopeY1 = t.bw3 - t.bw2;
                const int bzSlopeY0 = t.bz1 - t.bz0;
                const int bzSlopeY1 = t.bz3 - t.bz2;
                const int buSlopeY0 = t.bu1 - t.bu0;
                const int buSlopeY1 = t.bu3 - t.bu2;
                const int bvSlopeY0 = t.bv1 - t.bv0;
                const int bvSlopeY1 = t.bv3 - t.bv2;
                //Accumulators (actual interpolated value) for y
                int bwSlopeYAccum0 = t.bw0 << Q;
                int bwSlopeYAccum1 = t.bw2 << Q;
                int bzSlopeYAccum0 = t.bz0 << Q;
                int bzSlopeYAccum1 = t.bz2 << Q;
                int buSlopeYAccum0 = t.bu0 << Q;
                int buSlopeYAccum1 = t.bu2 << Q;
                int bvSlopeYAccum0 = t.bv0 << Q;
                int bvSlopeYAccum1 = t.bv2 << Q;
                int col = y*width;
                for(int iy = y; iy < y+q; ++iy) {
                    //Gradients for x interpolation
                    const int bwSlopeX0 = bwSlopeYAccum1 - bwSlopeYAccum0;
                    const int bzSlopeX0 = bzSlopeYAccum1 - bzSlopeYAccum0;
                    const int buSlopeX0 = buSlopeYAccum1 - buSlopeYAccum0;
                    const int bvSlopeX0 = bvSlopeYAccum1 - bvSlopeYAccum0;
                    //Accumulators (actual interpolated value) for x
                    int bwSlopeXAccum0 = bwSlopeYAccum0 << Q;
                    int bzSlopeXAccum0 = bzSlopeYAccum0 << Q;
                    int buSlopeXAccum0 = buSlopeYAccum0 << Q;
                    int bvSlopeXAccum0 = bvSlopeYAccum0 << Q;
                    int fbIndex = x+col;
                    if(skipZTest) {
                        for(int ix = x; ix < x+q; ++ix) {
                            unsigned short z = bzSlopeXAccum0 >> (Q*2);
                            depthbuffer[fbIndex] = z;
                            int uw = buSlopeXAccum0>>(Q*2);
                            int vw = bvSlopeXAccum0>>(Q*2);
                            int w = bwSlopeXAccum0>>(Q*2);
                            int u = ((long long)uw*w*(iTw - 1)) >> (coeff_precision_base * 2);
                            int v = ((long long)vw*w*(iTh - 1)) >> (coeff_precision_base * 2);
                            u = clamp((int)u, 0, iTw-1);
                            v = clamp((int)v, 0, iTh-1);
                            int idxTex = u + v*iTw;
                            colorbuffer[fbIndex] = tbuf[idxTex];
                            ++fbIndex;
                            bwSlopeXAccum0 += bwSlopeX0;
                            bzSlopeXAccum0 += bzSlopeX0;
                            buSlopeXAccum0 += buSlopeX0;
                            bvSlopeXAccum0 += bvSlopeX0;
                        }
                    } else {
                        for(int ix = x; ix < x+q; ++ix) {
                            unsigned short z = bzSlopeXAccum0 >> (Q*2);
                            if(z < depthbuffer[fbIndex]) {
                                depthbuffer[fbIndex] = z;
                                int uw = buSlopeXAccum0>>(Q*2);
                                int vw = bvSlopeXAccum0>>(Q*2);
                                int w = bwSlopeXAccum0>>(Q*2);
                                int u = ((long long)uw*w*(iTw - 1)) >> (coeff_precision_base * 2);
                                int v = ((long long)vw*w*(iTh - 1)) >> (coeff_precision_base * 2);
                                u = clamp((int)u, 0, iTw-1);
                                v = clamp((int)v, 0, iTh-1);
                                int idxTex = u + v*iTw;
                                colorbuffer[fbIndex] = tbuf[idxTex];
                            }
                            ++fbIndex;
                            bwSlopeXAccum0 += bwSlopeX0;
                            bzSlopeXAccum0 += bzSlopeX0;
                            buSlopeXAccum0 += buSlopeX0;
                            bvSlopeXAccum0 += bvSlopeX0;
                        }
                    }
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
            }
        }
    }
    wc_colorbuffer->Unlock();
}

void BlitTiles()
{
    unsigned int* colorbuffer = wc_colorbuffer->Lock();
    unsigned short* depthbuffer = wc_depthbuffer->Ptr();
    const unsigned int width = wc_colorbuffer->w;
    const unsigned int height = wc_colorbuffer->h;
    const unsigned int numTilesX = (width >> Q);
    const unsigned int numTilesY = (height >> Q);
    const int iTw = wc_texture0->width;
    const int iTh = wc_texture0->height;
    const unsigned int* tbuf = &wc_texture0->texels[0];

    for(int y = 0; y < height; y+= q) {
        for(int x = 0; x < width; x+= q) {
            const int tileX = x >> Q;
            const int tileY = y >> Q;
            TileSet& tileSet = wc_tileList[tileX + tileY*numTilesX];
            if(tileSet.tiles.empty()) continue;
            for(int i = 0; i < tileSet.count; ++i) {
                Tile& t = tileSet.tiles[i];
                //Gradients for y interpolation
                const int bwSlopeY0 = t.bw1 - t.bw0;
                const int bwSlopeY1 = t.bw3 - t.bw2;
                const int bzSlopeY0 = t.bz1 - t.bz0;
                const int bzSlopeY1 = t.bz3 - t.bz2;
                const int buSlopeY0 = t.bu1 - t.bu0;
                const int buSlopeY1 = t.bu3 - t.bu2;
                const int bvSlopeY0 = t.bv1 - t.bv0;
                const int bvSlopeY1 = t.bv3 - t.bv2;
                const int FDY12 = t.FDY12;
                const int FDY23 = t.FDY23;
                const int FDY31 = t.FDY31;
                const int FDX12 = t.FDX12;
                const int FDX23 = t.FDX23;
                const int FDX31 = t.FDX31;
                //Accumulators (actual interpolated value) for y
                int bwSlopeYAccum0 = t.bw0 << Q;
                int bwSlopeYAccum1 = t.bw2 << Q;
                int bzSlopeYAccum0 = t.bz0 << Q;
                int bzSlopeYAccum1 = t.bz2 << Q;
                int buSlopeYAccum0 = t.bu0 << Q;
                int buSlopeYAccum1 = t.bu2 << Q;
                int bvSlopeYAccum0 = t.bv0 << Q;
                int bvSlopeYAccum1 = t.bv2 << Q;
                int CY1 = t.CY1;
                int CY2 = t.CY2;
                int CY3 = t.CY3;
                int col = y*width;
                for(int iy = y; iy < y+q; ++iy) {
                    int CX1 = CY1;
                    int CX2 = CY2;
                    int CX3 = CY3;
                    //Gradients for x interpolation
                    const int bwSlopeX0 = bwSlopeYAccum1 - bwSlopeYAccum0;
                    const int bzSlopeX0 = bzSlopeYAccum1 - bzSlopeYAccum0;
                    const int buSlopeX0 = buSlopeYAccum1 - buSlopeYAccum0;
                    const int bvSlopeX0 = bvSlopeYAccum1 - bvSlopeYAccum0;
                    //Accumulators (actual interpolated value) for x
                    int bwSlopeXAccum0 = bwSlopeYAccum0 << Q;
                    int bzSlopeXAccum0 = bzSlopeYAccum0 << Q;
                    int buSlopeXAccum0 = buSlopeYAccum0 << Q;
                    int bvSlopeXAccum0 = bvSlopeYAccum0 << Q;
                    int fbIndex = x + col;
                    for(int ix = x; ix < x+q; ++ix) {
                        if(CX1 > 0 && CX2 > 0 && CX3 > 0) {
                            unsigned short z = bzSlopeXAccum0 >> (Q*2);
                            if(z < depthbuffer[fbIndex]) {
                                depthbuffer[fbIndex] = z;
                                int uw = buSlopeXAccum0>>(Q*2);
                                int vw = bvSlopeXAccum0>>(Q*2);
                                int w = bwSlopeXAccum0>>(Q*2);
                                int u = ((long long)uw*w*(iTw - 1)) >> (coeff_precision_base * 2);
                                int v = ((long long)vw*w*(iTh - 1)) >> (coeff_precision_base * 2);
                                u = clamp((int)u, 0, iTw-1);
                                v = clamp((int)v, 0, iTh-1);
                                int idxTex = u + v*iTw;
                                colorbuffer[fbIndex] = tbuf[idxTex];
                            }
                        }
                        ++fbIndex;
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
    wc_colorbuffer->Unlock();
}


/************************************************ Testing use of a mini z-bufffer ***************************/
/*
int zmin, zmax
zmin = 65535
zmax = 0

for(every y)
  for(every x)

 */
void BlitTilesCombined()
{
    // 16x16 z-buffer! :-)
    // We only need the buffer to be equal to the tile size
    unsigned short zbuffer[q*q];
    unsigned int* colorbuffer = wc_colorbuffer->Lock();
    unsigned short* depthbuffer = zbuffer;
    const unsigned int width = wc_colorbuffer->w;
    const unsigned int height = wc_colorbuffer->h;
    const unsigned int numTilesX = (width >> Q);
    const unsigned int numTilesY = (height >> Q);
    const int iTw = wc_texture0->width;
    const int iTh = wc_texture0->height;
    const unsigned int* tbuf = &wc_texture0->texels[0];
    for(int y = 0; y < height; y+= q) {
        for(int x = 0; x < width; x+= q) {
            //Here we reset the 16x16 z-buffer for every tile X,Y position
            memset(depthbuffer, 0xFF, sizeof(unsigned short) * q * q);
            const int tileX = x >> Q;
            const int tileY = y >> Q;
            TileSet& tileSet = wc_tileList[tileX + tileY*numTilesX];
            TileSet& tileSetFilled = wc_tileListFilled[tileX + tileY*numTilesX];
            for(int i = 0; i < tileSetFilled.count; ++i) {
                Tile& t = tileSetFilled.tiles[i];
                //Gradients for y interpolation
                const int bwSlopeY0 = t.bw1 - t.bw0;
                const int bwSlopeY1 = t.bw3 - t.bw2;
                const int bzSlopeY0 = t.bz1 - t.bz0;
                const int bzSlopeY1 = t.bz3 - t.bz2;
                const int buSlopeY0 = t.bu1 - t.bu0;
                const int buSlopeY1 = t.bu3 - t.bu2;
                const int bvSlopeY0 = t.bv1 - t.bv0;
                const int bvSlopeY1 = t.bv3 - t.bv2;
                const int FDY12 = t.FDY12;
                const int FDY23 = t.FDY23;
                const int FDY31 = t.FDY31;
                const int FDX12 = t.FDX12;
                const int FDX23 = t.FDX23;
                const int FDX31 = t.FDX31;
                //Accumulators (actual interpolated value) for y
                int bwSlopeYAccum0 = t.bw0 << Q;
                int bwSlopeYAccum1 = t.bw2 << Q;
                int bzSlopeYAccum0 = t.bz0 << Q;
                int bzSlopeYAccum1 = t.bz2 << Q;
                int buSlopeYAccum0 = t.bu0 << Q;
                int buSlopeYAccum1 = t.bu2 << Q;
                int bvSlopeYAccum0 = t.bv0 << Q;
                int bvSlopeYAccum1 = t.bv2 << Q;
                int CY1 = t.CY1;
                int CY2 = t.CY2;
                int CY3 = t.CY3;
                int col = y*width;
                for(int iy = 0; iy < q; ++iy) {
                    int CX1 = CY1;
                    int CX2 = CY2;
                    int CX3 = CY3;
                    //Gradients for x interpolation
                    const int bwSlopeX0 = bwSlopeYAccum1 - bwSlopeYAccum0;
                    const int bzSlopeX0 = bzSlopeYAccum1 - bzSlopeYAccum0;
                    const int buSlopeX0 = buSlopeYAccum1 - buSlopeYAccum0;
                    const int bvSlopeX0 = bvSlopeYAccum1 - bvSlopeYAccum0;
                    //Accumulators (actual interpolated value) for x
                    int bwSlopeXAccum0 = bwSlopeYAccum0 << Q;
                    int bzSlopeXAccum0 = bzSlopeYAccum0 << Q;
                    int buSlopeXAccum0 = buSlopeYAccum0 << Q;
                    int bvSlopeXAccum0 = bvSlopeYAccum0 << Q;
                    int fbIndex = x + col;
                    int zIndex = iy*q;
                    for(int ix = 0; ix < q; ++ix) {
                        unsigned short z = bzSlopeXAccum0 >> (Q*2);
                        if(z < depthbuffer[zIndex]) {
                            depthbuffer[zIndex] = z;
                            int uw = buSlopeXAccum0>>(Q*2);
                            int vw = bvSlopeXAccum0>>(Q*2);
                            int w = bwSlopeXAccum0>>(Q*2);
                            int u = ((long long)uw*w*(iTw - 1)) >> (coeff_precision_base * 2);
                            int v = ((long long)vw*w*(iTh - 1)) >> (coeff_precision_base * 2);
                            u = clamp((int)u, 0, iTw-1);
                            v = clamp((int)v, 0, iTh-1);
                            int idxTex = u + v*iTw;
                            colorbuffer[fbIndex] = tbuf[idxTex];
                        }
                        ++fbIndex;
                        ++zIndex;
                        bwSlopeXAccum0 += bwSlopeX0;
                        bzSlopeXAccum0 += bzSlopeX0;
                        buSlopeXAccum0 += buSlopeX0;
                        bvSlopeXAccum0 += bvSlopeX0;
                        CX1 -= FDY12;
                        CX2 -= FDY23;
                        CX3 -= FDY31;
                    } //x loop
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
                } //y loop
            } // tileset loop
            for(int i = 0; i < tileSet.count; ++i) {
                Tile& t = tileSet.tiles[i];
                //Gradients for y interpolation
                const int bwSlopeY0 = t.bw1 - t.bw0;
                const int bwSlopeY1 = t.bw3 - t.bw2;
                const int bzSlopeY0 = t.bz1 - t.bz0;
                const int bzSlopeY1 = t.bz3 - t.bz2;
                const int buSlopeY0 = t.bu1 - t.bu0;
                const int buSlopeY1 = t.bu3 - t.bu2;
                const int bvSlopeY0 = t.bv1 - t.bv0;
                const int bvSlopeY1 = t.bv3 - t.bv2;
                const int FDY12 = t.FDY12;
                const int FDY23 = t.FDY23;
                const int FDY31 = t.FDY31;
                const int FDX12 = t.FDX12;
                const int FDX23 = t.FDX23;
                const int FDX31 = t.FDX31;
                //Accumulators (actual interpolated value) for y
                int bwSlopeYAccum0 = t.bw0 << Q;
                int bwSlopeYAccum1 = t.bw2 << Q;
                int bzSlopeYAccum0 = t.bz0 << Q;
                int bzSlopeYAccum1 = t.bz2 << Q;
                int buSlopeYAccum0 = t.bu0 << Q;
                int buSlopeYAccum1 = t.bu2 << Q;
                int bvSlopeYAccum0 = t.bv0 << Q;
                int bvSlopeYAccum1 = t.bv2 << Q;
                int CY1 = t.CY1;
                int CY2 = t.CY2;
                int CY3 = t.CY3;
                int col = y*width;
                for(int iy = 0; iy < q; ++iy) {
                    int CX1 = CY1;
                    int CX2 = CY2;
                    int CX3 = CY3;
                    //Gradients for x interpolation
                    const int bwSlopeX0 = bwSlopeYAccum1 - bwSlopeYAccum0;
                    const int bzSlopeX0 = bzSlopeYAccum1 - bzSlopeYAccum0;
                    const int buSlopeX0 = buSlopeYAccum1 - buSlopeYAccum0;
                    const int bvSlopeX0 = bvSlopeYAccum1 - bvSlopeYAccum0;
                    //Accumulators (actual interpolated value) for x
                    int bwSlopeXAccum0 = bwSlopeYAccum0 << Q;
                    int bzSlopeXAccum0 = bzSlopeYAccum0 << Q;
                    int buSlopeXAccum0 = buSlopeYAccum0 << Q;
                    int bvSlopeXAccum0 = bvSlopeYAccum0 << Q;
                    int fbIndex = x + col;
                    int zIndex = iy*q;
                    for(int ix = 0; ix < q; ++ix) {
                        if(CX1 > 0 && CX2 > 0 && CX3 > 0) {
                            unsigned short z = bzSlopeXAccum0 >> (Q*2);
                            if(z < depthbuffer[zIndex]) {
                                depthbuffer[zIndex] = z;
                                int uw = buSlopeXAccum0>>(Q*2);
                                int vw = bvSlopeXAccum0>>(Q*2);
                                int w = bwSlopeXAccum0>>(Q*2);
                                int u = ((long long)uw*w*(iTw - 1)) >> (coeff_precision_base * 2);
                                int v = ((long long)vw*w*(iTh - 1)) >> (coeff_precision_base * 2);
                                u = clamp((int)u, 0, iTw-1);
                                v = clamp((int)v, 0, iTh-1);
                                int idxTex = u + v*iTw;
                                colorbuffer[fbIndex] = tbuf[idxTex];
                            }
                        }
                        ++fbIndex;
                        ++zIndex;
                        bwSlopeXAccum0 += bwSlopeX0;
                        bzSlopeXAccum0 += bzSlopeX0;
                        buSlopeXAccum0 += buSlopeX0;
                        bvSlopeXAccum0 += bvSlopeX0;
                        CX1 -= FDY12;
                        CX2 -= FDY23;
                        CX3 -= FDY31;
                    } //x loop
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
                } //y loop
            } // tileset loop
        } //x tile loop
    } //y tile loop
    wc_colorbuffer->Unlock();
}







void DrawTrianglesDeferred(unsigned int flags)
{
    using std::min;
    using std::max;

    const unsigned int width = wc_colorbuffer->w;
    const unsigned int height = wc_colorbuffer->h;
    const unsigned int numTilesX = (width >> Q);
    const unsigned int numTilesY = (height >> Q);
    const float fHalfWidthInv = 2.0f / (float)width;
    const float fHalfHeightInv = 2.0f / (float)height;
    const int NDC_x_step = fHalfWidthInv * f_ndc_precision; //1 subtracted later
    const int NDC_y_step = fHalfHeightInv * f_ndc_precision; //1 subtracted later

    if(wc_tileList.size() < (numTilesX * numTilesY))
        wc_tileList.resize(numTilesX * numTilesY);
    if(wc_tileListFilled.size() < (numTilesX * numTilesY))
        wc_tileListFilled.resize(numTilesX * numTilesY);

    for(int i = 0; i < wc_tileList.size(); ++i) {
        wc_tileList[i].count = 0;
    }
    for(int i = 0; i < wc_tileListFilled.size(); ++i) {
        wc_tileListFilled[i].count = 0;
    }

    const VectorPOD4f* vertices = &((*wc_vertices)[0]);
    const VectorPOD4f* tcoords = &((*wc_tcoords0)[0]);

    size_t len = wc_vertices->size();

    for(int i=0; i<wc_vertices->size(); i+=3) {
        const VectorPOD4f& v1 = vertices[i+0];
        const VectorPOD4f& v2 = vertices[i+2];
        const VectorPOD4f& v3 = vertices[i+1];

        const VectorPOD4f& tc1 = tcoords[i+0];
        const VectorPOD4f& tc2 = tcoords[i+1];
        const VectorPOD4f& tc3 = tcoords[i+2];

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

                int bwi0 = (((long long)Aw*bwx0 + Bw*bwy0) >> coeff_precision_base) + Cw; //top left
                int bwi1 = (((long long)Aw*bwx0 + Bw*bwy1) >> coeff_precision_base) + Cw; //bottom left
                int bwi2 = (((long long)Aw*bwx1 + Bw*bwy0) >> coeff_precision_base) + Cw; //top right
                int bwi3 = (((long long)Aw*bwx1 + Bw*bwy1) >> coeff_precision_base) + Cw; //bottom right

                int bz0 = (((long long)Az*bzx0 + Bz*bzy0) >> depth_precision_base) + Cz; //top left
                int bz1 = (((long long)Az*bzx0 + Bz*bzy1) >> depth_precision_base) + Cz; //bottom left
                int bz2 = (((long long)Az*bzx1 + Bz*bzy0) >> depth_precision_base) + Cz; //top right
                int bz3 = (((long long)Az*bzx1 + Bz*bzy1) >> depth_precision_base) + Cz; //bottom right

                //Compute u for the corners of the tile
                int bu0 = (((long long)Au*bwx0 + Bu*bwy0) >> coeff_precision_base) + Cu; //top left
                int bu1 = (((long long)Au*bwx0 + Bu*bwy1) >> coeff_precision_base) + Cu; //bottom left
                int bu2 = (((long long)Au*bwx1 + Bu*bwy0) >> coeff_precision_base) + Cu; //top right
                int bu3 = (((long long)Au*bwx1 + Bu*bwy1) >> coeff_precision_base) + Cu; //bottom right

                //Compute v for the corners of the tile
                int bv0 = (((long long)Av*bwx0 + Bv*bwy0) >> coeff_precision_base) + Cv; //top left
                int bv1 = (((long long)Av*bwx0 + Bv*bwy1) >> coeff_precision_base) + Cv; //bottom left
                int bv2 = (((long long)Av*bwx1 + Bv*bwy0) >> coeff_precision_base) + Cv; //top right
                int bv3 = (((long long)Av*bwx1 + Bv*bwy1) >> coeff_precision_base) + Cv; //bottom right

                int bw0, bw1, bw2, bw3;
                bw0 = bw1 = bw2 = bw3 = 0;

                if(bwi0) bw0 = (1<<(coeff_precision_base * 2)) / bwi0;
                if(bwi1) bw1 = (1<<(coeff_precision_base * 2)) / bwi1;
                if(bwi2) bw2 = (1<<(coeff_precision_base * 2)) / bwi2;
                if(bwi3) bw3 = (1<<(coeff_precision_base * 2)) / bwi3;

                Tile tile;

                tile.FDX12 = FDX12;
                tile.FDX23 = FDX23;
                tile.FDX31 = FDX31;
                tile.FDY12 = FDY12;
                tile.FDY23 = FDY23;
                tile.FDY31 = FDY31;
                tile.CY1 = C1 + DX12 * y0 - DY12 * x0;
                tile.CY2 = C2 + DX23 * y0 - DY23 * x0;
                tile.CY3 = C3 + DX31 * y0 - DY31 * x0;
                tile.bw0 = bw0;
                tile.bw1 = bw1;
                tile.bw2 = bw2;
                tile.bw3 = bw3;
                tile.bz0 = bz0;
                tile.bz1 = bz1;
                tile.bz2 = bz2;
                tile.bz3 = bz3;
                tile.bu0 = bu0;
                tile.bu1 = bu1;
                tile.bu2 = bu2;
                tile.bu3 = bu3;
                tile.bv0 = bv0;
                tile.bv1 = bv1;
                tile.bv2 = bv2;
                tile.bv3 = bv3;
                tile.zMin = std::min(std::min(std::min(bz0, bz1), bz2), bz3);
                tile.zMax = std::max(std::max(std::max(bz0, bz1), bz2), bz3);

                int tileIdx = (x >> Q) + (y >> Q) * numTilesX;
                // Accept whole block when totally covered
                if(a == 0xF && b == 0xF && c == 0xF) {
                    int tileListIdx = wc_tileListFilled[tileIdx].count;
                    wc_tileListFilled[tileIdx].tiles[tileListIdx] = tile;
                    wc_tileListFilled[tileIdx].count++;
                } else {
                    int tileListIdx = wc_tileList[tileIdx].count;
                    wc_tileList[tileIdx].tiles[tileListIdx] = tile;
                    wc_tileList[tileIdx].count++;
                }
            }
        }
    }
    //BlitTilesFilled();
    //BlitTiles();
    BlitTilesCombined();
}

