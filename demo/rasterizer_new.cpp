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
#include "shaders.h"


//#define PASSMODE //Fill-color blit-loop for testing
bool tTest = false;
//Fixedpoint base for coefficients
const int coeff_precision_base = 11;
//Fixedpoint base for NDC coordinates (we convert screen-space coords to NDCs later)
const int ndc_precision_base = 20;
//Need 16-bits precision for z-buffer
const int depth_precision_base = 16;
//The actual scalars for the bases above (precision = 1 << base)
const float f_coeff_precision = (float)(1 << coeff_precision_base);
const float f_ndc_precision = (float)(1 << ndc_precision_base);
const float f_depth_precision = (float)(1 << depth_precision_base);
const int i_coeff_precision = 1 << coeff_precision_base;
const int i_ndc_precision = 1 << ndc_precision_base;
const int i_depth_precision = 1 << depth_precision_base;
//Base to use when converting from NDC coordinates to coefficients
const int base_diff = (ndc_precision_base - coeff_precision_base);
//Base to use when converting from NDC coordinates to depth
const int base_diff_z = (ndc_precision_base - depth_precision_base);


struct Tile
{
  int FDX12, FDX23, FDX31;
  int FDY12, FDY23, FDY31;
  int CY1, CY2, CY3;
  int bw0, bw1, bw2, bw3; //w corner values
  int bz0, bz1, bz2, bz3; //w corner values
  int bu0, bu1, bu2, bu3; //w corner values
  int bv0, bv1, bv2, bv3; //w corner values
  int zMin, zMax; //for early z-culling
  bool operator<(const Tile& t) const{
	return zMin < t.zMin;
  }
};

struct TileSet
{
  TileSet() : count(0) {
	tiles.resize(500);
  }
  int count;
  std::vector<Tile> tiles;
};

static std::vector<TileSet> wc_tileListFilled; //completely filled tiles
static std::vector<TileSet> wc_tileList; //Partially filled

const int Q = 4;
const int q = (1<<Q);

int tileFillCount = 0;
int tileCount = 0;

void DrawTriangles(unsigned int flags)
{
  using std::min;
  using std::max;

  for(int i=0; i<wc_vertices->size(); i+=3){
    Vector4f& v1 = (*wc_vertices)[i+0];
    Vector4f& v2 = (*wc_vertices)[i+2];
    Vector4f& v3 = (*wc_vertices)[i+1];

    Vector4f& tc1 = (*wc_tcoords0)[i+0];
    Vector4f& tc2 = (*wc_tcoords0)[i+1];
    Vector4f& tc3 = (*wc_tcoords0)[i+2];

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

	// Block size, 16x16. Must be POT
	const int Q = 4;
	const int q = (1<<Q);

	// Start in corner of a 8x8 block alligned to block size
	minx &= ~(q - 1);
	miny &= ~(q - 1);
	maxx = (maxx + (q - 1)) & ~(q - 1);
	maxy = (maxy + (q - 1)) & ~(q - 1);

	unsigned int width = wc_colorbuffer.w;
	unsigned int height = wc_colorbuffer.h;

	// Half-edge constants
	int C1 = DY12 * X1 - DX12 * Y1;
	int C2 = DY23 * X2 - DX23 * Y2;
	int C3 = DY31 * X3 - DX31 * Y3;

	// Correct for fill convention
	if(DY12 < 0 || (DY12 == 0 && DX12 > 0)) C1++;
	if(DY23 < 0 || (DY23 == 0 && DX23 > 0)) C2++;
	if(DY31 < 0 || (DY31 == 0 && DX31 > 0)) C3++;

	// Loop through blocks
	for(int y = miny; y < maxy; y += q){
	  for(int x = minx; x < maxx; x += q){

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

		int bwi0 = (((long long)Aw*bwx0)>>coeff_precision_base) + (((long long)Bw*bwy0)>>coeff_precision_base) + Cw; //top left
		int bwi1 = (((long long)Aw*bwx0)>>coeff_precision_base) + (((long long)Bw*bwy1)>>coeff_precision_base) + Cw; //bottom left
		int bwi2 = (((long long)Aw*bwx1)>>coeff_precision_base) + (((long long)Bw*bwy0)>>coeff_precision_base) + Cw; //top right		
		int bwi3 = (((long long)Aw*bwx1)>>coeff_precision_base) + (((long long)Bw*bwy1)>>coeff_precision_base) + Cw; //bottom right

		int bz0 = (((long long)Az*bzx0)>>depth_precision_base) + (((long long)Bz*bzy0)>>depth_precision_base) + Cz; //top left
		int bz1 = (((long long)Az*bzx0)>>depth_precision_base) + (((long long)Bz*bzy1)>>depth_precision_base) + Cz; //bottom left
		int bz2 = (((long long)Az*bzx1)>>depth_precision_base) + (((long long)Bz*bzy0)>>depth_precision_base) + Cz; //top right		
		int bz3 = (((long long)Az*bzx1)>>depth_precision_base) + (((long long)Bz*bzy1)>>depth_precision_base) + Cz; //bottom right

		//Compute u for the corners of the tile
		int bu0 = (((long long)Au*bwx0)>>coeff_precision_base) + (((long long)Bu*bwy0)>>coeff_precision_base) + Cu; //top left
		int bu1 = (((long long)Au*bwx0)>>coeff_precision_base) + (((long long)Bu*bwy1)>>coeff_precision_base) + Cu; //bottom left
		int bu2 = (((long long)Au*bwx1)>>coeff_precision_base) + (((long long)Bu*bwy0)>>coeff_precision_base) + Cu; //top right		
		int bu3 = (((long long)Au*bwx1)>>coeff_precision_base) + (((long long)Bu*bwy1)>>coeff_precision_base) + Cu; //bottom right

		//Compute v for the corners of the tile
		int bv0 = (((long long)Av*bwx0)>>coeff_precision_base) + (((long long)Bv*bwy0)>>coeff_precision_base) + Cv; //top left
		int bv1 = (((long long)Av*bwx0)>>coeff_precision_base) + (((long long)Bv*bwy1)>>coeff_precision_base) + Cv; //bottom left
		int bv2 = (((long long)Av*bwx1)>>coeff_precision_base) + (((long long)Bv*bwy0)>>coeff_precision_base) + Cv; //top right		
		int bv3 = (((long long)Av*bwx1)>>coeff_precision_base) + (((long long)Bv*bwy1)>>coeff_precision_base) + Cv; //bottom right

		int bw0, bw1, bw2, bw3;
		bw0 = bw1 = bw2 = bw3 = 0;
		
		/* bw0, bw1 .. = Box W */
		if(bwi0)
		  bw0 = ((int)1<<(coeff_precision_base * 2)) / bwi0;
		if(bwi1)
		  bw1 = ((int)1<<(coeff_precision_base * 2)) / bwi1;
		if(bwi2)
		  bw2 = ((int)1<<(coeff_precision_base * 2)) / bwi2;
		if(bwi3)
		  bw3 = ((int)1<<(coeff_precision_base * 2)) / bwi3;
				
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

		// Accept whole block when totally covered
		if(a == 0xF && b == 0xF && c == 0xF){
		  unsigned int col = y*width;
		  const unsigned int* tbuf = &wc_texture0->texels[0];
		  int iTw = wc_texture0->width;
		  int iTh = wc_texture0->height;
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

		  for(int iy = y; iy < y + q; iy++){
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

			for(int ix = x; ix < x + q; ix++){
			  unsigned short z = bzSlopeXAccum0 >> (Q*2);
			  if(z < wc_depthbuffer[ix + col]){
				wc_depthbuffer[ix + col] = z;
                #ifndef PASSMODE
				int uw = buSlopeXAccum0>>(Q*2);
				int vw = bvSlopeXAccum0>>(Q*2);
				int w = bwSlopeXAccum0>>(Q*2);
				int u = ((long long)uw*w*(iTw - 1)) >> (coeff_precision_base * 2);
				int v = ((long long)vw*w*(iTh - 1)) >> (coeff_precision_base * 2);
				u = clamp((int)u, 0, iTw-1);
				v = clamp((int)v, 0, iTh-1);
				int idxTex = u + v*iTw;
				wc_colorbuffer[ix + col] = tbuf[idxTex];
				#else
				wc_colorbuffer[ix + col] = 0xFFFF0000;//tbuf[idxTex];
				#endif
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
		  const unsigned int* tbuf = &wc_texture0->texels[0];
		  int iTw = wc_texture0->width;
		  int iTh = wc_texture0->height;
		  int bwSlopeYAccum0 = bw0 << Q; //start vertical w left
		  int bwSlopeYAccum1 = bw2 << Q; //start vertical w right
		  int bzSlopeYAccum0 = bz0 << Q; //start vertical z left
		  int bzSlopeYAccum1 = bz2 << Q; //start vertical z right
		  int buSlopeYAccum0 = bu0 << Q; //start vertical u left
		  int buSlopeYAccum1 = bu2 << Q; //start vertical u right
		  int bvSlopeYAccum0 = bv0 << Q; //start vertical v left
		  int bvSlopeYAccum1 = bv2 << Q; //start vertical v right

		  for(int iy = y; iy < y + q; iy++){
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

			for(int ix = x; ix < x + q; ix++){
			  if(CX1 > 0 && CX2 > 0 && CX3 > 0){				
				unsigned short z = bzSlopeXAccum0 >> (Q*2);
				if(z < wc_depthbuffer[ix + col]){
				  wc_depthbuffer[ix + col] = z;
                  #ifndef PASSMODE
				  int uw = buSlopeXAccum0>>(Q*2);
				  int vw = bvSlopeXAccum0>>(Q*2);
				  int w = bwSlopeXAccum0>>(Q*2);				  
				  int u = ((long long)uw*w*(iTw - 1)) >> (coeff_precision_base * 2);
				  int v = ((long long)vw*w*(iTh - 1)) >> (coeff_precision_base * 2);
				  u = clamp((int)u, 0, iTw-1);
				  v = clamp((int)v, 0, iTh-1);
				  int idxTex = u + v*iTw;
				  wc_colorbuffer[ix + col] = tbuf[idxTex];
				  #else
				  wc_colorbuffer[ix + col] = 0xFF0000FF;//tbuf[idxTex];
				  #endif
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
}

void BlitTilesFilled()
{
#if 1
  unsigned int width = wc_colorbuffer.w;
  unsigned int height = wc_colorbuffer.h;
  unsigned int numTilesX = (width >> Q);
  unsigned int numTilesY = (height >> Q);
  int iTw = wc_texture0->width;
  int iTh = wc_texture0->height;
  const unsigned int* tbuf = &wc_texture0->texels[0];
  for(int y = 0; y < height; y+= q){
	for(int x = 0; x < width; x+= q){
	  int tileX = x >> Q;
	  int tileY = y >> Q;
	  TileSet& tileSet = wc_tileListFilled[tileX + tileY*numTilesX];
	  if(tileSet.tiles.empty()) continue;

	  int zMin0 = 0;
	  int zMax0 = 65535;
	  std::sort(tileSet.tiles.begin(), tileSet.tiles.begin() + tileSet.count);

	  for(int i = 0; i < tileSet.count; ++i){
		Tile& t = tileSet.tiles[i];
		int zMin1 = t.zMin;
		int zMax1 = t.zMax;
		bool skipZTest = false;

		// z < buffer (65535)
		if(zMin1 > zMax0){
		  // Occluded anyway, so skip
		  continue;
		} else if(zMax1 > zMax0 && zMin1 < zMax0){
		  // Intersecting
		  //zMax0 = zMax1;
		  zMin0 = std::min(zMin0, zMin1);
		  skipZTest = false;
		} else if(zMax1 < zMin0){
		  zMax0 = zMax1;
		  zMin0 = zMin1;
		  // totally at the front, so no need to z test
		  skipZTest = true;
		}

		//Gradients for y interpolation
		int bwSlopeY0 = t.bw1 - t.bw0;
		int bwSlopeY1 = t.bw3 - t.bw2;
		int bzSlopeY0 = t.bz1 - t.bz0;
		int bzSlopeY1 = t.bz3 - t.bz2;
		int buSlopeY0 = t.bu1 - t.bu0;
		int buSlopeY1 = t.bu3 - t.bu2;
		int bvSlopeY0 = t.bv1 - t.bv0;
		int bvSlopeY1 = t.bv3 - t.bv2;
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
		for(int iy = y; iy < y+q; ++iy){
		  //Gradients for x interpolation
		  int bwSlopeX0 = bwSlopeYAccum1 - bwSlopeYAccum0;
		  int bzSlopeX0 = bzSlopeYAccum1 - bzSlopeYAccum0;
		  int buSlopeX0 = buSlopeYAccum1 - buSlopeYAccum0;
		  int bvSlopeX0 = bvSlopeYAccum1 - bvSlopeYAccum0;
		  //Accumulators (actual interpolated value) for x
		  int bwSlopeXAccum0 = bwSlopeYAccum0 << Q;
		  int bzSlopeXAccum0 = bzSlopeYAccum0 << Q;
		  int buSlopeXAccum0 = buSlopeYAccum0 << Q;
		  int bvSlopeXAccum0 = bvSlopeYAccum0 << Q;
		  int fbIndex = col;
		  if(skipZTest){
			for(int ix = x; ix < x+q; ++ix){
			  unsigned short z = bzSlopeXAccum0 >> (Q*2);
			  wc_depthbuffer[ix + col] = z;
			  int uw = buSlopeXAccum0>>(Q*2);
			  int vw = bvSlopeXAccum0>>(Q*2);
			  int w = bwSlopeXAccum0>>(Q*2);
			  int u = ((long long)uw*w*(iTw - 1)) >> (coeff_precision_base * 2);
			  int v = ((long long)vw*w*(iTh - 1)) >> (coeff_precision_base * 2);
			  u = clamp((int)u, 0, iTw-1);
			  v = clamp((int)v, 0, iTh-1);
			  int idxTex = u + v*iTw;
			  wc_colorbuffer[ix + col] = tbuf[idxTex];
			}
			++fbIndex;
			bwSlopeXAccum0 += bwSlopeX0;
			bzSlopeXAccum0 += bzSlopeX0;
			buSlopeXAccum0 += buSlopeX0;
			bvSlopeXAccum0 += bvSlopeX0;
		  } else {
			for(int ix = x; ix < x+q; ++ix){
			  unsigned short z = bzSlopeXAccum0 >> (Q*2);
			  if(z < wc_depthbuffer[ix + col]){
				wc_depthbuffer[ix + col] = z;
				int uw = buSlopeXAccum0>>(Q*2);
				int vw = bvSlopeXAccum0>>(Q*2);
				int w = bwSlopeXAccum0>>(Q*2);
				int u = ((long long)uw*w*(iTw - 1)) >> (coeff_precision_base * 2);
				int v = ((long long)vw*w*(iTh - 1)) >> (coeff_precision_base * 2);
				u = clamp((int)u, 0, iTw-1);
				v = clamp((int)v, 0, iTh-1);
				int idxTex = u + v*iTw;
				wc_colorbuffer[ix + col] = tbuf[idxTex];
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
#endif
}

void BlitTiles()
{
#if 1
  unsigned int width = wc_colorbuffer.w;
  unsigned int height = wc_colorbuffer.h;
  unsigned int numTilesX = (width >> Q);
  unsigned int numTilesY = (height >> Q);
  int iTw = wc_texture0->width;
  int iTh = wc_texture0->height;
  const unsigned int* tbuf = &wc_texture0->texels[0];
  for(int y = 0; y < height; y+= q){
	for(int x = 0; x < width; x+= q){
	  int tileX = x >> Q;
	  int tileY = y >> Q;
	  TileSet& tileSet = wc_tileList[tileX + tileY*numTilesX];
	  if(tileSet.tiles.empty()) continue;
	  for(int i = 0; i < tileSet.count; ++i){
		Tile& t = tileSet.tiles[i];
		//Gradients for y interpolation
		int bwSlopeY0 = t.bw1 - t.bw0;
		int bwSlopeY1 = t.bw3 - t.bw2;
		int bzSlopeY0 = t.bz1 - t.bz0;
		int bzSlopeY1 = t.bz3 - t.bz2;
		int buSlopeY0 = t.bu1 - t.bu0;
		int buSlopeY1 = t.bu3 - t.bu2;
		int bvSlopeY0 = t.bv1 - t.bv0;
		int bvSlopeY1 = t.bv3 - t.bv2;
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
		for(int iy = y; iy < y+q; ++iy){
		  int CX1 = CY1;
		  int CX2 = CY2;
		  int CX3 = CY3;
		  //Gradients for x interpolation
		  int bwSlopeX0 = bwSlopeYAccum1 - bwSlopeYAccum0;
		  int bzSlopeX0 = bzSlopeYAccum1 - bzSlopeYAccum0;
		  int buSlopeX0 = buSlopeYAccum1 - buSlopeYAccum0;
		  int bvSlopeX0 = bvSlopeYAccum1 - bvSlopeYAccum0;
		  //Accumulators (actual interpolated value) for x
		  int bwSlopeXAccum0 = bwSlopeYAccum0 << Q;
		  int bzSlopeXAccum0 = bzSlopeYAccum0 << Q;
		  int buSlopeXAccum0 = buSlopeYAccum0 << Q;
		  int bvSlopeXAccum0 = bvSlopeYAccum0 << Q;
		  int fbIndex = col;
		  for(int ix = x; ix < x+q; ++ix){
			if(CX1 > 0 && CX2 > 0 && CX3 > 0){
			  unsigned short z = bzSlopeXAccum0 >> (Q*2);
			  if(z < wc_depthbuffer[ix + col]){
				wc_depthbuffer[ix + col] = z;
				int uw = buSlopeXAccum0>>(Q*2);
				int vw = bvSlopeXAccum0>>(Q*2);
				int w = bwSlopeXAccum0>>(Q*2);
				int u = ((long long)uw*w*(iTw - 1)) >> (coeff_precision_base * 2);
				int v = ((long long)vw*w*(iTh - 1)) >> (coeff_precision_base * 2);
				u = clamp((int)u, 0, iTw-1);
				v = clamp((int)v, 0, iTh-1);
				int idxTex = u + v*iTw;
				wc_colorbuffer[ix + col] = tbuf[idxTex];
			  }
			}
			++fbIndex;
			bwSlopeXAccum0 += bwSlopeX0;
			bzSlopeXAccum0 += bzSlopeX0;
			buSlopeXAccum0 += buSlopeX0;
			bvSlopeXAccum0 += bvSlopeX0;
			CX1 -= t.FDY12;
		    CX2 -= t.FDY23;
		    CX3 -= t.FDY31;
		  }
		  bwSlopeYAccum0 += bwSlopeY0;
		  bwSlopeYAccum1 += bwSlopeY1;
		  bzSlopeYAccum0 += bzSlopeY0;
		  bzSlopeYAccum1 += bzSlopeY1;
		  buSlopeYAccum0 += buSlopeY0;
		  buSlopeYAccum1 += buSlopeY1;
		  bvSlopeYAccum0 += bvSlopeY0;
		  bvSlopeYAccum1 += bvSlopeY1;
		  CY1 += t.FDX12;
		  CY2 += t.FDX23;
		  CY3 += t.FDX31;
		  col += width;
		}
	  }
	}
  }
#endif
}

void DrawTrianglesDeferred(unsigned int flags)
{
  using std::min;
  using std::max;

  unsigned int width = wc_colorbuffer.w;
  unsigned int height = wc_colorbuffer.h;
  unsigned int numTilesX = (width >> Q);
  unsigned int numTilesY = (height >> Q);
  /*
  wc_tileList.clear();
  wc_tileListFilled.clear();
  wc_tileList.resize(numTilesX * numTilesY);
  wc_tileListFilled.resize(numTilesX * numTilesY);
  */
  if(wc_tileList.size() < (numTilesX * numTilesY))
	wc_tileList.resize(numTilesX * numTilesY);
  if(wc_tileListFilled.size() < (numTilesX * numTilesY))
	wc_tileListFilled.resize(numTilesX * numTilesY);

  for(int i = 0; i < wc_tileList.size(); ++i){
	  wc_tileList[i].count = 0;
  }
  for(int i = 0; i < wc_tileListFilled.size(); ++i){
	  wc_tileListFilled[i].count = 0;
  }

  for(int i=0; i<wc_vertices->size(); i+=3){
    Vector4f& v1 = (*wc_vertices)[i+0];
    Vector4f& v2 = (*wc_vertices)[i+2];
    Vector4f& v3 = (*wc_vertices)[i+1];

    Vector4f& tc1 = (*wc_tcoords0)[i+0];
    Vector4f& tc2 = (*wc_tcoords0)[i+1];
    Vector4f& tc3 = (*wc_tcoords0)[i+2];

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

	// Block size, 16x16 must be POT
	const int Q = 4;
	const int q = (1<<Q);

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
	for(int y = miny; y < maxy; y += q){
	  for(int x = minx; x < maxx; x += q){

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

		int bwi0 = (((long long)Aw*bwx0)>>coeff_precision_base) + (((long long)Bw*bwy0)>>coeff_precision_base) + Cw; //top left
		int bwi1 = (((long long)Aw*bwx0)>>coeff_precision_base) + (((long long)Bw*bwy1)>>coeff_precision_base) + Cw; //bottom left
		int bwi2 = (((long long)Aw*bwx1)>>coeff_precision_base) + (((long long)Bw*bwy0)>>coeff_precision_base) + Cw; //top right		
		int bwi3 = (((long long)Aw*bwx1)>>coeff_precision_base) + (((long long)Bw*bwy1)>>coeff_precision_base) + Cw; //bottom right

		int bz0 = (((long long)Az*bzx0)>>depth_precision_base) + (((long long)Bz*bzy0)>>depth_precision_base) + Cz; //top left
		int bz1 = (((long long)Az*bzx0)>>depth_precision_base) + (((long long)Bz*bzy1)>>depth_precision_base) + Cz; //bottom left
		int bz2 = (((long long)Az*bzx1)>>depth_precision_base) + (((long long)Bz*bzy0)>>depth_precision_base) + Cz; //top right		
		int bz3 = (((long long)Az*bzx1)>>depth_precision_base) + (((long long)Bz*bzy1)>>depth_precision_base) + Cz; //bottom right

		//Compute u for the corners of the tile
		int bu0 = (((long long)Au*bwx0)>>coeff_precision_base) + (((long long)Bu*bwy0)>>coeff_precision_base) + Cu; //top left
		int bu1 = (((long long)Au*bwx0)>>coeff_precision_base) + (((long long)Bu*bwy1)>>coeff_precision_base) + Cu; //bottom left
		int bu2 = (((long long)Au*bwx1)>>coeff_precision_base) + (((long long)Bu*bwy0)>>coeff_precision_base) + Cu; //top right		
		int bu3 = (((long long)Au*bwx1)>>coeff_precision_base) + (((long long)Bu*bwy1)>>coeff_precision_base) + Cu; //bottom right

		//Compute v for the corners of the tile
		int bv0 = (((long long)Av*bwx0)>>coeff_precision_base) + (((long long)Bv*bwy0)>>coeff_precision_base) + Cv; //top left
		int bv1 = (((long long)Av*bwx0)>>coeff_precision_base) + (((long long)Bv*bwy1)>>coeff_precision_base) + Cv; //bottom left
		int bv2 = (((long long)Av*bwx1)>>coeff_precision_base) + (((long long)Bv*bwy0)>>coeff_precision_base) + Cv; //top right		
		int bv3 = (((long long)Av*bwx1)>>coeff_precision_base) + (((long long)Bv*bwy1)>>coeff_precision_base) + Cv; //bottom right

		int bw0, bw1, bw2, bw3;
		bw0 = bw1 = bw2 = bw3 = 0;
		
		if(bwi0) bw0 = ((int)1<<(coeff_precision_base * 2)) / bwi0;
		if(bwi1) bw1 = ((int)1<<(coeff_precision_base * 2)) / bwi1;
		if(bwi2) bw2 = ((int)1<<(coeff_precision_base * 2)) / bwi2;
		if(bwi3) bw3 = ((int)1<<(coeff_precision_base * 2)) / bwi3;			   

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
		tile.bw0 = bw0; tile.bw1 = bw1;
		tile.bw2 = bw2; tile.bw3 = bw3;
		tile.bz0 = bz0; tile.bz1 = bz1;
		tile.bz2 = bz2; tile.bz3 = bz3;
		tile.bu0 = bu0; tile.bu1 = bu1;
		tile.bu2 = bu2; tile.bu3 = bu3;
		tile.bv0 = bv0; tile.bv1 = bv1;
		tile.bv2 = bv2; tile.bv3 = bv3;
		tile.zMin = std::min(std::min(std::min(bz0, bz1), bz2), bz3);
		tile.zMax = std::max(std::max(std::max(bz0, bz1), bz2), bz3);

		int tileIdx = (x >> Q) + (y >> Q) * numTilesX;
		tTest = (a == 0xF && b == 0xF && c == 0xF);
		// Accept whole block when totally covered
		if(a == 0xF && b == 0xF && c == 0xF){
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
  BlitTilesFilled();
  BlitTiles();
}

bool ComputeCoeffMatrix(const Vector4f& v1, const Vector4f& v2, const Vector4f& v3, Matrix3f& m)
{
  //fesetexceptflag
  m = Matrix3f(Vector3f(v1.x, v1.y, v1.w),
			   Vector3f(v2.x, v2.y, v2.w),
			   Vector3f(v3.x, v3.y, v3.w));
  const float eps = (1.0f / 128.0f);
  float det = 
	m[0]*(m[4]*m[8] - m[5]*m[7]) +
	m[1]*(m[5]*m[6] - m[3]*m[8]) +
	m[2]*(m[3]*m[7] - m[4]*m[6]);

  //Degenerate or really small triangles
  //have really tiny determinants or equal
  //zero. We don't render these.
  if(std::abs(det) < 0.0125f){
	return false;
  }

  //Triangles with negative determinants
  //are backfaces, which don't need to be rendered
  if(det < 0.0f){
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
  m[0] = c00; m[1] = c01; m[2] = c02;
  m[3] = c10; m[4] = c11; m[5] = c12;
  m[6] = c20; m[7] = c21; m[8] = c22;

  det = 1.0f/det;
  m[0] *= det; m[1] *= det; m[2] *= det;
  m[3] *= det; m[4] *= det; m[5] *= det;
  m[6] *= det; m[7] *= det; m[8] *= det;

  return true;
}

inline static void SR_InterpTransform(float& f1, float& f2, float& f3, const Matrix3f& m)
{
  /* Vector made up by one scalar from each vertex */
  Vector3f v(f1, f2, f3);

  /* Multiply by coefficient matrix */
  v = m * v;
  /* Assign the transformed values back */

  f1 = v.x;
  f2 = v.y;
  f3 = v.z;
}

inline static void SR_InterpTransform(Vector4f& v1, Vector4f& v2, Vector4f& v3, const Matrix3f& m)
{
  SR_InterpTransform(v1.x, v2.x, v3.x, m);
  SR_InterpTransform(v1.y, v2.y, v3.y, m);
  SR_InterpTransform(v1.z, v2.z, v3.z, m);
  SR_InterpTransform(v1.w, v2.w, v3.w, m);
}

void SR_Render(unsigned int flags)
{
  //Matrix4f modelviewProjection = wc_projection * wc_modelview;
  size_t prevSize = wc_vertices->size();
  /*
  for(int i = 0; i < prevSize; i+=3){
    (*wc_vertices)[i+0] = modelviewProjection * (*wc_vertices)[i+0];
    (*wc_vertices)[i+1] = modelviewProjection * (*wc_vertices)[i+1];
    (*wc_vertices)[i+2] = modelviewProjection * (*wc_vertices)[i+2];
  }
  */
  size_t oldSize = wc_vertices->size();
  for(int i = 0; i < oldSize; i+=3){
	/* Compute [a,b,c] coefficients */
	Matrix3f m;
	bool b = ComputeCoeffMatrix((*wc_vertices)[i+0], (*wc_vertices)[i+1], (*wc_vertices)[i+2], m);
	//Skip degenerate triangles, small triangles and backfaces
	if(!b) continue;

	//Project() :
	//Compute screen space coordinates for x and y
	//Normalize z into [0.0f, 1.0f> half-range, Q0.16 fixedpoint
	(*wc_vertices)[i+0] = project((*wc_vertices)[i+0], wc_colorbuffer.w, wc_colorbuffer.h);
	(*wc_vertices)[i+1] = project((*wc_vertices)[i+1], wc_colorbuffer.w, wc_colorbuffer.h);
	(*wc_vertices)[i+2] = project((*wc_vertices)[i+2], wc_colorbuffer.w, wc_colorbuffer.h);
	
	//Must interpolate z linearly in screenspace!
	//To get the coefficients required for an affine interpolation, simply multiply z with w
	(*wc_vertices)[i+0].z *= (*wc_vertices)[i+0].w;
	(*wc_vertices)[i+1].z *= (*wc_vertices)[i+1].w;
	(*wc_vertices)[i+2].z *= (*wc_vertices)[i+2].w;
	//SR_InterpTransform(wc_vertices[i+0].z, wc_vertices[i+1].z, wc_vertices[i+2].z, m);

	Vector3f zv((*wc_vertices)[i+0].z, (*wc_vertices)[i+1].z, (*wc_vertices)[i+2].z);
	zv = m * zv;
	(*wc_vertices)[i+0].z = zv.x;
	(*wc_vertices)[i+1].z = zv.y;
	(*wc_vertices)[i+2].z = zv.z;

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

  switch(flags){
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
  case 0:
	{

	}
  }
  //switch not done yet, so just call the test rasterizer
  //DrawTriangles(flags);
  DrawTrianglesDeferred(flags);
}
