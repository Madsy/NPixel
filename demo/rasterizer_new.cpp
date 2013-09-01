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

inline int iround(float f)
{
  int r;

  asm volatile ("fistpl %[output]\n"
	   : [output] "=m" (r)
	   : [input] "t" (f)
	   : "st(0)");

  return r;
}

inline int fpceil15(int fp)
{
  return (fp & 32767) ? ((fp & ~32767) + 32768) : fp;
}

inline int fpceil10(int fp)
{
  return (fp & 1023) ? ((fp & ~1023) + 1024) : fp;
}

//int zmin, zmax;

template<class T> void DrawTriangles(T fs)
{
  using std::min;
  using std::max;

  //zmin = 65536;
  //zmax = -65536;

  for(int i=0; i<wc_vertices.size(); i+=3){
    Vector4f& v1 = wc_vertices[i+0];
    Vector4f& v2 = wc_vertices[i+2];
    Vector4f& v3 = wc_vertices[i+1];

    Vector4f& tc1 = wc_tcoords0[i+0];
    Vector4f& tc2 = wc_tcoords0[i+1];
    Vector4f& tc3 = wc_tcoords0[i+2];

	Coeff wCoeff(wc_vertices[i+0].w, wc_vertices[i+1].w, wc_vertices[i+2].w);
	Coeff texCoordCoeffU(tc1.x, tc2.x, tc3.x);
	Coeff texCoordCoeffV(tc1.y, tc2.y, tc3.y);

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

	// Block size, standard 8x8 (must be power of two)
	const int Q = 3;
	const int q = (1<<Q);

	// Start in corner of 8x8 block
	minx &= ~(q - 1);
	miny &= ~(q - 1);

	unsigned int width = wc_colorbuffer.w;
	unsigned int height = wc_colorbuffer.h;

	//buffer += (miny*width);
	int	colOffs = miny;
	int colTileOffs = 0;
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
		bool boundTest1 = true;
		
		bool px0min = x0 > 0;
		bool px0max = x0 < width;
		bool py0min = y0 > 0;
		bool py0max = y0 < height;

		bool px1min = x1 > 0;
		bool px1max = x1 < width;
		bool py1min = y1 > 0;
		bool py1max = y1 < height;

		int pflags0 = (px0min << 3) | (px1min << 2) | (px0max << 1) | px1max;
		int pflags1 = (py0min << 3) | (py1min << 2) | (py0max << 1) | py1max;
		
		if(pflags0 == 0xF && pflags1 == 0xF){
		  //Completely inside the frustum
		  boundTest1 = false;
		} else if(pflags0 == 0x3 || pflags0 == 0xC ||
				  pflags1 == 0x3 || pflags1 == 0xC){
		  // pflags0 == 0x3: outside right X plane
		  // pflags0 == 0xC: outside left X plane
		  // pflags1 == 0x3: outside bottom Y plane
		  // pflags1 == 0xC: outside top Y plane
		  continue;
		} else {
		  boundTest1 = true;
		}
		
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

		const int coeff_precision_base = 11;
		const int ndc_precision_base = 20;
		const int depth_precision_base = 16;
		const float f_coeff_precision = (float)(1 << coeff_precision_base);
		const float f_ndc_precision = (float)(1 << ndc_precision_base);
		const float f_depth_precision = (float)(1 << depth_precision_base);
		const int i_coeff_precision = 1 << coeff_precision_base;
		const int i_ndc_precision = 1 << ndc_precision_base;
		const int i_depth_precision = 1 << depth_precision_base;
		const int base_diff = (ndc_precision_base - coeff_precision_base);
		const int base_diff_z = (ndc_precision_base - depth_precision_base);
		/* s = Ax + By + C */
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

		// screenspace -> NDC space
		int NDC_x_step = fHalfWidthInv * f_ndc_precision;
		int NDC_y_step = fHalfHeightInv * f_ndc_precision;
		int NDC_x0 = x * NDC_x_step;  //halfWidthInv * f_ndc_precision;
		int NDC_y0 = y * NDC_y_step; //halfHeighthInv * f_ndc_precision;
		int NDC_x1 = (x + q - 1) * NDC_x_step;
		int NDC_y1 = (y + q - 1) * NDC_y_step;

		//compute 1/w for corner points here
		int bwx0 = (NDC_x0 - i_ndc_precision) >> base_diff;
		int bwx1 = (NDC_x1 - i_ndc_precision) >> base_diff;
		int bwy0 = (NDC_y0 - i_ndc_precision) >> base_diff;
		int bwy1 = (NDC_y1 - i_ndc_precision) >> base_diff;

		int bwi0 = (((long long)Aw*bwx0)>>coeff_precision_base) + (((long long)Bw*bwy0)>>coeff_precision_base) + Cw; //left
		int bwi1 = (((long long)Aw*bwx0)>>coeff_precision_base) + (((long long)Bw*bwy1)>>coeff_precision_base) + Cw; //left
		int bwi2 = (((long long)Aw*bwx1)>>coeff_precision_base) + (((long long)Bw*bwy0)>>coeff_precision_base) + Cw; //right		
		int bwi3 = (((long long)Aw*bwx1)>>coeff_precision_base) + (((long long)Bw*bwy1)>>coeff_precision_base) + Cw; //right

		int bw0, bw1, bw2, bw3;
		bw0 = bw1 = bw2 = bw3 = 0;

		if(bwi0)
		  bw0 = ((int)1<<(coeff_precision_base * 2)) / bwi0;
		if(bwi1)
		  bw1 = ((int)1<<(coeff_precision_base * 2)) / bwi1;
		if(bwi2)
		  bw2 = ((int)1<<(coeff_precision_base * 2)) / bwi2;
		if(bwi3)
		  bw3 = ((int)1<<(coeff_precision_base * 2)) / bwi3;

		int bwSlopeY0 = bw1 - bw0; //vertical slope left
		int bwSlopeY1 = bw3 - bw2; //vertical slope right

		// Accept whole block when totally covered
		if(a == 0xF && b == 0xF && c == 0xF){
		  unsigned int col = y*width;
		  const unsigned int* tbuf = &wc_texture0->texels[0];
		  int iTw = wc_texture0->width;
		  int iTh = wc_texture0->height;
		  int NDC_iy = NDC_y0; //current y, or iy in NDC space

		  int bwSlopeYAccum0 = bw0 << Q; //start vertical w left
		  int bwSlopeYAccum1 = bw2 << Q; //start vertical w right

		  for(int iy = y; iy < y + q; iy++){
			/* TODO: Ensure that -iy can never be larger than the tilesize q */
			if(boundTest1) {
			  if(iy < 0){
				int skip = -iy;
				NDC_iy += (NDC_y_step * skip);
			    bwSlopeYAccum0 += (bwSlopeY0 * skip);
				bwSlopeYAccum1 += (bwSlopeY1 * skip);
				col += (width * skip);
				continue;
			  } else if(iy >= height){
				break;
			  }
			}
			int NDC_ix = NDC_x0;
			int bwSlopeX0 = bwSlopeYAccum1 - bwSlopeYAccum0;
			int bwSlopeXAccum0 = bwSlopeYAccum0 << Q; //start horisontal w left
			for(int ix = x; ix < x + q; ix++){
			  if(boundTest1){
				if(ix < 0){
				  int skip = -ix;
				  NDC_ix += (NDC_x_step * skip);
				  bwSlopeXAccum0 += (bwSlopeX0 * skip);
				  continue;
				} else if(ix >= width){
				  break;
				}
			  }

			  int interpX  = (NDC_ix - i_ndc_precision) >> base_diff;
			  int interpY  = (NDC_iy - i_ndc_precision) >> base_diff;
			  int interpZX = (NDC_ix - i_ndc_precision) >> base_diff_z;
			  int interpZY = (NDC_iy - i_ndc_precision) >> base_diff_z;

			  //maybe need to cast later
			  unsigned short z =
				((((long long)Az*interpZX) + ((long long)Bz*interpZY))>>depth_precision_base) + Cz;


			  if(z < wc_depthbuffer[ix + col]){
				wc_depthbuffer[ix + col] = z;
				//int wi = (((long long)Aw*interpX)>>coeff_precision_base) + (((long long)Bw*interpY)>>coeff_precision_base) + Cw;
				int uw = (((long long)Au*interpX)>>coeff_precision_base) + (((long long)Bu*interpY)>>coeff_precision_base) + Cu;
				int vw = (((long long)Av*interpX)>>coeff_precision_base) + (((long long)Bv*interpY)>>coeff_precision_base) + Cv;
				//int w = ((int)1<<(coeff_precision_base * 2)) / wi;
				int w = bwSlopeXAccum0>>(Q*2);
				int u = ((long long)uw*w*iTw) >> (coeff_precision_base * 2);
				int v = ((long long)vw*w*iTh) >> (coeff_precision_base * 2);
				u = clamp((int)u, 0, iTw-1);
				v = clamp((int)v, 0, iTh-1);
				int idxTex = u + v*iTw;
				wc_colorbuffer[ix + col] = tbuf[idxTex];
			  }
			  NDC_ix += NDC_x_step;
			  bwSlopeXAccum0 += bwSlopeX0;
			}
			bwSlopeYAccum0 += bwSlopeY0;
			bwSlopeYAccum1 += bwSlopeY1;
			NDC_iy += NDC_y_step;
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
		  int NDC_iy = NDC_y0; //current y, or iy in NDC space
		  for(int iy = y; iy < y + q; iy++){
			if(boundTest1) {
			  if(iy < 0){
				int skip = -iy;
				NDC_iy += (NDC_y_step * skip);
				col += (width * skip);
				CY1 += (FDX12 * skip);
				CY2 += (FDX23 * skip);
				CY3 += (FDX31 * skip);
				continue;
			  } else if(iy >= height){
				break;
			  }
			}
			int CX1 = CY1;
			int CX2 = CY2;
			int CX3 = CY3;
			int NDC_ix = NDC_x0;
			for(int ix = x; ix < x + q; ix++){
			  if(boundTest1){
				if(ix < 0){
				  int skip = -ix;
				  NDC_ix += (NDC_x_step * skip);
				  CX1 -= (FDY12 * skip);
				  CX2 -= (FDY23 * skip);
				  CX3 -= (FDY31 * skip);
				  continue;
				} else if(ix >= width){
				  break;
				}
			  }
			  if(CX1 > 0 && CX2 > 0 && CX3 > 0){

				const int base_diff = (ndc_precision_base - coeff_precision_base);
				const int base_diff_z = (ndc_precision_base - depth_precision_base);

				int interpX = (NDC_ix - i_ndc_precision) >> base_diff;
				int interpY = (NDC_iy - i_ndc_precision) >> base_diff;
				int interpZX = (NDC_ix - i_ndc_precision) >> base_diff_z;
				int interpZY = (NDC_iy - i_ndc_precision) >> base_diff_z;

				//maybe need to cast later
				unsigned short z =
				(((Az*interpZX) + (Bz*interpZY))>>depth_precision_base) + Cz;

				if(z < wc_depthbuffer[ix + col]){
				  wc_depthbuffer[ix + col] = z;

				  int wi = (((long long)Aw*interpX)>>coeff_precision_base) + (((long long)Bw*interpY)>>coeff_precision_base) + Cw;
				  int uw = (((long long)Au*interpX)>>coeff_precision_base) + (((long long)Bu*interpY)>>coeff_precision_base) + Cu;
				  int vw = (((long long)Av*interpX)>>coeff_precision_base) + (((long long)Bv*interpY)>>coeff_precision_base) + Cv;
				  int w = ((int)1<<(coeff_precision_base * 2)) / wi;
				  int u = ((long long)uw*w*iTw) >> (coeff_precision_base * 2);
				  int v = ((long long)vw*w*iTh) >> (coeff_precision_base * 2);
				  u = clamp((int)u, 0, iTw-1);
				  v = clamp((int)v, 0, iTh-1);
				  int idxTex = u + v*iTw;
				  wc_colorbuffer[ix + col] = tbuf[idxTex];
				}				
			  }
			  NDC_ix += NDC_x_step;
			  CX1 -= FDY12;
			  CX2 -= FDY23;
			  CX3 -= FDY31;
			}
			NDC_iy += NDC_y_step;
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

  if(std::abs(det) < 0.0125f){
	return false;
  }

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
  //float width = wc_colorbuffer.w;
  //float height = wc_colorbuffer.h;

  //f1 = v.x*(2.0f / width);
  //f2 = v.y*(2.0f / height);
  //f3 = v.z - v.x - v.y;

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
  Matrix4f modelviewProjection = wc_projection * wc_modelview;
  size_t prevSize = wc_vertices.size();
  for(int i = 0; i < prevSize; i+=3){
	wc_vertices[i+0] = modelviewProjection * wc_vertices[i+0];
	wc_vertices[i+1] = modelviewProjection * wc_vertices[i+1];
	wc_vertices[i+2] = modelviewProjection * wc_vertices[i+2];
  }

  size_t oldSize = wc_vertices.size();
  for(int i = 0; i < oldSize; i+=3){
	/* Compute [a,b,c] coefficients */
	Matrix3f m;
	bool b = ComputeCoeffMatrix(wc_vertices[i+0], wc_vertices[i+1], wc_vertices[i+2], m);
	if(!b) continue;

	//Project() :
	//Compute screen space coordinates for x and y
	//Normalize z into [0.0f, 1.0f> half-range
	wc_vertices[i+0] = project(wc_vertices[i+0], wc_colorbuffer.w, wc_colorbuffer.h);
	wc_vertices[i+1] = project(wc_vertices[i+1], wc_colorbuffer.w, wc_colorbuffer.h);
	wc_vertices[i+2] = project(wc_vertices[i+2], wc_colorbuffer.w, wc_colorbuffer.h);
	
	//Must interpolate z linearly in screenspace!
	//To get the coefficients required for an affine interpolation, simply multiply z with w
	wc_vertices[i+0].z *= wc_vertices[i+0].w;
	wc_vertices[i+1].z *= wc_vertices[i+1].w;
	wc_vertices[i+2].z *= wc_vertices[i+2].w;
	//SR_InterpTransform(wc_vertices[i+0].z, wc_vertices[i+1].z, wc_vertices[i+2].z, m);

	Vector3f zv(wc_vertices[i+0].z, wc_vertices[i+1].z, wc_vertices[i+2].z);
	zv = m * zv;
	wc_vertices[i+0].z = zv.x;
	wc_vertices[i+1].z = zv.y;
	wc_vertices[i+2].z = zv.z;

	// To get "1.0f / w", multiply the 3D Vector [1,1,1] with the coefficient matrix.
	// We don't need w anymore after this point. It's stored as a part of the matrix.
	wc_vertices[i+0].w = wc_vertices[i+1].w = wc_vertices[i+2].w = 1.0f;
	SR_InterpTransform(wc_vertices[i+0].w, wc_vertices[i+1].w, wc_vertices[i+2].w, m);

	/* Compute the coefficients for the rest of the per-vertex data */
	if(flags & SR_TEXCOORD0)
	  SR_InterpTransform(wc_tcoords0[i+0], wc_tcoords0[i+1], wc_tcoords0[i+2], m);
	if(flags & SR_TEXCOORD1)
	  SR_InterpTransform(wc_tcoords1[i+0], wc_tcoords1[i+1], wc_tcoords1[i+2], m);
	if(flags & SR_LIGHTING)
	  SR_InterpTransform(wc_normals[i+0], wc_normals[i+1], wc_normals[i+2], m);
	if(flags & SR_COLOR)
	  SR_InterpTransform(wc_colors[i+0], wc_colors[i+1], wc_colors[i+2], m);

	wc_vertices.push_back(wc_vertices[i+0]);
	wc_vertices.push_back(wc_vertices[i+1]);
	wc_vertices.push_back(wc_vertices[i+2]);
	wc_tcoords0.push_back(wc_tcoords0[i+0]);
	wc_tcoords0.push_back(wc_tcoords0[i+1]);
	wc_tcoords0.push_back(wc_tcoords0[i+2]);
  }
  wc_vertices.erase(wc_vertices.begin(), wc_vertices.begin() + oldSize);
  wc_tcoords0.erase(wc_tcoords0.begin(), wc_tcoords0.begin() + oldSize);

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
  DrawTriangles(fs_basic);
}
