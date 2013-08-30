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

inline int fpceil15(int fp)
{
  return (fp & 32767) ? ((fp & ~32767) + 32768) : fp;
}

inline int fpceil10(int fp)
{
  return (fp & 1023) ? ((fp & ~1023) + 1024) : fp;
}


  /*
	wc_colorbuffer
	wc_depthbuffer
	wc_vertices
	wc_tcoords0
	wc_tcoords1
	wc_normals
	wc_colors
	wc_modelview;
	wc_projection;
  */


template<class T> void DrawTriangles(T fs)
{
  using std::min;
  using std::max;

  for(int i=0; i<wc_vertices.size(); i+=3){
    Vector4f& v1 = wc_vertices[i+0];
    Vector4f& v2 = wc_vertices[i+2];
    Vector4f& v3 = wc_vertices[i+1];

    //Vector4f& tc_0_1 = wc_tcoords0[i+0];
    //Vector4f& tc_0_2 = wc_tcoords0[i+2];
    //Vector4f& tc_0_3 = wc_tcoords0[i+1];

	Coeff wCoeff(wc_vertices[i+0].w, wc_vertices[i+1].w, wc_vertices[i+2].w);
	Coeff texCoordCoeff0U(wc_tcoords0[i+0].x, wc_tcoords0[i+1].x, wc_tcoords0[i+2].x);
	Coeff texCoordCoeff0V(wc_tcoords0[i+0].y, wc_tcoords0[i+1].y, wc_tcoords0[i+2].y);
	//Coeff texCoordCoeff1U(wc_tcoords1[i+0].x, wc_tcoords1[i+1].x, wc_tcoords1[i+2].x);
	//Coeff texCoordCoeff1V(wc_tcoords1[i+0].y, wc_tcoords1[i+1].y, wc_tcoords1[i+2].y);


	// 28.4 fixed-point coordinates
	const int Y1 = (int)(16.0f * v1.y);
	const int Y2 = (int)(16.0f * v2.y);
	const int Y3 = (int)(16.0f * v3.y);
	const int X1 = (int)(16.0f * v1.x);
	const int X2 = (int)(16.0f * v2.x);
	const int X3 = (int)(16.0f * v3.x);

	bool doOnce = true; //debug
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
	const int q = 8;

	// Start in corner of 8x8 block
	minx &= ~(q - 1);
	miny &= ~(q - 1);

	const unsigned int stride = wc_colorbuffer.w;
	//buffer += (miny*stride);
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
		int x0 = x << 4;
		int x1 = (x + q - 1) << 4;
		int y0 = y << 4;
		int y1 = (y + q - 1) << 4;

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

		// Accept whole block when totally covered
		if(a == 0xF && b == 0xF && c == 0xF){
		  for(int iy = y; iy < y + q; iy++){
			for(int ix = x; ix < x + q; ix++){
			  //wc_colorbuffer[ix + iy*stride] = 0xFFFF00FF;
			  fs(ix, iy, &wc_colorbuffer[ix + iy*stride], wCoeff, texCoordCoeff0U, texCoordCoeff0V);
			}
		  }
		} else { // Partially covered
		  int CY1 = C1 + DX12 * y0 - DY12 * x0;
		  int CY2 = C2 + DX23 * y0 - DY23 * x0;
		  int CY3 = C3 + DX31 * y0 - DY31 * x0;
		  for(int iy = y; iy < y + q; iy++){
			int CX1 = CY1;
			int CX2 = CY2;
			int CX3 = CY3;
			for(int ix = x; ix < x + q; ix++){
			  if(CX1 > 0 && CX2 > 0 && CX3 > 0){
				//wc_colorbuffer[ix + iy*stride] = 0xFFFFFFFF;
				fs(ix, iy, &wc_colorbuffer[ix + iy*stride], wCoeff, texCoordCoeff0U, texCoordCoeff0V);
			  }
			  CX1 -= FDY12;
			  CX2 -= FDY23;
			  CX3 -= FDY31;
			}
			CY1 += FDX12;
			CY2 += FDX23;
			CY3 += FDX31;
		  }
		}
	  }
	}
  }
}

Matrix3f ComputeCoeffMatrix(const Vector4f& v1, const Vector4f& v2, const Vector4f& v3)
{
  //fesetexceptflag
/*
  Matrix3f m(Vector3f(v1.x, v2.x, v3.x),
			 Vector3f(v1.y, v2.y, v3.y),
			 Vector3f(v1.w, v2.w, v3.w));
*/
  Matrix3f m(Vector3f(v1.x, v1.y, v1.w),
			 Vector3f(v2.x, v2.y, v2.w),
			 Vector3f(v3.x, v3.y, v3.w));
  const float eps = 0.0001f;
  float det = 
	m[0]*(m[4]*m[8] - m[5]*m[7]) +
	m[1]*(m[5]*m[6] - m[3]*m[8]) +
	m[2]*(m[3]*m[7] - m[4]*m[6]);

  //ASSERT(std::abs(det) > eps);
  if(std::abs(det) < eps)
	return Matrix3f();

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

  return m;
}

inline static void SR_InterpTransform(float& f1, float& f2, float& f3, const Matrix3f& m)
{
  /* Vector made up by one scalar from each vertex */
  Vector3f v(f1, f2, f3);

  /* Multiply by coefficient matrix */
  v = m * v;
  /* Assign the transformed values back */
  f1 = v.x*(2.0f / 640.0f);
  f2 = v.y*(2.0f / 360.0f);
  f3 = v.z - v.x - v.y;
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
  for(int i = 0; i < wc_vertices.size(); i+=3){
	wc_vertices[i+0] = modelviewProjection * wc_vertices[i+0];
	wc_vertices[i+1] = modelviewProjection * wc_vertices[i+1];
	wc_vertices[i+2] = modelviewProjection * wc_vertices[i+2];
  }

  clip_triangle(Vector4f(-1.0f,  0.0f, 0.0f, 1.0f), flags);
  clip_triangle(Vector4f( 1.0f,  0.0f, 0.0f, 1.0f), flags);
  clip_triangle(Vector4f( 0.0f,  1.0f, 0.0f, 1.0f), flags);
  clip_triangle(Vector4f( 0.0f, -1.0f, 0.0f, 1.0f), flags);
  clip_triangle(Vector4f( 0.0f,  0.0f,-1.0f, 1.0f), flags);
  clip_triangle(Vector4f( 0.0f,  0.0f, 1.0f, 1.0f), flags);

  for(int i = 0; i < wc_vertices.size(); i+=3){
	/* Compute [a,b,c] coefficients */
	Matrix3f m = ComputeCoeffMatrix(wc_vertices[i+0], wc_vertices[i+1], wc_vertices[i+2]);

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
	SR_InterpTransform(wc_vertices[i+0].z, wc_vertices[i+1].z, wc_vertices[i+2].z, m);

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
  }

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
