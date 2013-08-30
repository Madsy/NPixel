#include <vector>
#include <cstdio>
#include <SDL/SDL.h>
#include <linealg.h>
#include <fixedpoint.h>
#include "rasterizer.h"
#include "framebuffer.h"
#include "texture.h"
#include "myassert.h"

inline int fpceil15(int fp)
{
  return (fp & 32767) ? ((fp & ~32767) + 32768) : fp;
}

inline int fpceil10(int fp)
{
  return (fp & 1023) ? ((fp & ~1023) + 1024) : fp;
}


static void drawScanLine(unsigned int* cbuffer,
		  int width,
		  int height,
		  int y,
		  int x1, int x2,
		  int z1, int z2,
		  int w1, int w2,
		  int s1, int s2,
		  int t1, int t2)
{
  unsigned short* zbuffer;
  unsigned short z;
  int deltaX, deltaZ, deltaW, deltaS, deltaT;
  int slopeZ, slopeW, slopeS, slopeT;
  int zStart, zEnd, wStart, wEnd, sStart, sEnd, tStart, tEnd;
  int xError;
  /* new */
  int wCounter, wInterp0, wInterp1, wInterpDelta, wInterpSlope;
  int w, s, t;
  int texWidth, texHeight;
  int xStart, xEnd;
  int col;

  col = y*width;

  if(x1 > x2){
    std::swap(x1, x2);
    std::swap(z1, z2);
    std::swap(w1, w2);
    std::swap(s1, s2);
    std::swap(t1, t2);
  }

  xStart = fpceil15(x1);
  xEnd = fpceil15(x2) - 32768;
  xError = xStart - x1;
  xStart >>= 15;
  xEnd >>= 15;

  deltaX = x2 - x1;
  deltaZ = z2 - z1;
  deltaW = w2 - w1;
  deltaS = s2 - s1;
  deltaT = t2 - t1;

  slopeZ = slopeW = slopeS = slopeT = 0;

  if(deltaX != 0){
    deltaX = reci15(deltaX);
    slopeZ = ((long long)deltaZ * deltaX)>>15;
    slopeW = ((long long)deltaW * deltaX)>>15;
    slopeS = ((long long)deltaS * deltaX)>>15;
    slopeT = ((long long)deltaT * deltaX)>>15;
  } else {
    return;
  }
  
  /* start interpolants */
  sStart = s1;
  tStart = t1;
  zStart = z1;
  wStart = w1;

  /* Correct for the new x position */
  zStart += ((long long)slopeZ * xError)>>15;
  wStart += ((long long)slopeW * xError)>>15;
  sStart += ((long long)slopeS * xError)>>15;
  tStart += ((long long)slopeT * xError)>>15;  

  zbuffer = &depthbuffer.data[0];
  const unsigned int* texture = &currentTexture->color[0];
  texWidth = currentTexture->width;
  texHeight = currentTexture->height;
  
  int indexDst = xStart + col;    
  int indexSrc;
  /* obs, just rearranging [xStart xEnd] to start from 0 here */

  xEnd = xEnd-xStart;
  xStart = 0;
  int xRemain = xEnd;
  for(; xStart <= xEnd; ++xStart){
    z = zStart;
    /* update interpolants whenever xPos mod 16 is 0 */
    const int OPTIMIZE_ROW_INTERP_LEN = 8;
#if 0
    if(xRemain >= OPTIMIZE_ROW_INTERP_LEN){
      if(!(xStart & (OPTIMIZE_ROW_INTERP_LEN - 1))){
	wInterp0 = wStart + slopeW;
	wInterp1 = wStart + slopeW * OPTIMIZE_ROW_INTERP_LEN;
	if(wInterp0 != 0) wInterp0 = reci15(wInterp0);
	else wInterp0 = 0;
	if(wInterp1 != 0) wInterp1 = reci15(wInterp1);
	else wInterp1 = 0;
	wInterpDelta = wInterp1-wInterp0;
	wInterpSlope = wInterpDelta / OPTIMIZE_ROW_INTERP_LEN;
	w = wInterp0;
      }
    } else {
      if(wStart != 0) w = reci15(wStart);
      else w = 0;
    }
    --xRemain;
#endif
    w = wStart;

    if(z < zbuffer[indexDst]){
      zbuffer[indexDst] = z;
      if(w) w = reci15(w);
      s = ((long long)w * sStart)>>15;
      t = ((long long)w * tStart)>>15;
      s = ((long long)s * (texWidth - 1));
      t = ((long long)t * (texHeight - 1));
      s >>= 15;
      t >>= 15;
 
      if(s < 0) s = 0;
      if(s >= texWidth) s = texWidth-1;
      if(t < 0) t = 0;
      if(t >= texHeight) t = texHeight-1;
 
      ASSERT(s >= 0 && s < texWidth);
      ASSERT(t >= 0 && t < texHeight);
      
      indexSrc = s + t*texWidth;
      cbuffer[indexDst] = texture[indexSrc];      
    }
    ++indexDst;
    zStart += slopeZ;
    wStart += slopeW;
    /* lerp w over two correct samples instead of
       lerping 1/w */
    w += wInterpSlope;
    sStart += slopeS;
    tStart += slopeT;
  }
}

void DrawTriangle(std::vector<Vector4f>& vertexData,
		  std::vector<Vector4f>& textureData,
		  unsigned int* buffer,
		  unsigned int width,
		  unsigned int height
		  )
{
  //  const float eps = 0.00001f;
  for(int i=0; i<vertexData.size(); i+=3){
    Vector4f& v1 = vertexData[i+0];
    Vector4f& v2 = vertexData[i+1];
    Vector4f& v3 = vertexData[i+2];

    Vector4f& tc1 = textureData[i+0];
    Vector4f& tc2 = textureData[i+1];
    Vector4f& tc3 = textureData[i+2];

    /* deltas below are always positive due to this sorting.
       v1 = top, v2 = middle, v3 = bottom */
    if(v1.y > v2.y){
      std::swap(v1, v2);
      std::swap(tc1, tc2);
    }
    if(v2.y > v3.y){
      std::swap(v2, v3);
      std::swap(tc2, tc3);
    }
    if(v1.y > v2.y){
      std::swap(v1, v2);
      std::swap(tc1, tc2);
    }

    /* Q15.16 fixedpoint values*/
    Vector4i v1fp(v1.x * 32768.0f, v1.y * 32768.0f, v1.z * 65535.0f, v1.w * 32768.0f);
    Vector4i v2fp(v2.x * 32768.0f, v2.y * 32768.0f, v2.z * 65535.0f, v2.w * 32768.0f);
    Vector4i v3fp(v3.x * 32768.0f, v3.y * 32768.0f, v3.z * 65535.0f, v3.w * 32768.0f);

    Vector4i tc1fp(tc1.x * 32768.0f, tc1.y * 32768.0f, 0.0f, 0.0f);
    Vector4i tc2fp(tc2.x * 32768.0f, tc2.y * 32768.0f, 0.0f, 0.0f);
    Vector4i tc3fp(tc3.x * 32768.0f, tc3.y * 32768.0f, 0.0f, 0.0f);

    Vector4i delta1PTfp = v2fp - v1fp;
    Vector4i delta2PTfp = v3fp - v1fp;
    Vector4i delta3PTfp = v3fp - v2fp;

    Vector4i delta1TCfp = tc2fp - tc1fp;
    Vector4i delta2TCfp = tc3fp - tc1fp;
    Vector4i delta3TCfp = tc3fp - tc2fp;

    /* w not copied properly in operator- for Vector4<T>*/
    delta1PTfp.w = v2fp.w - v1fp.w;
    delta2PTfp.w = v3fp.w - v1fp.w;
    delta3PTfp.w = v3fp.w - v2fp.w;

    /* Slopes for x,z,w (vertices) and u,v (tcoords) */
    int slope1X, slope2X, slope3X;
    int slope1Z, slope2Z, slope3Z;
    int slope1W, slope2W, slope3W;
    int slope1S, slope2S, slope3S;
    int slope1T, slope2T, slope3T;

    slope1X = slope2X = slope3X = 0;
    slope1Z = slope2Z = slope3Z = 0;
    slope1W = slope2W = slope3W = 0;
    slope1S = slope2S = slope3S = 0;
    slope1T = slope2T = slope3T = 0;

    if(delta1PTfp.y != 0){
      slope1X = ((long long)delta1PTfp.x << 15) / delta1PTfp.y;
      slope1Z = ((long long)delta1PTfp.z << 15) / delta1PTfp.y;
      slope1W = ((long long)delta1PTfp.w << 15) / delta1PTfp.y;
      slope1S = ((long long)delta1TCfp.x << 15) / delta1PTfp.y;
      slope1T = ((long long)delta1TCfp.y << 15) / delta1PTfp.y;
    }

    if(delta2PTfp.y != 0){
      slope2X = ((long long)delta2PTfp.x << 15) / delta2PTfp.y;
      slope2Z = ((long long)delta2PTfp.z << 15) / delta2PTfp.y;
      slope2W = ((long long)delta2PTfp.w << 15) / delta2PTfp.y;
      slope2S = ((long long)delta2TCfp.x << 15) / delta2PTfp.y;
      slope2T = ((long long)delta2TCfp.y << 15) / delta2PTfp.y;
    }

    if(delta3PTfp.y != 0){
      slope3X = ((long long)delta3PTfp.x << 15) / delta3PTfp.y;
      slope3Z = ((long long)delta3PTfp.z << 15) / delta3PTfp.y;
      slope3W = ((long long)delta3PTfp.w << 15) / delta3PTfp.y;
      slope3S = ((long long)delta3TCfp.x << 15) / delta3PTfp.y;
      slope3T = ((long long)delta3TCfp.y << 15) / delta3PTfp.y;
    }

    int y1, y2;
    int x1, x2;
    int z1, z2;
    int w1, w2;
    int s1, s2;
    int t1, t2;
    int yError, yErrorInv;

    y1 = fpceil15(v1fp.y);
    y2 = fpceil15(v2fp.y) - 32768;
    yError = y1 - v1fp.y;
    x1 = x2 = v1fp.x;
    z1 = z2 = v1fp.z;
    w1 = w2 = v1fp.w;
    s1 = s2 = tc1fp.x;
    t1 = t2 = tc1fp.y;

    /* Correct for the new y position */    
    x1 += ((long long)slope1X * yError) >> 15;
    x2 += ((long long)slope2X * yError) >> 15;
    z1 += ((long long)slope1Z * yError) >> 15;
    z2 += ((long long)slope2Z * yError) >> 15;
    w1 += ((long long)slope1W * yError) >> 15;
    w2 += ((long long)slope2W * yError) >> 15;
    s1 += ((long long)slope1S * yError) >> 15;
    s2 += ((long long)slope2S * yError) >> 15;
    t1 += ((long long)slope1T * yError) >> 15;
    t2 += ((long long)slope2T * yError) >> 15;

    y1 >>= 15;
    y2 >>= 15;
    /* Skipped if delta1f.y < 1 */
    for(; y1<=y2; ++y1){
      drawScanLine(buffer, width, height, y1, x1, x2, z1, z2, w1, w2, s1, s2, t1, t2);
      z1 += slope1Z;
      z2 += slope2Z;
      w1 += slope1W;
      w2 += slope2W;
      s1 += slope1S;
      s2 += slope2S;
      t1 += slope1T;
      t2 += slope2T;
      x1 += slope1X; /* middle - top */
      x2 += slope2X; /* bottom - top */
    }

    /* Next triangle part */
    y1 = fpceil15(v2fp.y);
    y2 = fpceil15(v3fp.y) - 32768;
    yError = y1 - v2fp.y;
    x1 = v2fp.x;
    z1 = v2fp.z;
    w1 = v2fp.w;
    s1 = tc2fp.x;
    t1 = tc2fp.y;

    /* Interpolate to find this point*/
    x2 = v1fp.x;
    z2 = v1fp.z;
    w2 = v1fp.w;
    s2 = tc1fp.x;
    t2 = tc1fp.y;
    x2 += ((long long)slope2X * delta1PTfp.y) >> 15;
    z2 += ((long long)slope2Z * delta1PTfp.y) >> 15;
    w2 += ((long long)slope2W * delta1PTfp.y) >> 15;
    s2 += ((long long)slope2S * delta1PTfp.y) >> 15;
    t2 += ((long long)slope2T * delta1PTfp.y) >> 15;

    /* Correct for the new y position */
    x1 += ((long long)slope3X * yError) >> 15;
    x2 += ((long long)slope2X * yError) >> 15;
    z1 += ((long long)slope3Z * yError) >> 15;
    z2 += ((long long)slope2Z * yError) >> 15;
    w1 += ((long long)slope3W * yError) >> 15;
    w2 += ((long long)slope2W * yError) >> 15;
    s1 += ((long long)slope3S * yError) >> 15;
    s2 += ((long long)slope2S * yError) >> 15;
    t1 += ((long long)slope3T * yError) >> 15;
    t2 += ((long long)slope2T * yError) >> 15;
    
    y1 >>= 15;
    y2 >>= 15;
    /* Never iterated if delta3f.y < 1 */
    for(; y1<=y2; ++y1){
      drawScanLine(buffer, width, height, y1, x1, x2, z1, z2, w1, w2, s1, s2, t1, t2);
      z1 += slope3Z;
      z2 += slope2Z;
      w1 += slope3W;
      w2 += slope2W;
      s1 += slope3S;
      s2 += slope2S;
      t1 += slope3T;
      t2 += slope2T;
      x1 += slope3X; /* bottom - middle */
      x2 += slope2X; /* bottom - top */
    }
  }
}
