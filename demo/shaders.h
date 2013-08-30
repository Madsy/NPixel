#ifndef SHADERS_H_GUARD
#define SHADERS_H_GUARD

#include "texture.h"

struct Coeff
{
  float A, B, C;
  Coeff(float a, float b, float c) : A(a), B(b), C(c){}
};


inline void fs_basic(unsigned int x, unsigned int y, unsigned int* pixel,
					 const Coeff& wCoeff, const Coeff& texCoordCoeff0U, const Coeff& texCoordCoeff0V)
{
  float wi, ui, vi, u,v,w;
  wi = wCoeff.A*x + wCoeff.B*y + wCoeff.C;
  ui = texCoordCoeff0U.A*x + texCoordCoeff0U.B*y + texCoordCoeff0U.C;
  vi = texCoordCoeff0V.A*x + texCoordCoeff0V.B*y + texCoordCoeff0V.C;
  if(wi != 0.0f) w = 1.0f / wi;
  else w = 0.0f;

  u = ui*w;
  v = vi*w;

  unsigned char a = 255;
  unsigned char r = u*255.0f;
  unsigned char g = v*255.0f;
  unsigned char b = 0;

  int ty = v*(float)wc_texture0->height;
  int tx = u*(float)wc_texture0->width;

  unsigned int texel = wc_texture0->texels[tx + ty * wc_texture0->width];

  //  *pixel = (unsigned int)(a << 24) | (r << 16) | (g << 8) | b;
  *pixel = texel;
}

#endif
