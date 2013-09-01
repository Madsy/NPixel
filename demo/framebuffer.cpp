#include <algorithm>
#include "framebuffer.h"

Buffer2D<unsigned int> wc_colorbuffer;
Buffer2D<unsigned short> wc_depthbuffer;


void SR_InitBuffers(unsigned int width, unsigned int height)
{
    wc_colorbuffer = Buffer2D<unsigned int>(width, height);
    wc_depthbuffer = Buffer2D<unsigned short>(width, height);
    return;
}

void SR_ClearBuffer(unsigned int type)
{
  if(type & SR_COLOR_BUFFER)
	std::fill(wc_colorbuffer.data.begin(),
			  wc_colorbuffer.data.end(), 0xFFFFFFFF);
  if(type & SR_DEPTH_BUFFER)
	std::fill(wc_depthbuffer.data.begin(),
			  wc_depthbuffer.data.end(), 65535);
  return;
}

void SR_Flip()
{
  SDL_Surface* s = SDL_GetVideoSurface();
  SDL_LockSurface(s);
  unsigned int* p = static_cast<unsigned int*>(s->pixels);
#if 1
  //std::copy(wc_colorbuffer.data.begin(), wc_colorbuffer.data.end(), p);
  memcpy(p, &wc_colorbuffer.data[0], sizeof(unsigned int) * wc_colorbuffer.data.size());
#else
  for(int y = 0; y < wc_depthbuffer.h; ++y){
	for(int x = 0; x < wc_depthbuffer.w; ++x){
	  unsigned short d = wc_depthbuffer[x + y*wc_depthbuffer.w];
	  unsigned char a = 255;
	  unsigned char b = 0;
	  unsigned char r = (d>>8);
	  unsigned char g = 0;//d&0xFF;
	  p[x + y*wc_colorbuffer.w] = (unsigned int)((a<<24) | (r<<16) | (r<<8) | r);
	}
  }
#endif
  SDL_UnlockSurface(s);
  SDL_Flip(s);
}
