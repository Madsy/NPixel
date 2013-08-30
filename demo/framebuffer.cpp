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
			  wc_colorbuffer.data.end(), 0);
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
  std::copy(wc_colorbuffer.data.begin(), wc_colorbuffer.data.end(), p);
  SDL_UnlockSurface(s);
  SDL_Flip(s);
}
