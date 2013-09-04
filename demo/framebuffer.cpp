#include <algorithm>
#include "framebuffer.h"

Buffer2D<unsigned int> wc_screenbuffer; //actual pointer to HW framebuffer (default)
Buffer2D<unsigned short> wc_screendepthbuffer; //default depth buffer
Buffer2D<unsigned int>* wc_colorbuffer; //points to current buffer
Buffer2D<unsigned short>* wc_depthbuffer; //points to current


void SR_InitBuffers(unsigned int width, unsigned int height)
{
	wc_screenbuffer = Buffer2D<unsigned int>(width, height, 0, true);
	wc_screendepthbuffer = Buffer2D<unsigned short>(width, height, 0, false);
	SR_BindDefaultBuffers();
	return;
}

void SR_ClearBuffer(unsigned int type)
{
	if((type & SR_COLOR_BUFFER) && wc_colorbuffer != 0) {
		unsigned int* p = wc_colorbuffer->Lock();
		//std::fill(wc_colorbuffer.data.begin(),
		//		  wc_colorbuffer.data.end(), 0xFFFFFFFF);
		unsigned int len = wc_colorbuffer->w * wc_colorbuffer->h * sizeof(unsigned int);
		memset(p, 0, len);
		wc_colorbuffer->Unlock();
	}
	if((type & SR_DEPTH_BUFFER) && wc_depthbuffer != 0) {
		//std::fill(wc_depthbuffer.data.begin(),
		//		  wc_depthbuffer.data.end(), 65535);
		unsigned short* p = wc_depthbuffer->Ptr();
		unsigned int len = wc_depthbuffer->w * wc_depthbuffer->h * sizeof(unsigned short);
		memset(p, 255, len);
	}
	return;
}

void SR_BindDefaultBuffers()
{
	wc_colorbuffer = &wc_screenbuffer;
	wc_depthbuffer = &wc_screendepthbuffer;
}

/* This allows you to render to textures */
void SR_BindBuffers(Buffer2D<unsigned int>* colorbuffer, Buffer2D<unsigned short>* depthbuffer)
{
	wc_colorbuffer = colorbuffer;
	wc_depthbuffer = depthbuffer;
}

void SR_Flip()
{
	/*
	SDL_LockSurface(s);
	unsigned int* p = static_cast<unsigned int*>(s->pixels);
	*/
	SDL_Surface* s = SDL_GetVideoSurface();
	//std::copy(wc_colorbuffer.data.begin(), wc_colorbuffer.data.end(), p);
	//wc_colorbuffer->Lock();
	//unsigned int* p = wc_colorbuffer->Ptr();
	//memcpy(p, &wc_colorbuffer.data[0], sizeof(unsigned int) * wc_colorbuffer.data.size());
	//wc_colorbuffer->Unlock();
	/*
	SDL_UnlockSurface(s);
	*/
	SDL_Flip(s);
}
