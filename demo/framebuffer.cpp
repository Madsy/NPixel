#include <algorithm>
#include "framebuffer.h"

/* Width and height of framebuffer */
unsigned int wc_width = 640;
unsigned int wc_height = 480;
/* bits per pixel */
unsigned int wc_bpp = 16;
/* alpha not supported if alpha shift is 0 */
bool wc_alphaSupported = false;
/* Masks for R, G, B, A. Support 15, 16 and 24 bit modes */
unsigned int wc_aMask = 0xFF000000;
unsigned int wc_rMask = 0x00FF0000;
unsigned int wc_gMask = 0x0000FF00;
unsigned int wc_bMask = 0x000000FF;
/* Shifts for R, G, B, A. Support 15, 16 and 24 bit modes */
unsigned int wc_aShift = 24;
unsigned int wc_rShift = 16;
unsigned int wc_gShift = 8;
unsigned int wc_bShift = 0;
/* Shifts for when converting from 24bpp or 32bpp to the internal mode
 The "loss" is the number of bits lost per channel*/
unsigned int wc_aLoss = 0;
unsigned int wc_rLoss = 0;
unsigned int wc_gLoss = 0;
unsigned int wc_bLoss = 0;

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
        void* p = wc_colorbuffer->Lock();
        //std::fill(wc_colorbuffer.data.begin(),
        //		  wc_colorbuffer.data.end(), 0xFFFFFFFF);
        int bytesperpixel = (wc_bpp + 7) / 8;
        unsigned int len = wc_colorbuffer->w * wc_colorbuffer->h * bytesperpixel;
        memset(p, 0, len);
        wc_colorbuffer->Unlock();
    }
    if((type & SR_DEPTH_BUFFER) && wc_depthbuffer != 0) {
        //std::fill(wc_depthbuffer.data.begin(),
        //		  wc_depthbuffer.data.end(), 65535);
        void* p = wc_depthbuffer->Ptr();
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
    //SDL_Flip(s);
    SDL_UpdateRect(s, 0, 0, 0, 0);
}
