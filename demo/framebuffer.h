#ifndef FRAMEBUFFER_H_GUARD
#define FRAMEBUFFER_H_GUARD
#include <SDL/SDL.h>
#include <vector>
#include "buffer.h"

/* Width and height of framebuffer */
extern unsigned int wc_width;
extern unsigned int wc_height;
/* bits per pixel */
extern unsigned int wc_bpp;
/* alpha not supported if alpha shift is 0 */
extern bool wc_alphaSupported;
/* Masks for R, G, B, A. Support 15, 16 and 24 bit modes */
extern unsigned int wc_aMask;
extern unsigned int wc_rMask;
extern unsigned int wc_gMask;
extern unsigned int wc_bMask;
/* Shifts for R, G, B, A. Support 15, 16 and 24 bit modes */
extern unsigned int wc_aShift;
extern unsigned int wc_rShift;
extern unsigned int wc_gShift;
extern unsigned int wc_bShift;
/* Shifts for when converting from 24bpp or 32bpp to the internal mode
 The "loss" is the number of bits lost per channel*/
extern unsigned int wc_aLoss;
extern unsigned int wc_rLoss;
extern unsigned int wc_gLoss;
extern unsigned int wc_bLoss;


extern Buffer2D<unsigned int> wc_screenbuffer; //actual pointer to HW framebuffer (default)
extern Buffer2D<unsigned short> wc_screendepthbuffer; //default depth buffer
extern Buffer2D<unsigned int>* wc_colorbuffer;
extern Buffer2D<unsigned short>* wc_depthbuffer;

const unsigned int SR_COLOR_BUFFER=1;
const unsigned int SR_DEPTH_BUFFER=2;

void SR_InitBuffers(unsigned int width, unsigned int height);
void SR_ClearBuffer(unsigned int type);
void SR_BindDefaultBuffers();
void SR_BindBuffers(Buffer2D<unsigned int>* colorbuffer, Buffer2D<unsigned short>* depthbuffer);
void SR_Flip();
#endif
