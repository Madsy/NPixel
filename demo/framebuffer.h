#ifndef FRAMEBUFFER_H_GUARD
#define FRAMEBUFFER_H_GUARD
#include <SDL/SDL.h>
#include <vector>
#include "buffer.h"

extern Buffer2D<unsigned int> wc_colorbuffer;
extern Buffer2D<unsigned short> wc_depthbuffer;

const unsigned int SR_COLOR_BUFFER=1;
const unsigned int SR_DEPTH_BUFFER=2;

void SR_InitBuffers(unsigned int width, unsigned int height);
void SR_ClearBuffer(unsigned int type);
void SR_Flip();
#endif
