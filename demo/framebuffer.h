#ifndef FRAMEBUFFER_H_GUARD
#define FRAMEBUFFER_H_GUARD
#include <SDL/SDL.h>
#include <vector>

template<typename T> struct Buffer2D
{
    Buffer2D(unsigned int width, unsigned int height) : w(width), h(height), data(width*height)
    {
	  data.resize(width*height);
	}
    Buffer2D() : w(0), h(0), data(){}
    inline T& operator[](size_t index){ return data[index]; }
    inline const T& operator[](size_t index) const { return data[index]; }

    unsigned int w;
    unsigned int h;
    std::vector<T> data;
};

extern Buffer2D<unsigned int> wc_colorbuffer;
extern Buffer2D<unsigned short> wc_depthbuffer;

const unsigned int SR_COLOR_BUFFER=1;
const unsigned int SR_DEPTH_BUFFER=2;

void SR_InitBuffers(unsigned int width, unsigned int height);
void SR_ClearBuffer(unsigned int type);
void SR_Flip();
#endif
