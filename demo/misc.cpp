#include <SDL/SDL.h>
#include "misc.h"
#include "texture.h"
#include "framebuffer.h"

void SR_Init(int width, int height)
{
	SDL_Init(SDL_INIT_VIDEO);
	SDL_SetVideoMode(width, height, 32,
	                 SDL_DOUBLEBUF | SDL_HWSURFACE);
	SR_BindTexture0(NULL);
	SR_BindTexture1(NULL);
	SR_InitBuffers(width, height);
}

void SR_SetCaption(const std::string& title)
{
	SDL_WM_SetCaption(title.c_str(), NULL);
}
