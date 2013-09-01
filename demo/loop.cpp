#include "misc.h"
#include "framebuffer.h"
#include <SDL/SDL.h>

static void SR_Quit()
{
  SDL_Quit();
  exit(0);
}

void SR_MainLoop(void (*cb_loop)(void*),
				 void (*cb_clean)(void*),
				 void* data)
{
  bool running = true;
  SDL_Event event;
  //SDL_ShowCursor(0);
  while(running){
    while(SDL_PollEvent(&event)){
      switch(event.type){
		case SDL_KEYDOWN:
		  if(event.key.keysym.sym == SDLK_ESCAPE){
			running = false;
		  }
		  break;
		case SDL_QUIT:
		  running = false;
		  break;
	  }
    }
	cb_loop(data);
  }
  cb_clean(data);
  //  SDL_ShowCursor(1);
  SR_Quit();
}

