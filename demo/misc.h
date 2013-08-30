#ifndef MISC_H_GUARD
#define MISC_H_GUARD
#include <string>

void SR_Init(int width, int height);
void SR_SetCaption(const std::string& title);
void SR_MainLoop(void (*cb_loop)(void*),
				 void (*cb_clean)(void*),
				 void* data);
#endif
