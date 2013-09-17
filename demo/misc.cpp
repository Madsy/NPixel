#include <SDL/SDL.h>
#include "misc.h"
#include "texture.h"
#include "framebuffer.h"
#include "myassert.h"

static void printVideoInfo()
{
    const SDL_VideoInfo* info = SDL_GetVideoInfo();
    printf("HW available? %d\n", info->hw_available);
    printf("WM available? %d\n", info->wm_available);
    printf("HW blitting supported? %d\n", info->blit_hw);
    printf("HW alpha blitting accelerated? %d\n", info->blit_hw_A);
    printf("SW -> HW blits accelerated? %d\n", info->blit_sw);
    printf("SW -> HW alpha blits accelerated? %d\n", info->blit_sw_A);
    if(info->hw_available)
        printf("VRAM available: %d KiB\n", info->video_mem);
    printf("Width x Height: %d X %d\n", info->current_w, info->current_h);
    if(info->vfmt) {
        const SDL_PixelFormat* fmt = info->vfmt;
        printf("-----Video Format:-----\n");
        printf("Bits per pixel: %d\n", fmt->BitsPerPixel);
        printf("Bytes per pixel: %d\n", fmt->BytesPerPixel);
        printf("RGBA loss: %d %d %d %d\n", fmt->Rloss, fmt->Gloss, fmt->Bloss, fmt->Aloss);
        printf("RGBA shift: %d %d %d %d\n", fmt->Rshift, fmt->Gshift, fmt->Bshift, fmt->Ashift);
        printf("RGBA mask: %X %X %X %X\n", fmt->Rmask, fmt->Gmask, fmt->Bmask, fmt->Amask);
        //printf("Colorkey: %X\n". fmt->colorkey);
    }
}

static int getSupportedResolutions(SDL_Rect*** rects)
{
    const SDL_VideoInfo* vinfo = SDL_GetVideoInfo();
    *rects = SDL_ListModes(vinfo->vfmt, SDL_SWSURFACE | SDL_FULLSCREEN);
    if((int)(*rects) == -1) {
        *rects = 0;
        return -1;
    }
    if((*rects) == 0) return 0;
    int len = 0;
    for(; (*rects)[len]; ++len);
    return len;
}

/* A small misnormer. It selects the resolution, but also fills in color masks and color shift values.
   See framebuffer.h for wc_ globals. */
static bool selectResolution()
{
    /*
    Raspberry Pi modes: (RGB565 16-bit, single-buffering)

    [0] 1600 X 1200 (broken)
    [1] 1280 X 1024 (broken)
    [2] 1024 X 1024 (broken)
    [3] 1280 X 960  (broken)
    [4] 1152 X 864 (broken)
    [5] 1024 X 768 (broken)
    [6] 1184 X 624 (works)
    [7] 800 X 600 (works)
    [8] 768 X 576 (works)
    [9] 640 X 480 (works)
     */
    SDL_Rect** rects;
    int len = getSupportedResolutions(&rects);
    if(len == -1) {
        printf("Using default HD resolution.\n");
        return true;
    }
    if(len == 0) {
        printf("No screen modes supported!\n");
        return false;
    }

    bool legalChoice = false;
    unsigned int choice = 0;
    while(!legalChoice) {
        printf("Pick a screen mode: \n");
        for(int i = 0; i < len; ++i) {
            printf("[%d]: %d X %d\n", i, rects[i]->w, rects[i]->h);
        }
        printf("Choice: ");
        int numItems = scanf("%u", &choice);
        if(numItems == 1 && choice < len) legalChoice = true;
        else printf("Illegal choice. Entries range from 0 to %d\n", len-1);
    }
    /* Set globals for width and height */
    wc_width = rects[choice]->w;
    wc_height = rects[choice]->h;

    const SDL_VideoInfo* vinfo = SDL_GetVideoInfo();
    ASSERT(vinfo != 0);
    const SDL_PixelFormat* fmt = vinfo->vfmt;
    ASSERT(fmt != 0);

    /* Bits per pixel. Use 32 (or 24?) for Desktop, 16 for GP2X and Raspberry Pi */
    wc_bpp = fmt->BitsPerPixel;
    /* alpha not supported if alpha shift is 0 */
    wc_alphaSupported =
        (fmt->Ashift == 0) &&
        ((fmt->Rshift == 0) || (fmt->Gshift == 0) || (fmt->Bshift == 0));

    wc_aLoss = fmt->Aloss;
    wc_rLoss = fmt->Rloss;
    wc_gLoss = fmt->Gloss;
    wc_bLoss = fmt->Bloss;

    /* Masks for R, G, B, A. Support 15, 16 and 24 bit modes */
    wc_aMask = fmt->Amask;
    wc_rMask = fmt->Rmask;
    wc_gMask = fmt->Gmask;
    wc_bMask = fmt->Bmask;
    /* Shifts for R, G, B, A. Support 15, 16 and 24 bit modes */
    wc_aShift = fmt->Ashift;
    wc_rShift = fmt->Rshift;
    wc_gShift = fmt->Gshift;
    wc_bShift = fmt->Bshift;

    return true;
}


void SR_Init()
{
    if(SDL_Init(SDL_INIT_VIDEO) < 0){
        printf("Failed to initialize SDL!\n");
        exit(1);
    }

    #ifdef DEBUG
        printVideoInfo();
    #endif
    /* TODO: Change selectResolution to use a GUI on platforms that don't start it
       from the terminal */
    if(!selectResolution()) {
        SDL_Quit();
        exit(1);
    }

    SDL_SetVideoMode(wc_width, wc_height, wc_bpp,
                     SDL_FULLSCREEN | SDL_SWSURFACE);
    SR_BindTexture0(NULL);
    SR_BindTexture1(NULL);
    SR_InitBuffers(wc_width, wc_height);
}

void SR_SetCaption(const std::string& title)
{
    SDL_WM_SetCaption(title.c_str(), NULL);
}
