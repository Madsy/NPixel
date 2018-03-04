#include "texture.h"
#include <vector>
#include <IL/il.h>
#include <IL/ilu.h>

const Texture* wc_texture0;
const Texture* wc_texture1;

void SR_BindTexture0(const Texture* texture)
{
	wc_texture0 = texture;
}
void SR_BindTexture1(const Texture* texture)
{
	wc_texture1 = texture;
}


const Texture* ReadPNG(const std::string& name)
{
    Texture *texture;
    ILuint img;
    ilGenImages(1, &img);
    ilBindImage(img);
    if(!ilLoadImage(name.c_str())){
	    ilDeleteImages(1, &img);
        return nullptr;
    }
    iluFlipImage();

    texture = new Texture;
    texture->width = ilGetInteger(IL_IMAGE_WIDTH);
    texture->height = ilGetInteger(IL_IMAGE_HEIGHT);
    texture->texels.resize(texture->width * texture->height);
    ilCopyPixels(0, 0, 0, texture->width, texture->height, 1, IL_RGBA, IL_UNSIGNED_BYTE, &texture->texels[0]);
    ilDeleteImages(1, &img);

    return texture;
}

