#ifndef TEXTURE_H_GUARD
#define TEXTURE_H_GUARD
#include <string>
#include <vector>

struct Texture
{
    std::vector<unsigned int> texels;
    unsigned int width;
    unsigned int height;
};


const Texture* ReadPNG(const std::string& name);
void SR_BindTexture0(const Texture* texture);
void SR_BindTexture1(const Texture* texture);

extern const Texture* wc_texture0;
extern const Texture* wc_texture1;

#endif
