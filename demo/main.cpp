#include <SDL/SDL.h>
#include <vector>
#include <iterator>
#include <algorithm>
#include <linealg.h>
#include "framebuffer.h"
#include "rasterizer.h"
#include "vertexdata.h"
#include "meshgen.h"
#include "texture.h"
#include "myassert.h"
#include "misc.h"
#include "objloader.h"

struct Mesh {
    std::vector<VectorPOD4f> vertexData; // Our original mesh
    std::vector<VectorPOD4f> tcoordData; // Our original mesh
    float rotationSpeed;
    VectorPOD4f position;
};

unsigned int printAccum = 0;
int frameCount = 0;
MatrixPOD4f clipMatrix;

std::vector<VectorPOD4f> projVerts;
std::vector<VectorPOD4f> projTex;

float camRotX = 0.0f;
float camRotY = 0.0f;
VectorPOD4f camPosition = {0.0f, 0.0f, -20.0f, 1.0f};
const float camTransVelocity = 2.0f;
const float camRotVelocity = 12.0f;

static void updateCamera(MatrixPOD4f& camera, float t)
{
    MatrixPOD4f rotX, rotY, trans;
    VectorPOD4f forward, right;
    float deltaX, deltaY;
    Uint8 *keystate;

    translate(trans, camPosition);
    rotateX(rotX, camRotX);
    rotateY(rotY, camRotY);
    Mat4Mat4Mul(camera, rotX, rotY);
    Mat4Mat4Mul(camera, trans, camera);

#if 0
    forward.x = camera[2]; //2
    forward.y = camera[6]; //6
    forward.z = camera[10]; //10
    forward.w = 1.0f;
    right.x = camera[0]; //0
    right.y = camera[4]; //4
    right.z = camera[8]; //8
    right.w = 1.0f;
#else
    forward.x = 0.0f;
    forward.y = 0.0f;
    forward.z = -1.0f;
    forward.w = 1.0f;
    right.x = 1.0f;
    right.y = 0.0f;
    right.z = 0.0f;
    right.w = 1.0f;
#endif

    keystate = SDL_GetKeyState(NULL);
    t = 0.125f;
    if(keystate['w']) {
        camPosition.x -= forward.x * t * camTransVelocity;
        camPosition.y -= forward.y * t * camTransVelocity;
        camPosition.z -= forward.z * t * camTransVelocity;
        //printf("pos %3f %3f %3f\n", camPosition.x, camPosition.y, camPosition.z);
    } else if(keystate['s']) {
        camPosition.x += forward.x * t * camTransVelocity;
        camPosition.y += forward.y * t * camTransVelocity;
        camPosition.z += forward.z * t * camTransVelocity;
    }
    if(keystate['a']) {
        camPosition.x += right.x * t * camTransVelocity;
        camPosition.y += right.y * t * camTransVelocity;
        camPosition.z += right.z * t * camTransVelocity;
    } else if(keystate['d']) {
        camPosition.x -= right.x * t * camTransVelocity;
        camPosition.y -= right.y * t * camTransVelocity;
        camPosition.z -= right.z * t * camTransVelocity;
    }
    /*
    int x, y;
    (void)SDL_GetMouseState(&x, &y);
    int width = wc_colorbuffer->w;
    int height = wc_colorbuffer->h;
    camRotX = (float)(y - 240) / 240.0f * 90.0f;
    camRotY = (float)(x - 320) / 320.0f * 180.0f;
    camRotX = camRotX;
    */
    if(keystate[SDLK_UP]) {
        camRotX += (t * camRotVelocity);
    } else if(keystate[SDLK_DOWN]) {
        camRotX -= (t * camRotVelocity);
    }
    if(keystate[SDLK_LEFT]) {
        camRotY -= (t * camRotVelocity);
    } else if(keystate[SDLK_RIGHT]) {
        camRotY += (t * camRotVelocity);
    }
}

static void loop(void* data)
{
    int t = SDL_GetTicks();
    static int lt = SDL_GetTicks();
    float time_elapsed = static_cast<float>(t) * 0.001f;
    float tDelta = (float)(lt-t) * 0.001f;
    lt = t;
    if(tDelta > 1.0) tDelta = 1.0f;
    Mesh* mesh = static_cast<Mesh*>(data);
    const float fStep = 1.0f / (float)(1<<8);

    SR_ClearBuffer(SR_COLOR_BUFFER | SR_DEPTH_BUFFER );

    float rt = time_elapsed;

    projVerts.clear();
    projTex.clear();

    MatrixPOD4f worldMatrix, modelviewProjection;

    updateCamera(worldMatrix, tDelta);
    Mat4Mat4Mul(modelviewProjection, clipMatrix, worldMatrix);

    for(int j = 0; j < mesh->vertexData.size(); j+=3) {
        projVerts.push_back(Mat4Vec4Mul(modelviewProjection, mesh->vertexData[j + 0]));
        projVerts.push_back(Mat4Vec4Mul(modelviewProjection, mesh->vertexData[j + 1]));
        projVerts.push_back(Mat4Vec4Mul(modelviewProjection, mesh->vertexData[j + 2]));
        projTex.push_back(mesh->tcoordData[j + 0]);
        projTex.push_back(mesh->tcoordData[j + 1]);
        projTex.push_back(mesh->tcoordData[j + 2]);
    }

    SR_SetVertices(&projVerts);
    SR_SetTexCoords0(&projTex);
    SR_Render(SR_TEXCOORD0);
    SR_Flip();

#ifdef DEBUG
    float t2 = (float)SDL_GetTicks() * 0.001f;
    float tdelta = t2 - time_elapsed;
    float fps = 0.0f;
    if(tdelta != 0.0f) fps = 1.0f / tdelta;
    if(printAccum == 63) {
        printf("time: %f, fps: %f\n", tdelta, fps);
        printAccum = 0;
    }
    printAccum++;
#endif
    ++frameCount;
}

static void quit(void* data)
{

}

static float rnd_min_max(float mn, float mx)
{
    return mn + ((float)rand() / (float)RAND_MAX) * (mx - mn);
}

int main(int argc, char* argv[])
{
    (void)argc;
    (void)argv;
    const int width = 640;
    const int height = 480;
    const int depth = 32;
    perspective(clipMatrix, 90.0f, (float)width/(float)height, 0.1f, 40.0f);
    Mesh subway;
#if 1
    if(!importWaveFrontObjModel("subway2_3_HQ.obj", subway.vertexData, subway.tcoordData)) {
        printf("Couldn't load subway model\n");
        return 1;
    }
#else
    if(!importWaveFrontObjModel("testcube.obj", subway.vertexData, subway.tcoordData)) {
        printf("Couldn't load subway model\n");
        return 1;
    }
#endif

    ASSERT(subway.vertexData.size() == subway.tcoordData.size());
    printf("Rendering %d triangles.\n", subway.vertexData.size() / 3);

    SR_Init();
    SR_SetCaption("Tile-Rasterizer Test");

#if 1
    const Texture* tex = ReadPNG("subway2_2_baked.png");
#else
    const Texture* tex = ReadPNG("testcube.png");
#endif
    if(!tex) {
        printf("Couldn't load subway texture\n");
        return 1;
    }
    unsigned int timeStart = SDL_GetTicks();
    SR_BindTexture0(tex);
    SR_MainLoop(loop, quit, (void*)&subway);
    unsigned int timeEnd = SDL_GetTicks();
    printf("Rendered %d frames in %d seconds\n", frameCount, (timeEnd - timeStart) / 1000);
    return 0;
}


