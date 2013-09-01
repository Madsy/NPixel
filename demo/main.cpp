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

struct Mesh
{
  std::vector<Vector4f> vertexData; // Our original mesh
  std::vector<Vector4f> tcoordData; // Our original mesh
  float rotationSpeed;
  Vector4f position;
};

//float time_elapsed = 8.472656f;
unsigned int printAccum = 0;
const int NUM_MESHES = 100;

static void loop(void* data)
{
  unsigned int t = SDL_GetTicks();
  float time_elapsed = static_cast<float>(t) * 0.001f;
  Mesh* mesh = static_cast<Mesh*>(data);
  const float fStep = 1.0f / (float)(1<<8);
  if(SDL_GetKeyState(NULL)[SDLK_LEFT]){
	time_elapsed -= fStep;
	printf("t %f\n", time_elapsed);
  }
  if(SDL_GetKeyState(NULL)[SDLK_RIGHT]){
	time_elapsed += fStep;
	printf("t %f\n", time_elapsed);
  }

  SR_ClearBuffer(SR_COLOR_BUFFER | SR_DEPTH_BUFFER);

  for(int i = 0; i < NUM_MESHES; ++i){
	float xOffset = 1.8f * std::sin(2.0f * M_PI * time_elapsed * mesh[i].rotationSpeed);
	
	Matrix4f worldMatrix =
	  translate(Vector4f(xOffset, 0.0f, 0.0f, 1.0f)) * 
	  translate(mesh[i].position) * 
	  rotateX(45.0f * time_elapsed * mesh[i].rotationSpeed) *
	  rotateY(60.0f * time_elapsed * mesh[i].rotationSpeed) *
	  rotateZ(20.0f * time_elapsed * mesh[i].rotationSpeed);

	SR_SetModelViewMatrix(worldMatrix);
	SR_SetVertices(mesh[i].vertexData);
	SR_SetTexCoords0(mesh[i].tcoordData);
	SR_Render(SR_TEXCOORD0);
  }
  SR_Flip();
  float t2 = (float)SDL_GetTicks() * 0.001f;
  float tdiff = t2 - time_elapsed;
  if(tdiff != 0.0f) tdiff = 1.0f / tdiff;
  if(printAccum == 63){
	printf("fps: %f\n", tdiff);
	printAccum = 0;
  }
  printAccum++;
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
  Mesh mesh[NUM_MESHES];

  srand(time(NULL));

  for(int i = 0; i < NUM_MESHES; ++i){
	mesh[i].rotationSpeed = rnd_min_max(0.0f, 0.25f);
	mesh[i].position = Vector4f(rnd_min_max(-1.0f, 1.0f), rnd_min_max(-1.0f, 1.0f), rnd_min_max(-2.5f, -30.0), 1.0f);
  }
  SR_Init(width, height);
  SR_SetCaption("Tile-Rasterizer Test");

  const Texture* tex = ReadPNG("texture0.png");
  SR_BindTexture0(tex);

  Matrix4f clipMatrix = perspective(60.0f, (float)width/(float)height, 1.0f, 40.0f);
  SR_SetProjectionMatrix(clipMatrix);

  for(int i=0; i<NUM_MESHES; ++i)
    makeMeshCube(mesh[i].vertexData, mesh[i].tcoordData, 1.0f);

  SR_MainLoop(loop, quit, (void*)&mesh[0]);
}
