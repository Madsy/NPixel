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

  for(int i = 0; i < 3; ++i){
	float xOffset = 1.8f * std::sin(2.0f * M_PI * time_elapsed * mesh[i].rotationSpeed);
	
	Matrix4f worldMatrix =
	  translate(Vector4f(xOffset, 0.0f, 0.0f, 1.0f)) * 
	  translate(mesh[i].position) * 
	  rotateX(0.0f *   time_elapsed) *
	  rotateY(60.0f *   time_elapsed) *
	  rotateZ(0.0f * time_elapsed);

	SR_SetModelViewMatrix(worldMatrix);
	SR_SetVertices(mesh[i].vertexData);
	SR_SetTexCoords0(mesh[i].tcoordData);
	SR_Render(SR_TEXCOORD0);
  }
  SR_Flip();
  float t2 = (float)SDL_GetTicks() * 0.001f;
  float tdiff = t2 - time_elapsed;
  if(tdiff != 0.0f) tdiff = 1.0f / tdiff;
  if(printAccum == 0xFF){
	printf("fps: %f\n", tdiff);
	printAccum = 0;
  }
  printAccum++;
}

static void quit(void* data)
{

}

int main(int argc, char* argv[])
{
  (void)argc;
  (void)argv;
  const int NUM_MESHES = 3;
  const int width = 640;
  const int height = 360;
  const int depth = 32;
  Mesh mesh[NUM_MESHES];

  mesh[0].rotationSpeed = 0.1f;
  mesh[1].rotationSpeed = 0.03f;
  mesh[2].rotationSpeed = 0.07f;
  mesh[0].position = Vector4f(0.0f, 0.0f, -2.5f, 1.0f);
  mesh[1].position = Vector4f(0.1f, 0.0f, -2.0f, 1.0f);
  mesh[2].position = Vector4f(0.25f, 0.0f,-2.0f, 1.0f);

  SR_Init(width, height);
  SR_SetCaption("Tile-Rasterizer Test");

  const Texture* tex = ReadPNG("texture0.png");
  SR_BindTexture0(tex);

  Matrix4f clipMatrix = perspective(45.0f, 16.0f/9.0f, 1.0f, 100.0f);
  SR_SetProjectionMatrix(clipMatrix);

  for(int i=0; i<NUM_MESHES; ++i)
    makeMeshCube(mesh[i].vertexData, mesh[i].tcoordData, 1.0f);

  SR_MainLoop(loop, quit, (void*)&mesh[0]);
}


/*
int main(int argc, char* argv[])
{
  (void)argc;
  (void)argv;
  const int NUM_MESHES = 3;
  const int width = 640;
  const int height = 360;
  const int depth = 32;
  bool running = true;
  SDL_Event event;
  Mesh mesh[NUM_MESHES];
  std::vector<unsigned int> texBuf;   // RGBA Texture image

  mesh[0].rotationSpeed = 0.1f;
  mesh[1].rotationSpeed = 0.03f;
  mesh[2].rotationSpeed = 0.07f;
  mesh[0].position = Vector4f(0.0f, 0.0f, -2.5f, 1.0f);
  mesh[1].position = Vector4f(0.1f, 0.0f, -2.0f, 1.0f);
  mesh[2].position = Vector4f(0.25f, 0.0f,-2.0f, 1.0f);

  
  SR_Init(width, height);
  SR_SetCaption("TRSI / Primitive Coop-Demo");
  for(int i=0; i<NUM_MESHES; ++i)  
    makeMeshCube(mesh[i].vertexData, mesh[i].tcoordData, 1.0f);
  const Texture* texture = ReadPNG("texture0.png");
  if(!texture){
    printf("Couldn't load one or more texture maps.\n");
    return -1;
  }
  SR_BindTexture0(texture);

  int tricount = 0;
  int tricountFrame = 0;
  int frames = 0;
  float taccum = 0;

  while(running){
    while(SDL_PollEvent(&event)){
      switch(event.type)
		{
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
    float time_elapsed = (float)SDL_GetTicks() * 0.001f;
    SDL_LockSurface(screen);
    unsigned int* pixels = static_cast<unsigned int*>(screen->pixels);
	
    // perspective function is in linealg.h under /include
    Matrix4f clipMatrix = perspective(45.0f, 16.0f/9.0f, 1.0f, 10.0f);
    // Clear our depth buffer
    //SR_ClearBuffer(SR_DEPTH_BUFFER | SR_COLOR_BUFFER);
    memset(&depthbuffer.data[0], 65535, sizeof(Uint16) * width * height);
    // clear the screen to black
    memset(pixels, 0, sizeof(Uint32) * width * height);

    tricountFrame = 0;
    for(int k=0; k<NUM_MESHES; ++k){
      // world matrix transform
      float xOffset = 1.8f * std::sin(2.0f * M_PI * time_elapsed * mesh[k].rotationSpeed);
      Matrix4f worldMatrix =
		translate(Vector4f(xOffset, 0.0f, 0.0f, 1.0f)) * 
		translate(mesh[k].position) * 
		rotateX(0.0f *   time_elapsed) *
		rotateY(60.0f *   time_elapsed) *
		rotateZ(0.0f * time_elapsed);
      Matrix4f worldClipMatrix = clipMatrix * worldMatrix;

      // We need a new working copy every frame
      mesh[k].workingCopyVertex = mesh[k].vertexData;
      mesh[k].workingCopyTCoord = mesh[k].tcoordData;

      // Transform our points
      for(unsigned int i=0; i<mesh[k].workingCopyVertex.size(); i+=3){
		Vector4f& p1 = mesh[k].workingCopyVertex[i];
		Vector4f& p2 = mesh[k].workingCopyVertex[i+1];
		Vector4f& p3 = mesh[k].workingCopyVertex[i+2];
		p1 = worldClipMatrix * p1;
		p2 = worldClipMatrix * p2;
		p3 = worldClipMatrix * p3;
      }

      // Clip against the six frustum planes
      clip_triangle(mesh[k].workingCopyVertex,
					mesh[k].workingCopyTCoord,
					Vector4f(-1.0f,  0.0f, 0.0f, 1.0f));
      clip_triangle(mesh[k].workingCopyVertex,
					mesh[k].workingCopyTCoord,
					Vector4f( 1.0f,  0.0f, 0.0f, 1.0f));
      clip_triangle(mesh[k].workingCopyVertex,
					mesh[k].workingCopyTCoord,
					Vector4f( 0.0f,  1.0f, 0.0f, 1.0f));
      clip_triangle(mesh[k].workingCopyVertex,
					mesh[k].workingCopyTCoord,
					Vector4f( 0.0f, -1.0f, 0.0f, 1.0f));
      clip_triangle(mesh[k].workingCopyVertex,
					mesh[k].workingCopyTCoord,
					Vector4f( 0.0f,  0.0f,-1.0f, 1.0f));
      clip_triangle(mesh[k].workingCopyVertex,
					mesh[k].workingCopyTCoord,
					Vector4f( 0.0f,  0.0f, 1.0f, 1.0f));

      for(unsigned int j=0; j<mesh[k].workingCopyVertex.size(); ++j){
		float wInv = 1.0f / mesh[k].workingCopyVertex[j].w;
		// perspective divide for tcoords
		mesh[k].workingCopyTCoord[j].x *= wInv;
		mesh[k].workingCopyTCoord[j].y *= wInv;
	
		// project function is in linealg.h under /include
		//x and y is in screenspace
		//z is normalized into [0,1> range
		//w' = 1.0 / w

		mesh[k].workingCopyVertex[j] =
		  project(mesh[k].workingCopyVertex[j],
				  (float)width, (float)height);
      }
      

      //*******************************************************************
      //************* TODO: Clean up this backface-culling code! **********
	  //*******************************************************************

      std::vector<Vector4f> tmpV, tmpTC;
      tmpV.reserve(mesh[k].workingCopyVertex.size());
      tmpTC.reserve(mesh[k].workingCopyVertex.size());
      for(unsigned int j=0; j<mesh[k].workingCopyVertex.size(); j+=3){
		Vector4f& v1 = mesh[k].workingCopyVertex[j+0];
		Vector4f& v2 = mesh[k].workingCopyVertex[j+1];
		Vector4f& v3 = mesh[k].workingCopyVertex[j+2];
		Vector4f& tc1 = mesh[k].workingCopyTCoord[j+0];
		Vector4f& tc2 = mesh[k].workingCopyTCoord[j+1];
		Vector4f& tc3 = mesh[k].workingCopyTCoord[j+2];
		Vector4f d1 = v2 - v1;
		Vector4f d2 = v3 - v1;	
		if((d1.x * d2.y) - (d2.x * d1.y) > 0){
		  tmpV.push_back(v1);
		  tmpV.push_back(v2); 
		  tmpV.push_back(v3);
		  tmpTC.push_back(tc1);
		  tmpTC.push_back(tc2);
		  tmpTC.push_back(tc3);
		}
      }      
      mesh[k].workingCopyVertex = tmpV;
      mesh[k].workingCopyTCoord = tmpTC;

      //*******************************************************************
      //************* TODO: Clean up this backface-culling code! **********
      //*******************************************************************

      // Draw the triangles
      DrawTriangle(mesh[k].workingCopyVertex,
				   mesh[k].workingCopyTCoord,
				   pixels, width, height);
      tricountFrame += (mesh[k].workingCopyVertex.size() / 3);
    }

    tricount += tricountFrame;
    SDL_UnlockSurface(screen);
    SDL_Flip(screen);
    ++frames;

    float time_elapsed2 = (float)SDL_GetTicks() * 0.001f;
    float dt = time_elapsed2 - time_elapsed;
    taccum += dt;
    if(taccum >= 1.0f){
      printf("Rendered %d triangles in 1 second.\n", tricount);
      printf("Rendered %d frames in %f seconds.\n", frames, taccum);
      taccum = 0;
      frames = 0;
      tricount = 0;
    }
  }    
  SDL_Quit();
  return 0;
}
*/
