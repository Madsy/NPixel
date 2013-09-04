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
#define DEBUG
struct Mesh {
	std::vector<VectorPOD4f> vertexData; // Our original mesh
	std::vector<VectorPOD4f> tcoordData; // Our original mesh
	float rotationSpeed;
	VectorPOD4f position;
};

//float time_elapsed = 8.472656f;
unsigned int printAccum = 0;
const int NUM_MESHES = 100;

MatrixPOD4f clipMatrix;

std::vector<VectorPOD4f> projVerts;
std::vector<VectorPOD4f> projTex;

static void loop(void* data)
{
	static bool doOnce = true;
	unsigned int t = SDL_GetTicks();
	float time_elapsed = static_cast<float>(t) * 0.001f;
	Mesh* mesh = static_cast<Mesh*>(data);
	const float fStep = 1.0f / (float)(1<<8);
	if(SDL_GetKeyState(NULL)[SDLK_LEFT]) {
		time_elapsed -= fStep;
		printf("t %f\n", time_elapsed);
	}
	if(SDL_GetKeyState(NULL)[SDLK_RIGHT]) {
		time_elapsed += fStep;
		printf("t %f\n", time_elapsed);
	}

	SR_ClearBuffer(SR_COLOR_BUFFER | SR_DEPTH_BUFFER);

	float rt = time_elapsed;

	projVerts.clear();
	projTex.clear();

	for(int i = 0; i < NUM_MESHES; ++i) {
		float xOffset = 1.8f * std::sin(2.0f * M_PI * rt * mesh[i].rotationSpeed);
		VectorPOD4f offsetVec = {xOffset, 0.0f, 0.0f, 1.0f};

		MatrixPOD4f trans0, trans1, rotX, rotY, rotZ, worldMatrix, modelviewProjection;
		translate(trans0, offsetVec);
		translate(trans1, mesh[i].position);
		rotateX(rotX, 45.0f * rt * mesh[i].rotationSpeed);
		rotateY(rotY, 60.0f * rt * mesh[i].rotationSpeed);
		rotateZ(rotZ, 20.0f * rt * mesh[i].rotationSpeed);

		Mat4Mat4Mul(worldMatrix, rotY, rotZ);
		Mat4Mat4Mul(worldMatrix, rotX, worldMatrix);
		Mat4Mat4Mul(worldMatrix, trans1, worldMatrix);
		Mat4Mat4Mul(worldMatrix, trans0, worldMatrix);
		Mat4Mat4Mul(modelviewProjection, clipMatrix, worldMatrix);

		for(int j = 0; j < mesh[i].vertexData.size(); j+=3) {
			projVerts.push_back(Mat4Vec4Mul(modelviewProjection, mesh[i].vertexData[j + 0]));
			projVerts.push_back(Mat4Vec4Mul(modelviewProjection, mesh[i].vertexData[j + 1]));
			projVerts.push_back(Mat4Vec4Mul(modelviewProjection, mesh[i].vertexData[j + 2]));
			projTex.push_back(mesh[i].tcoordData[j + 0]);
			projTex.push_back(mesh[i].tcoordData[j + 1]);
			projTex.push_back(mesh[i].tcoordData[j + 2]);
		}
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

	//srand(time(NULL));
	perspective(clipMatrix, 60.0f, (float)width/(float)height, 1.0f, 40.0f);

	for(int i = 0; i < NUM_MESHES; ++i) {
		mesh[i].rotationSpeed = rnd_min_max(0.0f, 0.25f);
		mesh[i].position.x = rnd_min_max(-1.0f, 1.0f);
		mesh[i].position.y = rnd_min_max(-1.0f, 1.0f);
		mesh[i].position.z = rnd_min_max(-2.5f, -30.0f);
		mesh[i].position.w = 1.0f;
	}
	SR_Init(width, height);
	SR_SetCaption("Tile-Rasterizer Test");

	const Texture* tex = ReadPNG("texture0.png");
	SR_BindTexture0(tex);

	for(int i=0; i<NUM_MESHES; ++i)
		makeMeshCube(mesh[i].vertexData, mesh[i].tcoordData, 1.0f);

	SR_MainLoop(loop, quit, (void*)&mesh[0]);
}
