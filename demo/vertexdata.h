#ifndef VERTEXDATA_H_GUARD
#define VERTEXDATA_H_GUARD

#include <vector>
#include <linealg.h>

/* Render flags for SR_RENDER */
const unsigned int SR_TEXCOORD0 = 1;
const unsigned int SR_TEXCOORD1 = 2;
const unsigned int SR_LIGHTING = 4;
const unsigned int SR_COLOR = 8;

extern std::vector<Vector4f>* wc_vertices;
extern std::vector<Vector4f>* wc_tcoords0;
extern std::vector<Vector4f>* wc_tcoords1;
extern std::vector<Vector4f>* wc_normals;
extern std::vector<Vector4f>* wc_colors;

extern Matrix4f wc_modelview;
extern Matrix4f wc_projection;

void SR_SetVertices(std::vector<Vector4f>* vertices);
void SR_SetTexCoords0(std::vector<Vector4f>* tcoords0);
void SR_SetTexCoords1(std::vector<Vector4f>* tcoords1);
void SR_SetNormals(std::vector<Vector4f>* normals);
void SR_SetColors(std::vector<Vector4f>* normals);

void SR_SetModelViewMatrix(const Matrix4f& matrix);
void SR_SetProjectionMatrix(const Matrix4f& matrix);

#endif
