#include "vertexdata.h"

std::vector<Vector4f>* wc_vertices;
std::vector<Vector4f>* wc_tcoords0;
std::vector<Vector4f>* wc_tcoords1;
std::vector<Vector4f>* wc_normals;
std::vector<Vector4f>* wc_colors;

Matrix4f wc_modelview;
Matrix4f wc_projection;

void SR_SetVertices(std::vector<Vector4f>* vertices)
{
  wc_vertices = vertices;
}
void SR_SetTexCoords0(std::vector<Vector4f>* tcoords0)
{
  wc_tcoords0 = tcoords0;
}
void SR_SetTexCoords1(std::vector<Vector4f>& tcoords1)
{

}
void SR_SetNormals(std::vector<Vector4f>& normals)
{

}
void SR_SetColors(std::vector<Vector4f>& colors)
{

}

void SR_SetModelViewMatrix(const Matrix4f& matrix)
{
  wc_modelview = matrix;
}

void SR_SetProjectionMatrix(const Matrix4f& matrix)
{
  wc_projection = matrix;
}
