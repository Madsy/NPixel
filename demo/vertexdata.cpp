#include "vertexdata.h"

std::vector<VectorPOD4f>* wc_vertices;
std::vector<VectorPOD4f>* wc_tcoords0;
std::vector<VectorPOD4f>* wc_tcoords1;
std::vector<VectorPOD4f>* wc_normals;
std::vector<VectorPOD4f>* wc_colors;

Matrix4f wc_modelview;
Matrix4f wc_projection;

void SR_SetVertices(std::vector<VectorPOD4f>* vertices)
{
    wc_vertices = vertices;
}
void SR_SetTexCoords0(std::vector<VectorPOD4f>* tcoords0)
{
    wc_tcoords0 = tcoords0;
}
void SR_SetTexCoords1(std::vector<VectorPOD4f>* tcoords1)
{
    wc_tcoords1 = tcoords1;
}
void SR_SetNormals(std::vector<VectorPOD4f>* normals)
{
    wc_normals = normals;
}
void SR_SetColors(std::vector<VectorPOD4f>* colors)
{
    wc_colors = colors;
}

void SR_SetModelViewMatrix(const Matrix4f& matrix)
{
    wc_modelview = matrix;
}

void SR_SetProjectionMatrix(const Matrix4f& matrix)
{
    wc_projection = matrix;
}
