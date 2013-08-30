#include "vertexdata.h"

std::vector<Vector4f> wc_vertices;
std::vector<Vector4f> wc_tcoords0;
std::vector<Vector4f> wc_tcoords1;
std::vector<Vector4f> wc_normals;
std::vector<Vector4f> wc_colors;

Matrix4f wc_modelview;
Matrix4f wc_projection;

void SR_SetVertices(const std::vector<Vector4f>& vertices)
{
  if(wc_vertices.size() != vertices.size())
	wc_vertices.resize(vertices.size());
  if(!vertices.empty())
	std::copy(vertices.begin(), vertices.end(), wc_vertices.begin());
}
void SR_SetTexCoords0(const std::vector<Vector4f>& tcoords0)
{
  if(wc_tcoords0.size() != tcoords0.size())
	wc_tcoords0.resize(tcoords0.size());
  if(!tcoords0.empty())
	std::copy(tcoords0.begin(), tcoords0.end(), wc_tcoords0.begin());
}
void SR_SetTexCoords1(const std::vector<Vector4f>& tcoords1)
{
  if(wc_tcoords1.size() != tcoords1.size())
	wc_tcoords1.resize(tcoords1.size());
  if(!tcoords1.empty())
	std::copy(tcoords1.begin(), tcoords1.end(), wc_tcoords1.begin());
}
void SR_SetNormals(const std::vector<Vector4f>& normals)
{
  if(wc_normals.size() != normals.size())
	wc_normals.resize(normals.size());
  if(!normals.empty())
	std::copy(normals.begin(), normals.end(), wc_normals.begin());
}
void SR_SetColors(const std::vector<Vector4f>& colors)
{
  if(wc_colors.size() != colors.size())
	wc_colors.resize(colors.size());
  if(!colors.empty())
	std::copy(colors.begin(), colors.end(), wc_colors.begin());
}

void SR_SetModelViewMatrix(const Matrix4f& matrix)
{
  wc_modelview = matrix;
}

void SR_SetProjectionMatrix(const Matrix4f& matrix)
{
  wc_projection = matrix;
}
