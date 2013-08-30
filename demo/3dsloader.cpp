#include "3dsloader.hpp"
/*
struct Mesh
{
  std::vector<Vector4f> vertices;
  std::vector<Vector4f> tcoords;
  std::vector<Vector4f> normals;

  Vector4f ambient;
  Vector4f diffuse;
  Vector4f specular;
};
*/

Mesh load3dsMesh(const std::string& filename)
{
  Lib3dsFile* file = lib3ds_file_load(filename.c_str());
  if(!file) return;
  
  /* Really a linked list via Lib3dsMesh->next, but load only the first one
     for now.. */
  Lib3dsMesh* mesh = file->meshes;
}
