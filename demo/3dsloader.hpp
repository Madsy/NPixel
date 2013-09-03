#include <linealg.h>
#include <vector>
#include <string>

struct Mesh {
	std::vector<Vector4f> vertices;
	std::vector<Vector4f> tcoords;
	std::vector<Vector4f> normals;

	Vector4f ambient;
	Vector4f diffuse;
	Vector4f specular;
};

Mesh load3dsMesh(const std::string& filename);
