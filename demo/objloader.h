#ifndef OBJLOADER_H_GUARD
#define OBJLOADER_H_GUARD
#include <string>
#include <vector>
#include <linealg.h>

bool importWaveFrontObjModel(const std::string& path,
							 std::vector<VectorPOD4f>& vertexData,
							 std::vector<VectorPOD4f>& tcoordData);

#endif
