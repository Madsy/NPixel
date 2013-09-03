#ifndef MESHGEN_GUARD_H
#define MESHGEN_GUARD_H
#include <linealg.h>

void makeMeshSphere(std::vector<VectorPOD4f>& vertexData,
                    std::vector<VectorPOD4f>& tcoordData,
                    float radius);
void makeMeshCircle(std::vector<VectorPOD4f>& dst, float radius);
void makeMeshPlane(std::vector<VectorPOD4f>& vertexData,
                   std::vector<VectorPOD4f>& tcoordData,
                   float size);
void makeMeshCube(std::vector<VectorPOD4f>& vertexData,
                  std::vector<VectorPOD4f>& tcoordData,
                  float size);
#endif
