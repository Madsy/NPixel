#include <vector>
#include <linealg.h>
#include "meshgen.h"

void makeMeshSphere(std::vector<VectorPOD4f>& vertexData,
                    std::vector<VectorPOD4f>& tcoordData,
                    float radius)
{

    const int resolution = 100;
    const float halfPI = PI * 0.5f;
    float interp = 1.0f / (float)resolution;
    VectorPOD4f v0,v1,v2,v3;
    VectorPOD4f tc0,tc1,tc2,tc3;
    radius *= 0.5f;

    tc0.z = tc1.z = tc2.z = tc3.z = 0.0f;
    tc0.w = tc1.w = tc2.w = tc3.w = 0.0f;
    vertexData.reserve(resolution*resolution*6);
    tcoordData.reserve(resolution*resolution*6);

    for(int i=0; i<=resolution; ++i) {
        float theta0 = interp*(float)(i+0)*PI - halfPI;
        float theta1 = interp*(float)(i+1)*PI - halfPI;
        float z1 = std::sin(theta0);
        float z2 = std::sin(theta1);

        v0.z = z1 * radius;
        v1.z = z1 * radius;
        v2.z = z2 * radius;
        v3.z = z2 * radius;

        for(int j=0; j<=resolution; ++j) {
            float phi0 = interp*(float)(j+0)*2.0f*PI;
            float phi1 = interp*(float)(j+1)*2.0f*PI;
            float x1 = std::cos(theta0)*std::cos(phi0);
            float x2 = std::cos(theta0)*std::cos(phi1);
            float y1 = std::cos(theta0)*std::sin(phi0);
            float y2 = std::cos(theta0)*std::sin(phi1);
            float x3 = std::cos(theta1)*std::cos(phi0);
            float x4 = std::cos(theta1)*std::cos(phi1);
            float y3 = std::cos(theta1)*std::sin(phi0);
            float y4 = std::cos(theta1)*std::sin(phi1);

            float s1 = x1 * 0.5f + 0.5f;
            float s2 = x2 * 0.5f + 0.5f;
            float s3 = x3 * 0.5f + 0.5f;
            float s4 = x4 * 0.5f + 0.5f;
            float t1 = y1 * 0.5f + 0.5f;
            float t2 = y2 * 0.5f + 0.5f;
            float t3 = y3 * 0.5f + 0.5f;
            float t4 = y4 * 0.5f + 0.5f;

            v0.x = x1 * radius;
            v0.y = y1 * radius;
            v1.x = x2 * radius;
            v1.y = y2 * radius;
            v2.x = x3 * radius;
            v2.y = y3 * radius;
            v3.x = x4 * radius;
            v3.y = y4 * radius;

            tc0.x = s1;
            tc0.y = t1;
            tc1.x = s2;
            tc1.y = t2;
            tc2.x = s3;
            tc2.y = t3;
            tc3.x = s4;
            tc3.y = t4;

            vertexData.push_back(v0);
            vertexData.push_back(v1);
            vertexData.push_back(v2);
            vertexData.push_back(v2);
            vertexData.push_back(v1);
            vertexData.push_back(v3);
            tcoordData.push_back(tc0);
            tcoordData.push_back(tc1);
            tcoordData.push_back(tc2);
            tcoordData.push_back(tc2);
            tcoordData.push_back(tc1);
            tcoordData.push_back(tc3);
        }
    }
}

void makeMeshCircle(std::vector<VectorPOD4f>& dst, float radius)
{
    const int resolution = 64;
    float interp = 1.0f / (float)resolution;
    VectorPOD4f middle;
    VectorPOD4f v0, v1;

    middle.x = 0;
    middle.y = 0;
    middle.z = 0;
    middle.w = 1;

    dst.reserve(resolution*3);

    for(int i = 0; i<=resolution; ++i) {
        float alpha0 = (float)(i+0) * interp * PI * 2.0f;
        float alpha1 = (float)(i+1) * interp * PI * 2.0f;

        v0.x = sin(alpha0);
        v0.y = cos(alpha0);
        v1.x = sin(alpha1);
        v1.y = cos(alpha1);

        dst.push_back(v0);
        dst.push_back(middle);
        dst.push_back(v1);
    }
}

void makeMeshPlane(std::vector<VectorPOD4f>& vertexData,
                   std::vector<VectorPOD4f>& tcoordData,
                   float size)
{
    VectorPOD4f v0,v1,v2,v3;
    VectorPOD4f t0,t1,t2,t3;

    v0.x =  1.0f * size;
    v0.y =  1.0f * size;
    v0.z = 0.0f;
    v0.w = 1.0f;
    v1.x = -1.0f * size;
    v1.y =  1.0f * size;
    v1.z = 0.0f;
    v1.w = 1.0f;
    v2.x =  1.0f * size;
    v2.y = -1.0f * size;
    v2.z = 0.0f;
    v2.w = 1.0f;
    v3.x = -1.0f * size;
    v3.y = -1.0f * size;
    v3.z = 0.0f;
    v3.w = 1.0f;

    t0.x = 1.0f;
    t0.y = 1.0f;
    t0.z = 0.0f;
    t0.w = 1.0f;
    t1.x = 0.0f;
    t1.y = 1.0f;
    t1.z = 0.0f;
    t1.w = 1.0f;
    t2.x = 1.0f;
    t2.y = 0.0f;
    t2.z = 0.0f;
    t2.w = 1.0f;
    t3.x = 0.0f;
    t3.y = 0.0f;
    t3.z = 0.0f;
    t3.w = 1.0f;

    vertexData.push_back(v0);
    vertexData.push_back(v1);
    vertexData.push_back(v2);
    vertexData.push_back(v2);
    vertexData.push_back(v1);
    vertexData.push_back(v3);

    tcoordData.push_back(t0);
    tcoordData.push_back(t1);
    tcoordData.push_back(t2);
    tcoordData.push_back(t2);
    tcoordData.push_back(t1);
    tcoordData.push_back(t3);
}

void makeMeshCube(std::vector<VectorPOD4f>& vertexData,
                  std::vector<VectorPOD4f>& tcoordData,
                  float size)
{
    size *= 0.5f;

    VectorPOD4f v0,v1,v2,v3,v4,v5,v6,v7;
    VectorPOD4f t0,t1,t2,t3;

    v0.x =  1.0f * size;
    v0.y =  1.0f * size;
    v0.z = -1.0f * size;
    v0.w = 1.0f;
    v1.x = -1.0f * size;
    v1.y =  1.0f * size;
    v1.z = -1.0f * size;
    v1.w = 1.0f;
    v2.x =  1.0f * size;
    v2.y =  1.0f * size;
    v2.z =  1.0f * size;
    v2.w = 1.0f;
    v3.x = -1.0f * size;
    v3.y =  1.0f * size;
    v3.z =  1.0f * size;
    v3.w = 1.0f;
    v4.x =  1.0f * size;
    v4.y = -1.0f * size;
    v4.z = -1.0f * size;
    v4.w = 1.0f;
    v5.x = -1.0f * size;
    v5.y = -1.0f * size;
    v5.z = -1.0f * size;
    v5.w = 1.0f;
    v6.x =  1.0f * size;
    v6.y = -1.0f * size;
    v6.z =  1.0f * size;
    v6.w = 1.0f;
    v7.x = -1.0f * size;
    v7.y = -1.0f * size;
    v7.z =  1.0f * size;
    v7.w = 1.0f;

    t0.x = 1.0f;
    t0.y = 1.0f;
    t0.z = 0.0f;
    t0.w = 1.0f;
    t1.x = 0.0f;
    t1.y = 1.0f;
    t1.z = 0.0f;
    t1.w = 1.0f;
    t2.x = 1.0f;
    t2.y = 0.0f;
    t2.z = 0.0f;
    t2.w = 1.0f;
    t3.x = 0.0f;
    t3.y = 0.0f;
    t3.z = 0.0f;
    t3.w = 1.0f;

    //vertexData.reserve(6*6);
    //tcoordData.reserve(6*6);

    /* Top */
    vertexData.push_back(v0);
    vertexData.push_back(v1);
    vertexData.push_back(v2);
    vertexData.push_back(v2);
    vertexData.push_back(v1);
    vertexData.push_back(v3);

    /* Bottom */
    vertexData.push_back(v5);
    vertexData.push_back(v4);
    vertexData.push_back(v7);
    vertexData.push_back(v7);
    vertexData.push_back(v4);
    vertexData.push_back(v6);

    /* Right */
    vertexData.push_back(v0);
    vertexData.push_back(v2);
    vertexData.push_back(v4);
    vertexData.push_back(v4);
    vertexData.push_back(v2);
    vertexData.push_back(v6);

    /* Left */
    vertexData.push_back(v3);
    vertexData.push_back(v1);
    vertexData.push_back(v7);
    vertexData.push_back(v7);
    vertexData.push_back(v1);
    vertexData.push_back(v5);

    /* Back */
    vertexData.push_back(v1); //0
    vertexData.push_back(v0); //1
    vertexData.push_back(v5); //4
    vertexData.push_back(v5); //4
    vertexData.push_back(v0); //1
    vertexData.push_back(v4); //5

    /* Front */
    vertexData.push_back(v2);
    vertexData.push_back(v3);
    vertexData.push_back(v6);
    vertexData.push_back(v6);
    vertexData.push_back(v3);
    vertexData.push_back(v7);

    for(int i=0; i<6; ++i) {
        tcoordData.push_back(t0);
        tcoordData.push_back(t1);
        tcoordData.push_back(t2);
        tcoordData.push_back(t2);
        tcoordData.push_back(t1);
        tcoordData.push_back(t3);
    }
}
