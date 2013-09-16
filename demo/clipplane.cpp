#include <vector>
#include <algorithm>
#include <cstdio>
#include <linealg.h>
#include "clipplane.h"
#include "vertexdata.h"

/*
	wc_vertices
	wc_tcoords0
	wc_tcoords1
	wc_normals
	wc_colors
*/
void clip_triangle(Vector4f plane, unsigned int flags)
{
#if 1
	unsigned int edge0, edge1;
	unsigned int tri;

	size_t vertexListSize = wc_vertices->size();

	for(tri=0; tri<vertexListSize; tri+=3) {
		std::vector<VectorPOD4f> polygon; /* outputted vertices */
		std::vector<VectorPOD4f> polygon_tcoords0; /* interpolants */
		std::vector<VectorPOD4f> polygon_tcoords1;
		std::vector<VectorPOD4f> polygon_normals;
		std::vector<VectorPOD4f> polygon_colors;

		for(edge0=tri+2, edge1=tri; edge1 < tri+3; edge0 = edge1++) {
			VectorPOD4f point_current, point_next, tcoord0_current,
			         tcoord0_next, tcoord1_current, tcoord1_next, normal_current,
			         normal_next, color_current, color_next;

			point_current.x = (*wc_vertices)[edge0].x;
			point_current.y = (*wc_vertices)[edge0].y;
			point_current.z = (*wc_vertices)[edge0].z;
			point_current.w = (*wc_vertices)[edge0].w;
			point_next.x     = (*wc_vertices)[edge1].x;
			point_next.y     = (*wc_vertices)[edge1].y;
			point_next.z     = (*wc_vertices)[edge1].z;
			point_next.w     = (*wc_vertices)[edge1].w;

			if(flags & SR_TEXCOORD0) {
				tcoord0_current.x = (*wc_tcoords0)[edge0].x;
				tcoord0_current.y = (*wc_tcoords0)[edge0].y;
				tcoord0_current.z = (*wc_tcoords0)[edge0].z;
				tcoord0_current.w = (*wc_tcoords0)[edge0].w;
				tcoord0_next.x    = (*wc_tcoords0)[edge1].x;
				tcoord0_next.y    = (*wc_tcoords0)[edge1].y;
				tcoord0_next.z    = (*wc_tcoords0)[edge1].z;
				tcoord0_next.w    = (*wc_tcoords0)[edge1].w;
			}
			if(flags & SR_TEXCOORD1) {
				tcoord1_current.x = (*wc_tcoords1)[edge0].x;
				tcoord1_current.y = (*wc_tcoords1)[edge0].y;
				tcoord1_current.z = (*wc_tcoords1)[edge0].z;
				tcoord1_current.w = (*wc_tcoords1)[edge0].w;
				tcoord1_next.x    = (*wc_tcoords1)[edge1].x;
				tcoord1_next.y    = (*wc_tcoords1)[edge1].y;
				tcoord1_next.z    = (*wc_tcoords1)[edge1].z;
				tcoord1_next.w    = (*wc_tcoords1)[edge1].w;
			}
			if(flags & SR_LIGHTING) {
				normal_current.x = (*wc_normals)[edge0].x;
				normal_current.y = (*wc_normals)[edge0].y;
				normal_current.z = (*wc_normals)[edge0].z;
				normal_current.w = (*wc_normals)[edge0].w;
				normal_next.x    = (*wc_normals)[edge1].x;
				normal_next.y    = (*wc_normals)[edge1].y;
				normal_next.z    = (*wc_normals)[edge1].z;
				normal_next.w    = (*wc_normals)[edge1].w;
			}
			if(flags & SR_COLOR) {
				color_current.x = (*wc_colors)[edge0].x;
				color_current.y = (*wc_colors)[edge0].y;
				color_current.z = (*wc_colors)[edge0].z;
				color_current.w = (*wc_colors)[edge0].w;
				color_next.x    = (*wc_colors)[edge1].x;
				color_next.y    = (*wc_colors)[edge1].y;
				color_next.z    = (*wc_colors)[edge1].z;
				color_next.w    = (*wc_colors)[edge1].w;
			}

			/* point.w is positive*/
			/* v.n - d = 0 */
			//float wInv_current = 1.0f / point_current.w;
			//float wInv_next = 1.0f / point_next.w;

			float dot0 =
			  point_current.x*plane.x +
			  point_current.y*plane.y +
			  point_current.z*plane.z +
			  point_current.w*plane.w;

			float dot1 =
			  point_next.x*plane.x +
			  point_next.y*plane.y +
			  point_next.z*plane.z +
			  point_next.w*plane.w;

			bool inside0 = dot0 > 0.0f;
			bool inside1 = dot1 > 0.0f;

			/* If start is inside, output it */
			if(inside0) {
				polygon.push_back(point_current);
				if(flags & SR_TEXCOORD0)
					polygon_tcoords0.push_back(tcoord0_current);
				if(flags & SR_TEXCOORD1)
					polygon_tcoords1.push_back(tcoord1_current);
				if(flags & SR_LIGHTING)
					polygon_normals.push_back(normal_current);
				if(flags & SR_COLOR)
					polygon_colors.push_back(color_current);
			}
			if(inside0 != inside1) {
				/* We're clipping an edge */
				float t=0.0f;
				float diff = 0.0f;
				/* swap vertices  if point0 is outside and point1 is inside*/
				if(inside0 == false) {
					//std::swap(wInv_current, wInv_next);
					std::swap(point_current, point_next);
					std::swap(dot0, dot1);

					//add all data cases
					if(flags & SR_TEXCOORD0)
						std::swap(tcoord0_current, tcoord0_next);
					if(flags & SR_TEXCOORD1)
						std::swap(tcoord1_current, tcoord1_next);
					if(flags & SR_LIGHTING)
						std::swap(normal_current, normal_next);
					if(flags & SR_COLOR)
						std::swap(color_current, color_next);
				}
				diff = dot0 - dot1;
				//const float eps = 1.0f / 4096.0f;

				if(std::abs(diff) > 1e-7f) //1e-7
					t = dot0 / diff;
				if(std::abs(t) < 0.001f) t = 0.0f;

				VectorPOD4f offsetPoint, clipPoint;
				offsetPoint.x = (point_next.x - point_current.x) * t;
				offsetPoint.y = (point_next.y - point_current.y) * t;
				offsetPoint.z = (point_next.z - point_current.z) * t;
				offsetPoint.w = (point_next.w - point_current.w) * t;
				clipPoint.x = point_current.x + offsetPoint.x;
				clipPoint.y = point_current.y + offsetPoint.y;
				clipPoint.z = point_current.z + offsetPoint.z;
				clipPoint.w = point_current.w + offsetPoint.w;
				polygon.push_back(clipPoint);

				if(flags & SR_TEXCOORD0) {
				    VectorPOD4f offsetTCoord, clipTCoord;
					offsetTCoord.x = (tcoord0_next.x - tcoord0_current.x) * t;
					offsetTCoord.y = (tcoord0_next.y - tcoord0_current.y) * t;
					offsetTCoord.z = (tcoord0_next.z - tcoord0_current.z) * t;
					offsetTCoord.w = (tcoord0_next.w - tcoord0_current.w) * t;
					clipTCoord.x = tcoord0_current.x + offsetTCoord.x;
					clipTCoord.y = tcoord0_current.y + offsetTCoord.y;
					clipTCoord.z = tcoord0_current.z + offsetTCoord.z;
					clipTCoord.w = tcoord0_current.w + offsetTCoord.w;
					polygon_tcoords0.push_back(clipTCoord);
				}
				if(flags & SR_TEXCOORD1) {
				    VectorPOD4f offsetTCoord, clipTCoord;
					offsetTCoord.x = (tcoord1_next.x - tcoord1_current.x) * t;
					offsetTCoord.y = (tcoord1_next.y - tcoord1_current.y) * t;
					offsetTCoord.z = (tcoord1_next.z - tcoord1_current.z) * t;
					offsetTCoord.w = (tcoord1_next.w - tcoord1_current.w) * t;
					clipTCoord.x = tcoord1_current.x + offsetTCoord.x;
					clipTCoord.y = tcoord1_current.y + offsetTCoord.y;
					clipTCoord.z = tcoord1_current.z + offsetTCoord.z;
					clipTCoord.w = tcoord1_current.w + offsetTCoord.w;
					polygon_tcoords1.push_back(clipTCoord);
				}
				if(flags & SR_LIGHTING) {
				    VectorPOD4f offsetNormal, clipNormal;
					offsetNormal.x = (normal_next.x - normal_current.x) * t;
					offsetNormal.y = (normal_next.y - normal_current.y) * t;
					offsetNormal.z = (normal_next.z - normal_current.z) * t;
					offsetNormal.w = (normal_next.w - normal_current.w) * t;
					clipNormal.x = normal_current.x + offsetNormal.x;
					clipNormal.y = normal_current.y + offsetNormal.y;
					clipNormal.z = normal_current.z + offsetNormal.z;
					clipNormal.w = normal_current.w + offsetNormal.w;
					polygon_normals.push_back(clipNormal);
				}
				if(flags & SR_COLOR) {
				    VectorPOD4f offsetColor, clipColor;
					offsetColor.x = (color_next.x - color_current.x) * t;
					offsetColor.y = (color_next.y - color_current.y) * t;
					offsetColor.z = (color_next.z - color_current.z) * t;
					offsetColor.w = (color_next.w - color_current.w) * t;
					clipColor.x = color_current.x + offsetColor.x;
					clipColor.y = color_current.y + offsetColor.y;
					clipColor.z = color_current.z + offsetColor.z;
					clipColor.w = color_current.w + offsetColor.w;
					polygon_colors.push_back(clipColor);
				}
			}
		}

		/* Split the resulting polygon into triangles.
		   TODO: Turn this into a triangle strip pattern to avoid narrow triangles.*/

		size_t vcount = polygon.size();

		if(vcount < 3)
			continue;

		for(int p=1; p<vcount-1; ++p) {

			polygon.push_back(polygon[0]);
			polygon.push_back(polygon[p]);
			polygon.push_back(polygon[p+1]);

			if(flags & SR_TEXCOORD0) {
				polygon_tcoords0.push_back(polygon_tcoords0[0]);
				polygon_tcoords0.push_back(polygon_tcoords0[p]);
				polygon_tcoords0.push_back(polygon_tcoords0[p+1]);
			}
			if(flags & SR_TEXCOORD1) {
				polygon_tcoords1.push_back(polygon_tcoords1[0]);
				polygon_tcoords1.push_back(polygon_tcoords1[p]);
				polygon_tcoords1.push_back(polygon_tcoords1[p+1]);
			}
			if(flags & SR_LIGHTING) {
				polygon_normals.push_back(polygon_normals[0]);
				polygon_normals.push_back(polygon_normals[p]);
				polygon_normals.push_back(polygon_normals[p+1]);
			}
			if(flags & SR_COLOR) {
				polygon_colors.push_back(polygon_colors[0]);
				polygon_colors.push_back(polygon_colors[p]);
				polygon_colors.push_back(polygon_colors[p+1]);
			}
		}

		wc_vertices->reserve(wc_vertices->size() + (polygon.size() - vcount));
		std::copy(polygon.begin() + vcount, polygon.end(), std::back_inserter(*wc_vertices));

		if(flags & SR_TEXCOORD0) {
			wc_tcoords0->reserve(wc_tcoords0->size() + (polygon_tcoords0.size() - vcount));
			std::copy(polygon_tcoords0.begin() + vcount, polygon_tcoords0.end(), std::back_inserter(*wc_tcoords0));
		}
		if(flags & SR_TEXCOORD1) {
			wc_tcoords1->reserve(wc_tcoords1->size() + (polygon_tcoords1.size() - vcount));
			std::copy(polygon_tcoords1.begin() + vcount, polygon_tcoords1.end(), std::back_inserter(*wc_tcoords1));
		}
		if(flags & SR_LIGHTING) {
			wc_normals->reserve(wc_normals->size() + (polygon_normals.size() - vcount));
			std::copy(polygon_normals.begin() + vcount, polygon_normals.end(), std::back_inserter(*wc_normals));
		}
		if(flags & SR_COLOR) {
			wc_colors->reserve(wc_colors->size() + (polygon_colors.size() - vcount));
			std::copy(polygon_colors.begin() + vcount, polygon_colors.end(), std::back_inserter(*wc_colors));
		}
	}
	wc_vertices->erase(wc_vertices->begin(), wc_vertices->begin() + vertexListSize);
	if(flags & SR_TEXCOORD0)
		wc_tcoords0->erase(wc_tcoords0->begin(), wc_tcoords0->begin() + vertexListSize);
	if(flags & SR_TEXCOORD1)
		wc_tcoords1->erase(wc_tcoords1->begin(), wc_tcoords1->begin() + vertexListSize);
	if(flags & SR_LIGHTING)
		wc_normals->erase(wc_normals->begin(), wc_normals->begin() + vertexListSize);
	if(flags & SR_COLOR)
		wc_colors->erase(wc_colors->begin(), wc_colors->begin() + vertexListSize);
#endif
}


