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
  unsigned int edge0, edge1;
  unsigned int tri;

  size_t vertexListSize = wc_vertices.size();

  for(tri=0; tri<vertexListSize; tri+=3){
	std::vector<Vector4f> polygon; /* outputted vertices */
	std::vector<Vector4f> polygon_tcoords0; /* interpolants */
	std::vector<Vector4f> polygon_tcoords1;
	std::vector<Vector4f> polygon_normals;
	std::vector<Vector4f> polygon_colors;

	for(edge0=tri+2, edge1=tri; edge1 < tri+3; edge0 = edge1++){
	  Vector4f point_current, point_next, tcoord0_current,
		tcoord0_next, tcoord1_current, tcoord1_next, normal_current,
		normal_next, color_current, color_next;

	  point_current  = wc_vertices[edge0];
	  point_next     = wc_vertices[edge1];

#if 1
	  if(flags & SR_TEXCOORD0) {
		tcoord0_current = wc_tcoords0[edge0];
		tcoord0_next    = wc_tcoords0[edge1];
	  }
	  if(flags & SR_TEXCOORD1){
		tcoord1_current = wc_tcoords1[edge0];
		tcoord1_next    = wc_tcoords1[edge1];
	  }
	  if(flags & SR_LIGHTING){
		normal_current = wc_normals[edge0];
		normal_next    = wc_normals[edge1];
	  }
	  if(flags & SR_COLOR){
		color_current = wc_colors[edge0];
		color_next    = wc_colors[edge1];
	  }
#endif

	  /* point.w is positive*/
	  /* v.n - d = 0 */
	  //float wInv_current = 1.0f / point_current.w;
	  //float wInv_next = 1.0f / point_next.w;

	  float dot0 = dot(point_current, plane) + point_current.w * plane.w;
	  float dot1 = dot(point_next,    plane) + point_next.w    * plane.w;
		
	  bool inside0 = dot0 > 0.0f;
	  bool inside1 = dot1 > 0.0f;

	  /* If start is inside, output it */
	  if(inside0){
		polygon.push_back(point_current);
#if 1
		if(flags & SR_TEXCOORD0)
		  polygon_tcoords0.push_back(tcoord0_current);
		if(flags & SR_TEXCOORD1)
		  polygon_tcoords1.push_back(tcoord1_current);
		if(flags & SR_LIGHTING)
		  polygon_normals.push_back(normal_current);
		if(flags & SR_COLOR)
		  polygon_colors.push_back(color_current);
#endif

	  }
	  if(inside0 != inside1){
		/* We're clipping an edge */
		float t=0.0f;
		float diff = 0.0f;
		/* swap vertices  if point0 is outside and point1 is inside*/
		if(inside0 == false){
		  //std::swap(wInv_current, wInv_next);
		  std::swap(point_current, point_next);
		  std::swap(dot0, dot1);

		  //add all data cases
#if 1
		  if(flags & SR_TEXCOORD0)
			std::swap(tcoord0_current, tcoord0_next);
		  if(flags & SR_TEXCOORD1)
			std::swap(tcoord1_current, tcoord1_next);
		  if(flags & SR_LIGHTING)
			std::swap(normal_current, normal_next);
		  if(flags & SR_COLOR)
			std::swap(color_current, color_next);
#endif

		}
		diff = dot0 - dot1;
		//const float eps = 1.0f / 4096.0f;

		if(std::abs(diff) > 1e-7f) //1e-7
		  t = dot0 / diff;
		if(std::abs(t) < 0.001f) t = 0.0f;

		Vector4f offsetPoint = (point_next - point_current) * t;
		Vector4f clipPoint = point_current + offsetPoint;
		clipPoint.w = point_current.w + (point_next.w - point_current.w) * t;
		polygon.push_back(clipPoint);

#if 1
		if(flags & SR_TEXCOORD0){
		  Vector4f offsetTCoord = (tcoord0_next - tcoord0_current) * t;
		  Vector4f clipTCoord = tcoord0_current + offsetTCoord;
		  polygon_tcoords0.push_back(clipTCoord);
		}
		if(flags & SR_TEXCOORD1){
		  Vector4f offsetTCoord = (tcoord1_next - tcoord1_current) * t;
		  Vector4f clipTCoord = tcoord1_current + offsetTCoord;
		  polygon_tcoords1.push_back(clipTCoord);
		}
		if(flags & SR_LIGHTING){
		  Vector4f offsetNormal = (normal_next - normal_current) * t;
		  Vector4f clipNormal = normal_current + offsetNormal;
		  polygon_normals.push_back(clipNormal);
		}
		if(flags & SR_COLOR){
		  Vector4f offsetColor = (color_next - color_current) * t;
		  Vector4f clipColor = color_current + offsetColor;
		  polygon_colors.push_back(clipColor);
		}
#endif

	  }
	}

	/* Split the resulting polygon into triangles.
	   TODO: Turn this into a triangle strip pattern to avoid narrow triangles.*/

	size_t vcount = polygon.size();

	if(vcount < 3)
	  continue;

	for(int p=1; p<vcount-1; ++p){

	  polygon.push_back(polygon[0]);
	  polygon.push_back(polygon[p]);
	  polygon.push_back(polygon[p+1]);
#if 1
	  if(flags & SR_TEXCOORD0){
		polygon_tcoords0.push_back(polygon_tcoords0[0]);
		polygon_tcoords0.push_back(polygon_tcoords0[p]);
		polygon_tcoords0.push_back(polygon_tcoords0[p+1]);
	  }
	  if(flags & SR_TEXCOORD1){
		polygon_tcoords1.push_back(polygon_tcoords1[0]);
		polygon_tcoords1.push_back(polygon_tcoords1[p]);
		polygon_tcoords1.push_back(polygon_tcoords1[p+1]);
	  }
	  if(flags & SR_LIGHTING){
		polygon_normals.push_back(polygon_normals[0]);
		polygon_normals.push_back(polygon_normals[p]);
		polygon_normals.push_back(polygon_normals[p+1]);
	  }
	  if(flags & SR_COLOR){
		polygon_colors.push_back(polygon_colors[0]);
		polygon_colors.push_back(polygon_colors[p]);
		polygon_colors.push_back(polygon_colors[p+1]);
	  }
#endif
	}

	wc_vertices.reserve(wc_vertices.size() + (polygon.size() - vcount));
	std::copy(polygon.begin() + vcount, polygon.end(), std::back_inserter(wc_vertices));
#if 1
	if(flags & SR_TEXCOORD0){
	  wc_tcoords0.reserve(wc_tcoords0.size() + (polygon_tcoords0.size() - vcount));
	  std::copy(polygon_tcoords0.begin() + vcount, polygon_tcoords0.end(), std::back_inserter(wc_tcoords0));
	}
	if(flags & SR_TEXCOORD1){
	  wc_tcoords1.reserve(wc_tcoords1.size() + (polygon_tcoords1.size() - vcount));
	  std::copy(polygon_tcoords1.begin() + vcount, polygon_tcoords1.end(), std::back_inserter(wc_tcoords1));
	}
	if(flags & SR_LIGHTING){
	  wc_normals.reserve(wc_normals.size() + (polygon_normals.size() - vcount));
	  std::copy(polygon_normals.begin() + vcount, polygon_normals.end(), std::back_inserter(wc_normals));
	}
	if(flags & SR_COLOR){
	  wc_colors.reserve(wc_colors.size() + (polygon_colors.size() - vcount));
	  std::copy(polygon_colors.begin() + vcount, polygon_colors.end(), std::back_inserter(wc_colors));
	}
#endif
  }

  wc_vertices.erase(wc_vertices.begin(), wc_vertices.begin() + vertexListSize);
#if 1
  if(flags & SR_TEXCOORD0)
	wc_tcoords0.erase(wc_tcoords0.begin(), wc_tcoords0.begin() + vertexListSize);
  if(flags & SR_TEXCOORD1)
	wc_tcoords1.erase(wc_tcoords1.begin(), wc_tcoords1.begin() + vertexListSize);
  if(flags & SR_LIGHTING)
	wc_normals.erase(wc_normals.begin(), wc_normals.begin() + vertexListSize);
  if(flags & SR_COLOR)
	wc_colors.erase(wc_colors.begin(), wc_colors.begin() + vertexListSize);
#endif
}
