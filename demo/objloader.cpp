#include "objloader.h"
#include "myassert.h"
#include <fstream>
#include <sstream>

static inline bool str_starts_with(const std::string& token, const std::string& str)
{
  size_t pos = str.find(token);
  return (pos == 0);
}

static bool readVCoord(const std::string& line, VectorPOD4f& vout)
{
  std::istringstream strm;

  //skip 'v' and ' '
  strm.str(line.substr(2));
  std::vector<float> fv;
  float f;
  while((strm >> f)){
	fv.push_back(f);
  }
  if(fv.size() != 3){
	ASSERT(false);
	return false;
  }
  //Blender uses +z = up, +y = forwards
  //So correct the coordinate system to use +y = up, -z = forwards
  vout.x = -fv[0];
  vout.y = -fv[1];
  vout.z = fv[2];
  vout.w = 1.0f;
  return true;
}

static bool readTCoord(const std::string& line, VectorPOD4f& vout)
{
  std::istringstream strm;

  //skip 'vt' and ' '
  strm.str(line.substr(3));
  std::vector<float> fv;
  float f;
  while((strm >> f)){
	fv.push_back(f);
  }
  if(fv.size() != 2){
	ASSERT(false);
	return false;
  }
  vout.x = fv[0];
  vout.y = 1.0f - fv[1];
  vout.z = 0.0f;
  vout.w = 1.0f;
  return true;
}

static bool readIndices(const std::string& line,
						int& v1, int& v2, int& v3,
						int& t1, int& t2, int& t3)
{
  std::istringstream strm;

  //skip 'f' and ' '
  strm.str(line.substr(2));
  std::vector<int> iv;
  int idx;
  while((strm >> idx)){
	iv.push_back(idx);
	strm.ignore(); //ignore ' ' and '/'
  }
  if(iv.size() != 6){
	ASSERT(false);
	return false;
  }
  //subtract 1 to get zero-based indices
  //Wavefron OBJ counts from 1
  v1 = iv[0] - 1;
  v2 = iv[2] - 1;
  v3 = iv[4] - 1;
  t1 = iv[1] - 1;
  t2 = iv[3] - 1;
  t3 = iv[5] - 1;
  return true;
}
 
bool importWaveFrontObjModel(const std::string& path,
							 std::vector<VectorPOD4f>& vertexData,
							 std::vector<VectorPOD4f>& tcoordData)
{
  std::ifstream f(path.c_str());
  if(!f.is_open()) return false;

  std::string line;
  VectorPOD4f vout, tout;
  std::vector<VectorPOD4f> vertexList;
  std::vector<VectorPOD4f> tcoordList;
  std::vector<int> indexList;

  while(std::getline(f, line)){
	if(str_starts_with("v ", line)){
	  if(!readVCoord(line, vout)) return false;
	  vertexList.push_back(vout);
	} else if(str_starts_with("vt ", line)){
	  if(!readTCoord(line, tout)) return false;
	  tcoordList.push_back(tout);
	} else if(str_starts_with("f ", line)){
	  int v1,v2,v3,t1,t2,t3;
	  if(!readIndices(line, v1, v2, v3, t1, t2, t3))
		return false;
	  indexList.push_back(v1);
	  indexList.push_back(v2);
	  indexList.push_back(v3);
	  indexList.push_back(t1);
	  indexList.push_back(t2);
	  indexList.push_back(t3);
	}
  }
  if((indexList.size() % 6) != 0){
	ASSERT(false);
	return false;
  }

  for(int i = 0; i < indexList.size(); i += 6){
	vertexData.push_back(vertexList[indexList[i + 0]]);
	vertexData.push_back(vertexList[indexList[i + 1]]);
	vertexData.push_back(vertexList[indexList[i + 2]]);
	tcoordData.push_back(tcoordList[indexList[i + 3]]);
	tcoordData.push_back(tcoordList[indexList[i + 4]]);
	tcoordData.push_back(tcoordList[indexList[i + 5]]);
  }
  return true;
}
