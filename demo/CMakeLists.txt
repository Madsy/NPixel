PROJECT(demo)
SET( ${PROJECT_NAME}_SOURCES
  main.cpp
  clipplane.cpp
  meshgen.cpp
  rasterizer_new.cpp
  texture.cpp
  framebuffer.cpp
  reci.cpp
  vertexdata.cpp
  loop.cpp
  misc.cpp
)

ADD_EXECUTABLE(${PROJECT_NAME} ${${PROJECT_NAME}_SOURCES})
TARGET_LINK_LIBRARIES( ${PROJECT_NAME} ${SDL_LIBRARY} ${IL_LIBRARIES} ${ILU_LIBRARIES} ${ZLIB_LIBRARIES})
