cmake_minimum_required (VERSION 2.6)
project (PerpleX Fortran)

# The version number.
set (VERSION 6.8.3)

set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fPIC -O3")

include_directories(include)

add_library(perplex STATIC nlib.f olib.f tlib.f rlib.f flib.f clib.f dlib.f dumlib.f getxz1.f resub.f pslib.f pscom.f cont_lib.f)
set_target_properties(perplex PROPERTIES VERSION ${VERSION})

set(TARGET_SRC "actcor" "build" "fluids" "ctransf" "frendly" "htog" "meemum"
  "pstable" "pspts" "psvdraw" "psect" "pt2curv" "vertex" "werami" "unsplt" "convex")

FOREACH(_binary ${TARGET_SRC})
  add_executable(${_binary} ${_binary}.f)
  target_link_libraries(${_binary} perplex)
ENDFOREACH()

# Static library. Note no dlib (redefines getxz)
add_library(perplex_o STATIC nlib.f olib.f tlib.f rlib.f flib.f clib.f dumlib.f getxz1.f resub.f pslib.f pscom.f cont_lib.f)
set_target_properties(perplex_o PROPERTIES VERSION ${VERSION}
  CXX_VISIBILITY_PRESET hidden
  Fortran_VISIBILITY_PRESET hidden)


# Wrapper shared library
project (PerpleX CXX)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wl,--version-script=perplex.version") 

set(SHARED_LIB perplex_cpp)

add_library(${SHARED_LIB} SHARED perplex.cpp)
set_target_properties(${SHARED_LIB} PROPERTIES
  VERSION ${VERSION})
target_link_libraries(${SHARED_LIB} PRIVATE perplex_o)

set(INCLUDE_INSTALL_DIR ${CMAKE_INSTALL_PREFIX}/include)
set(LIB_INSTALL_DIR ${CMAKE_INSTALL_PREFIX}/lib/)

include(CMakePackageConfigHelpers)
configure_package_config_file(PerpleXConfig.cmake.in ${CMAKE_CURRENT_BINARY_DIR}/PerpleXConfig.cmake
  INSTALL_DESTINATION ${LIB_INSTALL_DIR}/PerpleX/cmake
  PATH_VARS INCLUDE_INSTALL_DIR )

write_basic_package_version_file(${CMAKE_CURRENT_BINARY_DIR}/PerpleXConfigVersion.cmake
  VERSION ${VERSION}
  COMPATIBILITY ExactVersion )
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/PerpleXConfig.cmake ${CMAKE_CURRENT_BINARY_DIR}/PerpleXConfigVersion.cmake
  DESTINATION ${LIB_INSTALL_DIR}/PerpleX/cmake )


INSTALL(TARGETS ${TARGET_SRC} RUNTIME DESTINATION bin COMPONENT RUNTIME)
INSTALL(TARGETS perplex_cpp LIBRARY DESTINATION lib COMPONENT RUNTIME)

# headers:
INSTALL(DIRECTORY include/ 
  DESTINATION include
  COMPONENT includes
  FILES_MATCHING PATTERN "*.h")
