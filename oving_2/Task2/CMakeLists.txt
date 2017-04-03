
cmake_minimum_required (VERSION 2.6)
project (Task2)
find_package(Armadillo REQUIRED)
link_directories("/usr/lib64")

include_directories(${ARMADILLO_INCLUDE_DIRS})
set(LIB ${ARMADILLO_LIBRARIES})

find_library( LAPACK lapack )
if ( NOT LAPACK )
  MESSAGE( FATAL_ERROR "Did not find lapack!")
endif()
MESSAGE( STATUS "Found lapack ${LAPACK}")

set( LIB ${LIB} ${LAPACK} )

find_package( GSL REQUIRED )
if( NOT GSL )
  MESSAGE( FATAL_ERROR "Did not find GSL!")
endif()
MESSAGE( STATUS "Found GSL ${GSL}")
set( LIB ${LIB} ${GSL_LIBRARIES} )

find_package(OpenMP)
if (OPENMP_FOUND)
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif()

add_executable(main main.cpp Lattice.cpp)
target_link_libraries(main ${LIB})

add_executable(buffalo BuffaloMain.cpp Buffalo.cpp)
target_link_libraries(buffalo ${LIB})
