IF(DEFINED CMAKE_BUILD_TYPE)
    SET(CMAKE_BUILD_TYPE ${CMAKE_BUILD_TYPE} CACHE STRING "Choose the type of build, options are: None(CMAKE_CXX_FLAGS or CMAKE_C_FLAGS used) Debug Release RelWithDebInfo MinSizeRel.")
ELSE()
    SET(CMAKE_BUILD_TYPE Release CACHE STRING "Choose the type of build, options are: None(CMAKE_CXX_FLAGS or CMAKE_C_FLAGS used) Debug Release RelWithDebInfo MinSizeRel.")
ENDIF()

cmake_minimum_required(VERSION 2.6)

project(lineage)

set(CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake/modules")

include_directories(include)

file(GLOB headers
    cmake/modules/*
    include/levinkov/*.hxx
    include/andres/*.hxx
    include/andres/graphics/*.hxx
    include/andres/graph/*.hxx
    include/andres/ilp/*.hxx
    include/andres/hdf5/*.hxx
    include/lineage/*.hxx
)
enable_testing()

##############################################################################
# GNU g++ specific settings
##############################################################################
if(CMAKE_COMPILER_IS_GNUCXX)
    add_definitions(-std=c++11)
endif()

##############################################################################
# OpenMP
##############################################################################
find_package(OpenMP)
if(OPENMP_FOUND)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    message(STATUS "OpenMP found")
else()
    message("OpenMP not found")
endif()

##############################################################################
# HDF5
##############################################################################
find_package(HDF5)
if(HDF5_FOUND)
    include_directories(${HDF5_INCLUDE_DIR})
endif()

##############################################################################
# GLUT
##############################################################################
find_package(GLUT)

##############################################################################
# GUROBI
##############################################################################
find_package(GUROBI)
if(GUROBI_FOUND)
    message(STATUS "Found Gurobi")
    include_directories(${GUROBI_INCLUDE_DIR})
endif()

##############################################################################
# OpenMP
##############################################################################
# find_package(OpenMP)
if(OPENMP_FOUND)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif()

##############################################################################
# targets: graphics
##############################################################################
if(GLUT_FOUND AND HDF5_FOUND)
    add_executable(viewer-opengl src/andres/graphics/viewer-opengl.cxx ${headers})
    target_link_libraries(viewer-opengl ${GLUT_LIBRARIES} GL ${HDF5_LIBRARIES})
endif()

##############################################################################
# targets: lineage
##############################################################################
if(GUROBI_FOUND)
    add_executable(track-ilp src/lineage/track-ilp.cxx)
    target_link_libraries(track-ilp ${GUROBI_LIBRARIES})

    add_executable(track-lp src/lineage/track-lp.cxx)
    target_link_libraries(track-lp ${GUROBI_LIBRARIES})
endif()

if(HDF5_FOUND)
    add_executable(draw src/lineage/draw.cxx)
    target_link_libraries(draw ${HDF5_LIBRARIES})
endif()

add_executable(plot-intermediate src/lineage/plot-intermediate.cxx)

add_executable(validate-solution src/lineage/validate-solution.cxx)
