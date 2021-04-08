# required cmake version
cmake_minimum_required(VERSION 3.10)
project(EikonalXX VERSION 0.1.0 LANGUAGES CXX)
enable_testing()

include(CheckCXXCompilerFlag)

include(FindOpenMP QUIET)
if (OPENMP_FOUND)
   message("OpenMP found")
   string(APPEND CMAKE_CXX_FLAGS " ${OpenMP_CXX_FLAGS}")
endif()

find_package(GTest REQUIRED)
find_package(Boost COMPONENTS graph)
set(FindTBB_DIR ${CMAKE_SOURCE_DIR}/CMakeModules)
find_package(FindTBB REQUIRED)

configure_file(${CMAKE_SOURCE_DIR}/include/private/version.hpp.in
               ${CMAKE_SOURCE_DIR}/include/eikonalxx/version.hpp)

include_directories(
   ${CMAKE_SOURCE_DIR}/include
   ${TBB_INCLUDE_DIR}
)

set(PUBLIC_HEADER_DIRECTORIES
    ${CMAKE_SOURCE_DIR}/include/eikonalxx)
set(PRIVATE_HEADER_DIRECTORIES
    ${CMAKE_SOURCE_DIR}/include/private)

set(LIBSRC
    src/version.cpp
    src/solverOptions.cpp
    src/geometry2d.cpp src/geometry3d.cpp
    src/graph3d.cpp
    src/model2d.cpp src/model3d.cpp
    src/source2d.cpp
    src/source3d.cpp
    src/solver2d.cpp
    src/solver3d.cpp
    src/station3d.cpp)

set(BUILD_SHARED_LIBS YES)
add_library(eikonalxx ${LIBSRC})
target_link_libraries(eikonalxx PRIVATE sycl ${TBB_LIBRARY})
if (${Boost_FOUND})
   target_compile_definitions(eikonalxx PRIVATE USE_BOOST)
   target_include_directories(eikonalxx PRIVATE ${Boost_INCLUDE_DIRS})
endif()
set_target_properties(eikonalxx PROPERTIES
                      CXX_STANDARD 17
                      CXX_STANDARD_REQUIRED YES
                      CXX_EXTENSIONS NO)


if (WRAP_PYTHON)
   find_package(pybind11 REQUIRED)
   add_library(pyEikonalXX MODULE
               python/pyEikonalXX.cpp
               python/pyGeometry2d.cpp python/pyGeometry3d.cpp)
   target_link_libraries(pyEikonalXX PRIVATE
                         pybind11::module eikonalxx)
   set_target_properties(pyEikonalXX PROPERTIES
                         PREFIX ""
                         CXX_STANDARD 17
                         CXX_STANDARD_REQUIRED YES
                         CXX_EXTENSIONS NO)
   target_include_directories(pyEikonalXX PRIVATE pybind11
                              ${CMAKE_SOURCE_DIR}/python)
endif()

##########################################################################################
#                                        Testing                                         #
##########################################################################################
add_executable(unitTests
               testing/main.cpp
               testing/geometry.cpp
               testing/source.cpp
               testing/model.cpp
               testing/solverOptions.cpp
               testing/graph.cpp
               testing/solver2d.cpp
               testing/solver3d.cpp
               )
target_link_libraries(unitTests PRIVATE eikonalxx ${GTEST_BOTH_LIBRARIES} ${TBB_LIBRARY} OpenCL pthread)
target_include_directories(unitTests PRIVATE ${GTEST_INCLUDE_DIRS})
set_target_properties(unitTests PROPERTIES
                      CXX_STANDARD_REQUIRED YES
                      CXX_EXTENSIONS NO)
add_test(NAME unitTest
         COMMAND unitTests)
