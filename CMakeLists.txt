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

include_directories(
   ${CMAKE_SOURCE_DIR}/include
)

set(PUBLIC_HEADER_DIRECTORIES
    ${CMAKE_SOURCE_DIR}/include/eikonalxx)
set(PRIVATE_HEADER_DIRECTORIES
    ${CMAKE_SOURCE_DIR}/include/private)

set(LIBSRC
    src/solverOptions.cpp
    src/geometry2d.cpp src/geometry3d.cpp
    src/graph3d.cpp
    src/model2d.cpp src/model3d.cpp
    src/solver2d.cpp)

set(BUILD_SHARED_LIBS YES)
add_library(eikonalxx ${LIBSRC})
target_link_libraries(eikonalxx PRIVATE sycl)# ${TBB_LIBRARY})
if (${Boost_FOUND})
   target_compile_definitions(eikonalxx PRIVATE USE_BOOST)
   target_include_directories(eikonalxx PRIVATE ${Boost_INCLUDE_DIRS})
endif()
set_target_properties(eikonalxx PROPERTIES
                      CXX_STANDARD 17
                      CXX_STANDARD_REQUIRED YES
                      CXX_EXTENSIONS NO)

add_executable(unitTests
               testing/main.cpp
               testing/geometry.cpp
               testing/model.cpp
               testing/solverOptions.cpp
               testing/graph.cpp
               testing/solver2d.cpp
               testing/solver3d.cpp
               )
target_link_libraries(unitTests PRIVATE eikonalxx ${GTEST_BOTH_LIBRARIES} OpenCL pthread)
target_include_directories(unitTests PRIVATE ${GTEST_INCLUDE_DIRS})
set_target_properties(unitTests PROPERTIES
                      CXX_STANDARD_REQUIRED YES
                      CXX_EXTENSIONS NO)

add_test(NAME unitTest
         COMMAND unitTests)
