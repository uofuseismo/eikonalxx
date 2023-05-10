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
find_package(spdlog REQUIRED)
set(FindTBB_DIR ${CMAKE_SOURCE_DIR}/cmake)
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
    src/station2d.cpp src/station3d.cpp
    src/analytic/homogeneous2d.cpp src/analytic/homogeneous3d.cpp
    src/analytic/linearGradient2d.cpp
    src/io/vtkRectilinearGrid2d.cpp src/io/vtkRectilinearGrid3d.cpp)

set(BUILD_SHARED_LIBS YES)
add_library(eikonalxx ${LIBSRC})
target_link_libraries(eikonalxx PRIVATE sycl ${TBB_LIBRARY} stdc++fs)
if (${Boost_FOUND})
   target_compile_definitions(eikonalxx PRIVATE USE_BOOST)
   target_include_directories(eikonalxx PRIVATE ${Boost_INCLUDE_DIRS})
endif()
set_target_properties(eikonalxx PROPERTIES
                      CXX_STANDARD 20
                      CXX_STANDARD_REQUIRED YES
                      CXX_EXTENSIONS NO)


if (WRAP_PYTHON)
   find_package(pybind11 REQUIRED)
   add_library(pyEikonalXX MODULE
               python/pyEikonalXX.cpp
               python/pySolverOptions.cpp
               python/pyGeometry2d.cpp python/pyGeometry3d.cpp
               python/pyModel2d.cpp
               python/pySource2d.cpp)
   target_link_libraries(pyEikonalXX PRIVATE
                         pybind11::module eikonalxx)
   set_target_properties(pyEikonalXX PROPERTIES
                         PREFIX ""
                         CXX_STANDARD 20
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
               testing/station.cpp
               testing/model.cpp
               testing/solverOptions.cpp
               testing/graph.cpp
               testing/analytic.cpp
               testing/solver2d.cpp
               testing/solver3d.cpp
               )
target_link_libraries(unitTests PRIVATE eikonalxx ${GTEST_BOTH_LIBRARIES} ${TBB_LIBRARY} OpenCL pthread)
target_include_directories(unitTests PRIVATE ${GTEST_INCLUDE_DIRS} ${CMAKE_SOURCE_DIR}/include)
set_target_properties(unitTests PROPERTIES
                      CXX_STANDARD 20 
                      CXX_STANDARD_REQUIRED YES
                      CXX_EXTENSIONS NO)
add_test(NAME unitTest
         COMMAND unitTests)

add_executable(converge2d
               testing/converge2d.cpp)
target_link_libraries(converge2d PRIVATE eikonalxx)
target_include_directories(converge2d PRIVATE ${CMAKE_SOURCE_DIR}/include)
set_target_properties(converge2d PROPERTIES
                      CXX_STANDARD 20
                      CXX_STANDARD_REQUIRED YES
                      CXX_EXTENSIONS NO)

##########################################################################################
#                                         Installation                                   #
##########################################################################################
include(GNUInstallDirs)
include(CMakePackageConfigHelpers)
configure_package_config_file(
    cmake/EikonalXXConfig.cmake.in "${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}Config.cmake"
    INSTALL_DESTINATION "${CMAKE_INSTALL_LIBDIR}/cmake/${PROJECT_NAME}")
write_basic_package_version_file(
    "${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}ConfigVersion.cmake"
    VERSION "${version}"
    COMPATIBILITY AnyNewerVersion
)
if (WRAP_PYTHON)
   install(TARGETS eikonalxx pyEikonalXX
           EXPORT ${PROJECT_NAME}-targets
           RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
           LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
           ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
           PUBLIC_HEADER DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
           COMPONENT Runtime)
else()
   install(TARGETS eikonalxx
           EXPORT ${PROJECT_NAME}-targets
           RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
           LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
           ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
           PUBLIC_HEADER DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
           COMPONENT Runtime)
endif()
install(DIRECTORY ${CMAKE_SOURCE_DIR}/include/eikonalxx
        DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
export(EXPORT ${PROJECT_NAME}-targets
       FILE "${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}Targets.cmake")
install(EXPORT ${PROJECT_NAME}-targets
        FILE ${PROJECT_NAME}Targets.cmake
        DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/${PROJECT_NAME})
##########################################################################################
#                                     CPACK Packaging                                    #
##########################################################################################
set(CPACK_PACKAGE_NAME ${PROJECT_NAME})
set(CPACK_PACKAGE_VENDOR "University of Utah Seismograph Stations")
set(CPACK_PACKAGE_CONTACT "ben.baker@utah.edu")
set(CPACK_PACKAGE_LICENSE "MIT")
set(CPACK_PACKAGE_DESCRIPTION_SUMMARY "A library for solving the eikonal equation.")
set(CPACK_PACKAGE_INSTALL_DIRECTORY ${CPACK_PACKAGE_NAME})
set(CPACK_VERBATIM_VARIABLES YES)
set(CPACK_PACKAGE_VERSION_MAJOR ${PROJECT_VERSION_MAJOR})
set(CPACK_PACKAGE_VERSION_MINOR ${PROJECT_VERSION_MINOR})
set(CPACK_PACKAGE_VERSION_PATCH ${PROJECT_VERSION_PATCH})
set(CPACK_RESOURCE_FILE_LICENSE "${CMAKE_CURRENT_SOURCE_DIR}/LICENSE")
set(CPACK_COMPONENTS_ALL libraries headers)
if (WIN32)
    set(CPACK_GENERATOR ZIP WIX)
elseif (APPLE)
    set(CPACK_GENERATOR TGZ productbuild)
elseif(CMAKE_SYSTEM_NAME STREQUAL "Linux")
    set(CPACK_GENERATOR TGZ RPM)
else()
    set(CPACK_GENERATOR TGZ)
endif()
set(CPACK_SOURCE_IGNORE_FILES
  /\\.git/
  \\.swp
  \\.orig
  /CMakeLists\\.txt\\.user
  /private/
)
include(CPack) # Put this last!
