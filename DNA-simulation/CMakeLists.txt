#----------------------------------------------------------------------------
# Setup the project
cmake_minimum_required(VERSION 3.16...3.21)
project(rbe)

execute_process(COMMAND git rev-parse --short HEAD
  WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
  OUTPUT_VARIABLE GIT_HASH OUTPUT_STRIP_TRAILING_WHITESPACE ERROR_QUIET)
if(NOT GIT_HASH)
  set(GIT_HASH "unversioned")
endif()
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/generated)
configure_file(cmake/git_version.cpp.in generated/git_version.cpp @ONLY)
configure_file(cmake/git_version.hh generated/git_version.hh COPYONLY)
add_library(git_version ${CMAKE_CURRENT_BINARY_DIR}/generated/git_version.cpp)
target_include_directories(git_version PUBLIC ${CMAKE_CURRENT_BINARY_DIR}/generated)

#----------------------------------------------------------------------------
# Find Geant4 package, activating all available Vis drivers by default
# You can set WITH_GEANT4_VIS to OFF via the command line or ccmake/cmake-gui
# to build a batch mode only executable
#
option(WITH_GEANT4_VIS "Build example with Geant4 Vis drivers" ON)
if(WITH_GEANT4_VIS)
  find_package(Geant4 11.1.3 EXACT REQUIRED ui_all vis_all)
else()
  find_package(Geant4 11.1.3 EXACT REQUIRED)
endif()

#----------------------------------------------------------------------------
# Setup Geant4 include directories and compile definitions
#
include(${Geant4_USE_FILE})

#----------------------------------------------------------------------------
#----------------------------------------------------------------------------
# Locate sources and headers for this project
#

include_directories(${PROJECT_SOURCE_DIR}/include
                    ${Geant4_INCLUDE_DIR})
file(GLOB sources ${PROJECT_SOURCE_DIR}/src/*.cc)
file(GLOB headers ${PROJECT_SOURCE_DIR}/include/*.hh)

#----------------------------------------------------------------------------
# Add the executable, and link it to the Geant4 libraries
#
add_executable(rbe rbe.cc ${sources} ${headers})
target_link_libraries(rbe ${Geant4_LIBRARIES} git_version)

#----------------------------------------------------------------------------
# Copy all scripts to the build directory, i.e. the directory in which we
# build rbe. This is so that we can run the executable directly because it
# relies on these scripts being in the current working directory.
#

file(GLOB MAC_FILES ${PROJECT_SOURCE_DIR}/*.mac 
     ${PROJECT_SOURCE_DIR}/*.in
     ${PROJECT_SOURCE_DIR}/*.py
     ${PROJECT_SOURCE_DIR}/*.C)

foreach(_script ${MAC_FILES})
  configure_file(
    ${_script}
    ${PROJECT_BINARY_DIR}/.
    COPYONLY
    )
endforeach()

#----------------------------------------------------------------------------
# Install the executable to 'bin' directory under CMAKE_INSTALL_PREFIX
#
install(TARGETS rbe DESTINATION bin)



target_compile_features(rbe PRIVATE cxx_std_17)
