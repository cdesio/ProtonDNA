cmake_minimum_required(VERSION 3.16...3.27)
project(TrackingGrid)

option(WITH_GEANT4_UIVIS "Build with Geant4 UI and visualization drivers" ON)

if(WITH_GEANT4_UIVIS)
  find_package(Geant4 11.1.3 EXACT REQUIRED ui_all vis_all)
else()
  find_package(Geant4 11.1.3 EXACT REQUIRED)
endif()

include(${Geant4_USE_FILE})

file(GLOB sources ${PROJECT_SOURCE_DIR}/src/*.cc)
file(GLOB headers ${PROJECT_SOURCE_DIR}/include/*.hh)

add_executable(TrackingGrid
  TrackingGrid.cc
  ${sources}
  ${headers}
)

target_include_directories(TrackingGrid PRIVATE
  ${PROJECT_SOURCE_DIR}/include
)
if(WITH_GEANT4_UIVIS)
  target_compile_definitions(TrackingGrid PRIVATE TRACKINGGRID_WITH_GEANT4_UIVIS=1)
endif()
target_link_libraries(TrackingGrid ${Geant4_LIBRARIES})
target_compile_features(TrackingGrid PRIVATE cxx_std_17)

file(GLOB macro_files ${PROJECT_SOURCE_DIR}/macros/*.mac)
foreach(_script ${macro_files})
  configure_file(${_script} ${PROJECT_BINARY_DIR}/. COPYONLY)
endforeach()

install(TARGETS TrackingGrid DESTINATION bin)
