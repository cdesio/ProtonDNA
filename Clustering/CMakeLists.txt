cmake_minimum_required (VERSION 3.16...3.27)
project(pyClustering)



include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/src)

find_package(Python3 REQUIRED COMPONENTS Interpreter Development)
execute_process(
  COMMAND "${Python3_EXECUTABLE}" -m pybind11 --cmakedir
  OUTPUT_VARIABLE pybind11_DIR
  OUTPUT_STRIP_TRAILING_WHITESPACE
)
find_package(pybind11 REQUIRED CONFIG)



file(GLOB sources ${PROJECT_SOURCE_DIR}/src/*.cc)
file(GLOB headers ${PROJECT_SOURCE_DIR}/include/*.hh)

file(GLOB MAC_FILES 
${PROJECT_SOURCE_DIR}/runSingle.py
${PROJECT_SOURCE_DIR}/sugarPosTest.csv
${PROJECT_SOURCE_DIR}/test.csv
${PROJECT_SOURCE_DIR}/testFile.root
${PROJECT_SOURCE_DIR}/__init__.py
)

foreach(_script ${MAC_FILES})
  configure_file(
    ${_script}
    ${PROJECT_BINARY_DIR}/.
    COPYONLY
    )
endforeach()

pybind11_add_module(clustering pyClustering.cc  ${sources} ${headers})
