cmake_minimum_required(VERSION 3.15)
project(algos_cpp)

set(CMAKE_CXX_STANDARD 11)

set(BIN_DESTINATION ${PROJECT_SOURCE_DIR}/bin)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${BIN_DESTINATION})
#set(CMAKE_CURRENT_SOURCE_DIR ${PROJECT_SOURCE_DIR}/src)
set(CMAKE_CXX_FLAGS "-Wall")

file(GLOB_RECURSE algos RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}" "${CMAKE_CURRENT_SOURCE_DIR}/__*.cpp")
foreach(algo IN LISTS algos)
    # Get file name without directory
    get_filename_component(algo ${algo} NAME_WE)
    add_executable(${algo} ${algo})
endforeach()

#add_executable(algos_cpp main.cpp __FastFourierTransform.cpp)
