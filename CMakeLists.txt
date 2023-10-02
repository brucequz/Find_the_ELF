cmake_minimum_required(VERSION 3.22.1)
project(elf VERSION 0.1.0 LANGUAGES C CXX)


set(CMAKE_CXX_STANDARD 14)  
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

include(CTest)
enable_testing()

# Append the -g flag to the compiler options
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fsanitize=address")

file(GLOB SOURCES "src/*.cpp")
include_directories(/Applications/MATLAB_R2022b.app/extern/include)
include_directories(${CMAKE_SOURCE_DIR}/include/*)

add_executable(elf ${SOURCES})

target_link_libraries(elf PRIVATE  /Applications/MATLAB_R2022b.app/bin/maci64/libmat.dylib)
target_link_libraries(elf PRIVATE  /Applications/MATLAB_R2022b.app/bin/maci64/libmx.dylib)

set(CPACK_PROJECT_NAME ${PROJECT_NAME})
set(CPACK_PROJECT_VERSION ${PROJECT_VERSION})
include(CPack)
