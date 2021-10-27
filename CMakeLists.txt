cmake_minimum_required(VERSION 3.16)
project(bench-inv)

set(CMAKE_CXX_STANDARD 17)

# ------ stuff in src, for main binary
file(GLOB_RECURSE SRC LIST_DIRECTORIES true src/*.hpp src/*.cpp)
set(SRC ${SRC})

add_executable(
    ${PROJECT_NAME}
    ${SRC}
)

find_package(OpenMP)
if(OpenMP_CXX_FOUND)
    target_link_libraries(${PROJECT_NAME} PUBLIC OpenMP::OpenMP_CXX)
endif()

include_directories(src)
include_directories(lib/eigen)

add_subdirectory(lib/eigen)