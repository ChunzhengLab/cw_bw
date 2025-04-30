cmake_minimum_required(VERSION 3.10)
project(cw_bw)

# Find ROOT for histograms and I/O
find_package(ROOT REQUIRED COMPONENTS Core RIO Hist Tree)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_EXTENSIONS OFF)

# 查找 yaml-cpp 库
find_package(yaml-cpp REQUIRED)

# Include ROOT headers and definitions
include_directories(${ROOT_INCLUDE_DIRS} ${CMAKE_SOURCE_DIR}/include)
add_definitions(${ROOT_CXX_FLAGS})

# 源文件
file(GLOB LIB_SOURCES ${CMAKE_SOURCE_DIR}/src/*.cxx)

# 库 target
add_library(${PROJECT_NAME}_lib ${LIB_SOURCES})
target_link_libraries(${PROJECT_NAME}_lib PRIVATE yaml-cpp::yaml-cpp ${ROOT_LIBRARIES})


# 可执行程序：事件生成器
add_executable(bwgen apps/bwgen.cpp)
target_link_libraries(bwgen PRIVATE ${PROJECT_NAME}_lib ${ROOT_LIBRARIES})
