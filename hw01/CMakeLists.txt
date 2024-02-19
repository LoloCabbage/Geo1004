cmake_minimum_required(VERSION 3.26)
project(hw01)

set(CMAKE_CXX_STANDARD 17)
find_package(CGAL REQUIRED)
add_executable(hw01 main.cpp)
