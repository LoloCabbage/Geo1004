cmake_minimum_required(VERSION 3.26)
project(hw01)


set(CMAKE_CXX_STANDARD 17)

find_package(CGAL REQUIRED)
find_library(GMP_LIBRARIES NAMES gmp PATHS /usr/lib/x86_64-linux-gnu)

add_executable(hw01 main.cpp)
target_link_libraries(hw01 CGAL::CGAL ${GMP_LIBRARIES})