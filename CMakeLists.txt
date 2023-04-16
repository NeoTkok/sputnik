cmake_minimum_required(VERSION 3.25.2)

project(Spacecraft_maneuver)

set(CMAKE_CXX_STANDARD 17)


add_subdirectory(tests)


target_include_directories(TESTs PRIVATE src)

#set_target_properties(test_3diag PROPERTIES PUBLIC_HEADER src/ClassMat.hpp)