cmake_minimum_required(VERSION 3.12)
project(DILI)

set(CMAKE_CXX_STANDARD 17)

# Define the macro ‘DEBUG' in the debug mode
if(CMAKE_BUILD_TYPE STREQUAL Debug)
    ADD_DEFINITIONS(-DDEBUG)
endif()

if(MSVC)
    set(CMAKE_CXX_FLAGS "/O2 /arch:AVX2 /W1 /EHsc")
elseif (CMAKE_CXX_COMPILER_ID STREQUAL "Intel")
    set(CMAKE_CXX_FLAGS "-O3 -xHost")
else()
    # clang and gcc
    set(CMAKE_CXX_FLAGS "-lpthread -mpopcnt -O3 -march=native -Wextra -mlzcnt -mcmodel=large -fpermissive")
endif()

add_library(parse src/parse.cpp)

add_subdirectory(src/global)
add_subdirectory(src/dili)
add_subdirectory(src/butree)

include_directories(src/global)
include_directories(src/butree)
include_directories(src/dili)
#include_directories(src/competitors)


#add_executable(simple src/compare.cpp)
#target_compile_definitions(simple PRIVATE -DALLOW_FAN2_NODE)
#target_compile_definitions(simple_dili_func PRIVATE -DALLOW_FAN2_NODE)
#target_link_libraries(simple simple_dili_func)
#target_link_libraries(simple parse)
#target_link_libraries(simple interval_func)
#target_link_libraries(simple global_func)



#add_executable(compare src/compare.cpp)
#target_link_libraries(compare dili_func)
#target_link_libraries(compare parse)
#target_link_libraries(compare butree_func)
#target_link_libraries(compare global_func)

add_executable(example src/example.cpp)
target_link_libraries(example dili_func)
target_link_libraries(example parse)
target_link_libraries(example butree_func)
target_link_libraries(example global_func)

