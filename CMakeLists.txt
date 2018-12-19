cmake_minimum_required (VERSION 3.5)

# project name
project (grid_precision)

# Sources
set (SRCS_EXE src/main.cpp src/import.cpp src/numerical_functions.cpp src/odes.cpp)

# Include directories
set (BOOST_ROOT "/usr/local/include/boost_1_69_0")
find_package(Boost REQUIRED)
INCLUDE_DIRECTORIES("${Boost_INCLUDE_DIRS}" include())

# Compiler options
set (DEBUG OFF)
if(DEBUG)
    add_compile_options(-Wall -O0 -g -std=c++17)
else()
    add_compile_options(-Wall -O3 -std=c++17)
endif(DEBUG)

# add the executable
set (EXE_NAME ${CMAKE_PROJECT_NAME})
add_executable (${EXE_NAME} ${SRCS_EXE})

# executable folder
set_target_properties(${EXE_NAME} 
    PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/bin )

target_link_libraries(${EXE_NAME} LINK_PUBLIC ${Boost_LIBRARIES})