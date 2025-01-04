cmake_minimum_required(VERSION 3.12)
project(NullSpaceExample)

# Set the C++ standard
set(CMAKE_CXX_STANDARD 20)

# Find the Eigen library (replace /path/to/eigen with the actual path)
find_package(Eigen3 REQUIRED)

# Add the executable
add_executable(NullSpaceExample main.cpp)

# Link Eigen3 to the executable
target_link_libraries(NullSpaceExample Eigen3::Eigen)

target_compile_features(NullSpaceExample PRIVATE cxx_std_20)
target_compile_options(NullSpaceExample PRIVATE
	$<$<CXX_COMPILER_ID:GNU>:     $<$<NOT:$<CONFIG:Debug>>:-O3>>
	$<$<CXX_COMPILER_ID:Clang>:   $<$<NOT:$<CONFIG:Debug>>:-O3>>
	$<$<CXX_COMPILER_ID:MSVC>:    $<$<NOT:$<CONFIG:Debug>>:/Ot>>)