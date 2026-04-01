cmake_minimum_required(VERSION 3.16.3)

# Set default build type if not specified
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
    message(STATUS "Build type not specified, defaulting to Release")
endif()
message(STATUS "Build type: ${CMAKE_BUILD_TYPE}")

project(
        BinnedDegreeCC
        VERSION 1.0
        LANGUAGES CXX
)

message(STATUS "Started CMake for ${PROJECT_NAME}...\n")

set(CMAKE_CXX_STANDARD 17)
set(PREPROCESS_DATASET_EXECUTABLE DataPreprocessing)
set(RUN_EXACT_EXECUTABLE RunExact)
set(RUN_ONE_PASS_ALGO_EXECUTABLE BOLIDE)


include(cmake/SetSourcesHeader.cmake)

add_compile_options(
        "-fexceptions"
        "$<$<CONFIG:Debug>:-g3;-ggdb>"
        "$<$<CONFIG:Release>:-O3>"          # Ensure O3 for release
        "$<$<CONFIG:RelWithDebInfo>:-O3>"   # Or use RelWithDebInfo build type
)

message(STATUS "Release flags: ${CMAKE_CXX_FLAGS_RELEASE}")
message(STATUS "Debug flags: ${CMAKE_CXX_FLAGS_DEBUG}")


add_executable(${RUN_EXACT_EXECUTABLE} ${header_files} ${sources})
add_executable(${RUN_ONE_PASS_ALGO_EXECUTABLE} ${header_files} ${sources})
add_executable(${PREPROCESS_DATASET_EXECUTABLE} ${header_files} ${sources})


target_include_directories(
        ${RUN_EXACT_EXECUTABLE}
        PUBLIC
        $<INSTALL_INTERFACE:include>
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
        PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}/src
)

target_include_directories(
        ${RUN_ONE_PASS_ALGO_EXECUTABLE}
        PUBLIC
        $<INSTALL_INTERFACE:include>
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
        PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}/src
)

target_include_directories(
        ${PREPROCESS_DATASET_EXECUTABLE}
        PUBLIC
        $<INSTALL_INTERFACE:include>
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
        PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}/src
)
