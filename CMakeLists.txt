cmake_minimum_required(VERSION 3.10)
cmake_policy(SET CMP0074 NEW)

project(cmuts LANGUAGES C CXX)
set(TESTS "cmuts-generate-tests")

# Set C++ standard
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Display install messages lazily
set(CMAKE_INSTALL_MESSAGE LAZY)

# Compiler-specific options
if (CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
    add_compile_options(-Wall)
endif()

# Find required packages
find_package(PkgConfig REQUIRED)
if (DEFINED MPI_BUILD AND MPI_BUILD)
    find_package(MPI REQUIRED)
endif()
find_package(HDF5 REQUIRED)
pkg_check_modules(HTSLIB REQUIRED htslib)

# Fetch dependencies if not already available
include(FetchContent)

function(FetchContent_MakeAvailable_If_Not_Already_Present name)
    if (NOT TARGET ${name}::${name})
        FetchContent_MakeAvailable(${name})
    endif()
endfunction()

FetchContent_Declare(
    argparse
    GIT_REPOSITORY https://github.com/p-ranav/argparse.git
)
FetchContent_Declare(
    xtl
    GIT_REPOSITORY https://github.com/xtensor-stack/xtl.git
    GIT_TAG master
)
FetchContent_Declare(
    xtensor
    GIT_REPOSITORY https://github.com/xtensor-stack/xtensor.git
    GIT_TAG master
)

FetchContent_MakeAvailable_If_Not_Already_Present(argparse)
FetchContent_MakeAvailable_If_Not_Already_Present(xtl)
FetchContent_MakeAvailable_If_Not_Already_Present(xtensor)

# Create a static library for shared files
add_library(cmuts_shared
    src/cmuts.cpp
    src/hts.cpp
    src/mpi.cpp
    src/hdf5.cpp
    src/utils.cpp
)

if (DEFINED MPI_BUILD AND MPI_BUILD)
    target_include_directories(cmuts_shared PRIVATE ${MPI_C_INCLUDE_DIRS})
    target_link_libraries(cmuts_shared PRIVATE MPI::MPI_CXX)
endif()

target_include_directories(cmuts_shared PRIVATE
    ${HTSLIB_INCLUDE_DIRS}
    ${HDF5_INCLUDE_DIRS}
    ${PROJECT_SOURCE_DIR}/deps/highfive/include
)

target_link_directories(cmuts_shared PRIVATE
    ${HTSLIB_LIBRARY_DIRS}
)

target_link_libraries(cmuts_shared PRIVATE
    ${HTSLIB_LIBRARIES}
    ${HDF5_LIBRARIES}
    argparse
    xtensor
)

target_compile_options(cmuts_shared PRIVATE -O3)
if (DEFINED MPI_BUILD AND MPI_BUILD)
    target_compile_options(cmuts_shared PRIVATE "-DMPI_BUILD")
endif()


# Function to set up common properties for a target
function(setup_target target_name)

    if (DEFINED MPI_BUILD AND MPI_BUILD)
        target_include_directories(${target_name} PRIVATE ${MPI_C_INCLUDE_DIRS})
        target_link_libraries(${target_name} PRIVATE MPI::MPI_CXX)
    endif()

    target_include_directories(${target_name} PRIVATE
        ${HTSLIB_INCLUDE_DIRS}
        ${HDF5_INCLUDE_DIRS}
        ${PROJECT_SOURCE_DIR}/deps/highfive/include
    )

    target_link_directories(${target_name} PRIVATE
        ${HTSLIB_LIBRARY_DIRS}
    )

    target_link_libraries(${target_name} PRIVATE
        ${HTSLIB_LIBRARIES}
        ${HDF5_LIBRARIES}
        argparse
        xtensor
        cmuts_shared
    )

    target_compile_options(${target_name} PRIVATE -O3)
endfunction()

# Add executable for main project
add_executable(${PROJECT_NAME}
    src/main.cpp
)

setup_target(${PROJECT_NAME})

# Add executable for tests
add_executable(${TESTS}
    src/tests.cpp
)

setup_target(${TESTS})

if (DEFINED MPI_BUILD AND MPI_BUILD)
    target_compile_options(${PROJECT_NAME} PRIVATE "-DMPI_BUILD")
    target_compile_options(${TESTS} PRIVATE "-DMPI_BUILD")
endif()

# Enable sanitizers in debug mode
if (DEFINED DEBUG AND DEBUG)
    set(SANITIZER_COMPILE_FLAGS "-fsanitize=address" "-fsanitize=undefined" "-fno-omit-frame-pointer" "-g")
    set(SANITIZER_LINK_FLAGS "-fsanitize=address" "-fsanitize=undefined" "-fno-omit-frame-pointer" "-g")

    target_compile_options(${PROJECT_NAME} PRIVATE ${SANITIZER_COMPILE_FLAGS})
    target_link_options(${PROJECT_NAME} PRIVATE ${SANITIZER_LINK_FLAGS})

    target_compile_options(${TESTS} PRIVATE ${SANITIZER_COMPILE_FLAGS})
    target_link_options(${TESTS} PRIVATE ${SANITIZER_LINK_FLAGS})
endif()

# Installation target
install(TARGETS ${PROJECT_NAME} ${TESTS}
    RUNTIME DESTINATION bin
)
