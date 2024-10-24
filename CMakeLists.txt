cmake_minimum_required(VERSION 3.10)
cmake_policy(SET CMP0074 NEW)
project(cmuts C)

# Set C standard
set(CMAKE_C_STANDARD 11)
set(H5_UTILS "h5utils")
set(HTSLIB_UTILS "htsutils")

# Check if using GNU or Clang compiler
if (CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
    add_compile_options(-Wall)
endif()

# Find MPI
find_package(MPI REQUIRED)

# Find HDF5
find_package(HDF5 REQUIRED)

# Find HTSlib
find_package(PkgConfig REQUIRED)
pkg_check_modules(HTSLIB REQUIRED htslib)

# Add executable
add_executable(${PROJECT_NAME} count/${PROJECT_NAME}.c count/${H5_UTILS}.c count/${HTSLIB_UTILS}.c)

# Include directories
target_include_directories(${PROJECT_NAME} PRIVATE
    ${HTSLIB_INCLUDE_DIRS}
    ${HDF5_INCLUDE_DIRS}
    ${MPI_C_INCLUDE_DIRS}
)

# Link directories
target_link_directories(${PROJECT_NAME} PRIVATE
    ${HTSLIB_LIBRARY_DIRS}
)

# Link libraries
target_link_libraries(${PROJECT_NAME} PRIVATE
    ${HTSLIB_LIBRARIES}
    ${HDF5_LIBRARIES}
    MPI::MPI_C
)

# Compiler options
target_compile_options(${PROJECT_NAME} PRIVATE -O3)

if(APPLE)
    # M1 Mac specific flags
    if(CMAKE_SYSTEM_PROCESSOR STREQUAL "arm64")
        target_include_directories(${PROJECT_NAME} PRIVATE /opt/homebrew/opt/libomp/include)
        target_link_directories(${PROJECT_NAME} PRIVATE /opt/homebrew/opt/libomp/lib)
    endif()
endif()

# Installation target
install(TARGETS ${PROJECT_NAME}
    RUNTIME DESTINATION bin
)
