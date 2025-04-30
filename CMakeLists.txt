cmake_minimum_required(VERSION 3.10)
cmake_policy(SET CMP0074 NEW)

project(cmuts CXX)
set(TESTS "cmuts-generate-tests")

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_INSTALL_MESSAGE LAZY)

if (CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
    add_compile_options(-Wall)
endif()

find_package(PkgConfig REQUIRED)
find_package(MPI REQUIRED)
find_package(HDF5 REQUIRED)
pkg_check_modules(HTSLIB REQUIRED htslib)

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

# Add object library for shared sources
add_executable(${PROJECT_NAME}
    src/main.cpp
    src/cmuts.cpp
    src/hts.cpp
    src/mpi.cpp
    src/hdf5.cpp
    src/utils.cpp
)

target_include_directories(${PROJECT_NAME} PRIVATE
    ${HTSLIB_INCLUDE_DIRS}
    ${HDF5_INCLUDE_DIRS}
    ${MPI_C_INCLUDE_DIRS}
    deps/highfive/include
)

target_link_directories(${PROJECT_NAME} PRIVATE
    ${HTSLIB_LIBRARY_DIRS}
)

target_link_libraries(${PROJECT_NAME} PRIVATE
    ${HTSLIB_LIBRARIES}
    ${HDF5_LIBRARIES}
    MPI::MPI_CXX
    argparse
    xtensor
)

target_compile_options(${PROJECT_NAME} PRIVATE -O2)




add_executable(${TESTS}
    src/tests.cpp
    src/cmuts.cpp
    src/hts.cpp
    src/mpi.cpp
    src/hdf5.cpp
    src/utils.cpp
)

target_include_directories(${TESTS} PRIVATE
    ${HTSLIB_INCLUDE_DIRS}
    ${HDF5_INCLUDE_DIRS}
    ${MPI_C_INCLUDE_DIRS}
    deps/highfive/include
)

target_link_directories(${TESTS} PRIVATE
    ${HTSLIB_LIBRARY_DIRS}
)

target_link_libraries(${TESTS} PRIVATE
    ${HTSLIB_LIBRARIES}
    ${HDF5_LIBRARIES}
    MPI::MPI_CXX
    argparse
    xtensor
)

target_compile_options(${TESTS} PRIVATE -O2)



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
