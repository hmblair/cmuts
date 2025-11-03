cmake_minimum_required(VERSION 3.29)
cmake_policy(SET CMP0074 NEW)


set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_COMPILER_CLANG_SCAN_DEPS ON)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
set(CMAKE_INSTALL_MESSAGE LAZY)

if(EXISTS "/opt/homebrew/opt/llvm")
    set(LLVM_ROOT "/opt/homebrew/opt/llvm")
    link_directories("${LLVM_ROOT}/lib/c++")
endif()

project(cmuts LANGUAGES C CXX)
set(CORE "cmuts-core")
set(TESTS "cmuts-generate-tests")

# set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -march=native")
# set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -march=native")
set(SANITIZER_COMPILE_FLAGS "-fsanitize=address" "-fsanitize=undefined" "-fno-omit-frame-pointer" "-g" "-O1")
set(SANITIZER_LINK_FLAGS "-fsanitize=address" "-fsanitize=undefined" "-fno-omit-frame-pointer" "-g" "-O1")

if (CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
    add_compile_options(-Wall)
endif()


#
# Find/fetch packages
#


find_package(PkgConfig REQUIRED)
if (DEFINED MPI_BUILD AND MPI_BUILD)
    find_package(MPI REQUIRED)
    set(HDF5_PREFER_PARALLEL TRUE)
endif()
find_package(HDF5 REQUIRED)
find_package(ZLIB REQUIRED)
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


#
# Shared static library
#


add_library(libcmuts STATIC
    src/cpp/fasta.cpp
    src/cpp/common.cpp
    src/cpp/cmuts.cpp
    src/cpp/bam.cpp
    src/cpp/cram.cpp
    src/cpp/mpi.cpp
    src/cpp/hdf5.cpp
    src/cpp/utils.cpp
)

set_target_properties(libcmuts PROPERTIES PREFIX "")

target_include_directories(libcmuts PRIVATE
    ${HTSLIB_INCLUDE_DIRS}
    ${HDF5_INCLUDE_DIRS}
    ${PROJECT_SOURCE_DIR}/deps/highfive/include
    ${PROJECT_SOURCE_DIR}/htscodecs/htscodecs
)

target_link_directories(libcmuts PRIVATE
    ${HTSLIB_LIBRARY_DIRS}
    ${PROJECT_SOURCE_DIR}/htscodecs/lib
)

target_link_libraries(libcmuts PRIVATE
    ${HTSLIB_LIBRARIES}
    ${HDF5_LIBRARIES}
    ${ZLIB_LIBRARIES}
    argparse
    xtensor
    htscodecs
)

target_compile_definitions(libcmuts PRIVATE VERSION="${VERSION}")

if (DEFINED MPI_BUILD AND MPI_BUILD)
    target_include_directories(libcmuts PRIVATE ${MPI_C_INCLUDE_DIRS})
    target_link_libraries(libcmuts PRIVATE MPI::MPI_CXX)
    target_compile_options(libcmuts PRIVATE "-DMPI_BUILD")
endif()

if (DEFINED DEBUG AND DEBUG)
    target_compile_options(libcmuts PRIVATE ${SANITIZER_COMPILE_FLAGS})
    target_link_options(libcmuts PRIVATE ${SANITIZER_LINK_FLAGS})
    target_link_options(libcmuts PRIVATE "-Wl,-v")
    target_compile_options(libcmuts PRIVATE "-DDEBUG")
else()
    target_compile_options(libcmuts PRIVATE -O3 -fno-math-errno)
endif()


#
# Executables and installation
#


function(setup_target _target)

    target_include_directories(${_target} PRIVATE
        ${HTSLIB_INCLUDE_DIRS}
        ${HDF5_INCLUDE_DIRS}
        ${PROJECT_SOURCE_DIR}/deps/highfive/include
        ${PROJECT_SOURCE_DIR}/htscodecs/htscodecs
    )

    target_link_directories(${_target} PRIVATE
        ${HTSLIB_LIBRARY_DIRS}
        ${PROJECT_SOURCE_DIR}/htscodecs/lib
    )

    target_link_libraries(${_target} PRIVATE
        ${HTSLIB_LIBRARIES}
        ${HDF5_LIBRARIES}
        ${ZLIB_LIBRARIES}
        argparse
        xtensor
        htscodecs
        libcmuts
        # "/opt/homebrew/opt/llvm/lib/c++/libc++.dylib"
        # "/opt/homebrew/opt/llvm/lib/c++/libc++abi.dylib"
    )

    set_target_properties(${_target} PROPERTIES
        INSTALL_RPATH "${PROJECT_SOURCE_DIR}/htscodecs/lib"
        BUILD_WITH_INSTALL_RPATH TRUE
    )

    target_compile_definitions(${_target} PRIVATE VERSION="${VERSION}")

    if (DEFINED MPI_BUILD AND MPI_BUILD)
        target_include_directories(${_target} PRIVATE ${MPI_C_INCLUDE_DIRS})
        target_link_libraries(${_target} PRIVATE MPI::MPI_CXX)
        target_compile_options(${_target} PRIVATE "-DMPI_BUILD")
    endif()

    if (DEFINED DEBUG AND DEBUG)
        target_compile_options(${_target} PRIVATE ${SANITIZER_COMPILE_FLAGS})
        target_link_options(${_target} PRIVATE ${SANITIZER_LINK_FLAGS})
        target_link_options(${_target} PRIVATE "-Wl,-v")
        target_compile_options(${_target} PRIVATE "-DDEBUG")
    else()
        target_compile_options(${_target} PRIVATE -O3 -fno-math-errno)
    endif()

endfunction()

add_executable(${CORE}
    src/cpp/main.cpp
)

add_executable(${TESTS}
    src/cpp/tests.cpp
)

setup_target(${CORE})
setup_target(${TESTS})

install(TARGETS ${CORE} ${TESTS}
    RUNTIME DESTINATION ${PROJECT_SOURCE_DIR}/bin
)
