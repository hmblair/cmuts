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
set(TESTS "_cmuts-generate-tests")

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
find_package(OpenMP REQUIRED)

include(FetchContent)
function(FetchContent_MakeAvailable_If_Not_Already_Present name)
    if (NOT TARGET ${name}::${name})
        FetchContent_MakeAvailable(${name})
    endif()
endfunction()

FetchContent_Declare(
    argparse
    GIT_REPOSITORY https://github.com/p-ranav/argparse.git
    GIT_TAG v3.2
)
FetchContent_Declare(
    xtl
    GIT_REPOSITORY https://github.com/xtensor-stack/xtl.git
    GIT_TAG 0.8.0
)
FetchContent_Declare(
    xtensor
    GIT_REPOSITORY https://github.com/xtensor-stack/xtensor.git
    GIT_TAG 0.27.0
)

FetchContent_MakeAvailable_If_Not_Already_Present(argparse)
FetchContent_MakeAvailable_If_Not_Already_Present(xtl)
FetchContent_MakeAvailable_If_Not_Already_Present(xtensor)


#
# Generate argument parsing headers from JSON
#

find_package(Python3 REQUIRED COMPONENTS Interpreter)

set(GENERATED_DIR ${CMAKE_CURRENT_SOURCE_DIR}/src/cpp/generated)
file(MAKE_DIRECTORY ${GENERATED_DIR})

add_custom_command(
    OUTPUT ${GENERATED_DIR}/cmuts_args.hpp
    COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/scripts/generate_args.py
            ${CMAKE_CURRENT_SOURCE_DIR}/config/cmuts_args.json
            ${GENERATED_DIR}/cmuts_args.hpp
            cmutsProgram
    DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/config/cmuts_args.json
            ${CMAKE_CURRENT_SOURCE_DIR}/scripts/generate_args.py
    COMMENT "Generating cmuts argument header"
)

add_custom_command(
    OUTPUT ${GENERATED_DIR}/tests_args.hpp
    COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/scripts/generate_args.py
            ${CMAKE_CURRENT_SOURCE_DIR}/config/tests_args.json
            ${GENERATED_DIR}/tests_args.hpp
            testsProgram
    DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/config/tests_args.json
            ${CMAKE_CURRENT_SOURCE_DIR}/scripts/generate_args.py
    COMMENT "Generating tests argument header"
)

add_custom_target(generate_args
    DEPENDS ${GENERATED_DIR}/cmuts_args.hpp ${GENERATED_DIR}/tests_args.hpp
)


#
# Shared static library
#


add_library(libcmuts STATIC
    src/cpp/io/fasta.cpp
    src/cpp/io/bam.cpp
    src/cpp/io/cram.cpp
    src/cpp/io/hdf5.cpp
    src/cpp/common.cpp
    src/cpp/core/cmuts.cpp
    src/cpp/infra/logging.cpp
    src/cpp/infra/mpi.cpp
    src/cpp/infra/utils.cpp
    src/cpp/infra/mutex.cpp
)

set_target_properties(libcmuts PROPERTIES PREFIX "")

target_include_directories(libcmuts PRIVATE
    ${HTSLIB_INCLUDE_DIRS}
    ${HDF5_INCLUDE_DIRS}
    ${PROJECT_SOURCE_DIR}/deps/highfive/include
    ${PROJECT_SOURCE_DIR}/htscodecs/htscodecs
    ${PROJECT_SOURCE_DIR}/src/cpp
    ${PROJECT_SOURCE_DIR}/src/cpp/compat
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
    OpenMP::OpenMP_CXX
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
    target_compile_options(libcmuts PRIVATE "-O3" "-fno-math-errno" "-ffast-math")
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
        ${PROJECT_SOURCE_DIR}/src/cpp
        ${PROJECT_SOURCE_DIR}/src/cpp/compat
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
        OpenMP::OpenMP_CXX
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
    src/cpp/app/main.cpp
)

add_executable(${TESTS}
    src/cpp/app/tests.cpp
)

setup_target(${CORE})
setup_target(${TESTS})

add_dependencies(${CORE} generate_args)
add_dependencies(${TESTS} generate_args)

set_target_properties(${CORE} ${TESTS} PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/bin
)


#
# Testing (CTest integration)
#


option(BUILD_TESTING "Build tests" ON)

if(BUILD_TESTING)
    enable_testing()

    # Python integration tests via pytest
    find_program(PYTEST_EXECUTABLE pytest)
    if(PYTEST_EXECUTABLE)
        # Quick tests (edge cases only, ~30 seconds)
        add_test(
            NAME python_integration_quick
            COMMAND ${PYTEST_EXECUTABLE}
                ${PROJECT_SOURCE_DIR}/tests/python/test_core.py::TestEdgeCases
                -v --tb=short
            WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
        )

        # Standard tests (edge cases + BAM fuzzing)
        add_test(
            NAME python_integration
            COMMAND ${PYTEST_EXECUTABLE}
                ${PROJECT_SOURCE_DIR}/tests/python/test_core.py
                -v --tb=short
                -m "not slow"
            WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
        )

        # CRAM-specific tests
        add_test(
            NAME python_integration_cram
            COMMAND ${PYTEST_EXECUTABLE}
                ${PROJECT_SOURCE_DIR}/tests/python/test_core.py::TestCramFormat
                -v --tb=short
            WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
        )

        # Full tests including stress tests (slow)
        add_test(
            NAME python_integration_full
            COMMAND ${PYTEST_EXECUTABLE}
                ${PROJECT_SOURCE_DIR}/tests/python
                -v --tb=short
            WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
        )

        # Set test properties
        set_tests_properties(python_integration_quick PROPERTIES
            LABELS "quick;python"
            TIMEOUT 120
        )
        set_tests_properties(python_integration PROPERTIES
            LABELS "python"
            TIMEOUT 600
        )
        set_tests_properties(python_integration_cram PROPERTIES
            LABELS "cram;python"
            TIMEOUT 300
        )
        set_tests_properties(python_integration_full PROPERTIES
            LABELS "slow;python"
            TIMEOUT 1800
        )
    else()
        message(WARNING "pytest not found, Python integration tests disabled")
    endif()
endif()
