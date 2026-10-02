cmake_minimum_required(VERSION 3.4...3.18)
project(PyMatching2)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

include_directories(src)
set(CMAKE_CXX_STANDARD 20 CACHE STRING "C++ version selection")
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(SIMD_WIDTH 128)
if (NOT(MSVC))
    if (CMAKE_SYSTEM_PROCESSOR MATCHES x86_64)
         set(ARCH_OPT "-O3" "-mno-avx2")
    else ()
         set(ARCH_OPT "-O3")
    endif ()
else ()
    set(ARCH_OPT "-O2")
    # https://stackoverflow.com/a/8591946 (RTC1 incompatible with O2 flag)
    STRING (REGEX REPLACE "/RTC(su|[1su])" "" CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}")
endif ()

include(FetchContent)
FetchContent_Declare(
        googletest
        GIT_REPOSITORY https://github.com/google/googletest.git
        GIT_TAG release-1.12.1
)

FetchContent_MakeAvailable(googletest)

FetchContent_Declare(stim
        GIT_REPOSITORY https://github.com/quantumlib/Stim.git 
        GIT_TAG 1320ad7eac7de34d2e9c70daa44fbc6d84174450)
FetchContent_MakeAvailable(stim)
if (NOT (MSVC))
    target_compile_options(libstim PRIVATE -fno-strict-aliasing -fPIC ${ARCH_OPT})
else ()
    target_compile_options(libstim PRIVATE -fPIC ${ARCH_OPT})
endif ()

set(SOURCE_FILES_NO_MAIN
        src/pymatching/sparse_blossom/driver/namespaced_main.cc
        src/pymatching/sparse_blossom/driver/mwpm_decoding.cc
        src/pymatching/sparse_blossom/flooder/graph.cc
        src/pymatching/sparse_blossom/flooder/detector_node.cc
        src/pymatching/sparse_blossom/flooder_matcher_interop/compressed_edge.cc
        src/pymatching/sparse_blossom/flooder/graph_fill_region.cc
        src/pymatching/sparse_blossom/flooder/match.cc
        src/pymatching/sparse_blossom/flooder/graph_flooder.cc
        src/pymatching/sparse_blossom/matcher/alternating_tree.cc
        src/pymatching/sparse_blossom/matcher/mwpm.cc
        src/pymatching/sparse_blossom/flooder_matcher_interop/region_edge.cc
        src/pymatching/sparse_blossom/flooder_matcher_interop/mwpm_event.cc
        src/pymatching/sparse_blossom/tracker/flood_check_event.cc
        src/pymatching/sparse_blossom/diagram/animation_main.cc
        src/pymatching/sparse_blossom/diagram/mwpm_diagram.cc
        src/pymatching/sparse_blossom/search/search_graph.cc
        src/pymatching/sparse_blossom/search/search_detector_node.cc
        src/pymatching/sparse_blossom/search/search_flooder.cc
        src/pymatching/sparse_blossom/driver/user_graph.cc
        src/pymatching/rand/rand_gen.cc
        )

# fastgap: certified complementary-gap lower bounds from the decoder's final dual state.
set(FASTGAP_SOURCE_FILES
        src/pymatching/fastgap/canonical.cc
        src/pymatching/fastgap/interior_graph.cc
        src/pymatching/fastgap/dual_state.cc
        src/pymatching/fastgap/gap_index.cc
        src/pymatching/fastgap/gap_decoder.cc
        src/pymatching/fastgap/exact.cc
        )

set(FASTGAP_TEST_FILES
        src/pymatching/fastgap/canonical.test.cc
        src/pymatching/fastgap/interior_graph.test.cc
        src/pymatching/fastgap/dual_state.test.cc
        src/pymatching/fastgap/gap_index.test.cc
        src/pymatching/fastgap/gap_decoder.test.cc
        src/pymatching/fastgap/exact.test.cc
        )

set(TEST_FILES
        src/pymatching/sparse_blossom/driver/namespaced_main.test.cc
        src/pymatching/sparse_blossom/driver/mwpm_decoding.test.cc
        src/pymatching/sparse_blossom/flooder_matcher_interop/varying.test.cc
        src/pymatching/sparse_blossom/flooder/graph.test.cc
        src/pymatching/sparse_blossom/flooder/detector_node.test.cc
        src/pymatching/sparse_blossom/flooder_matcher_interop/compressed_edge.test.cc
        src/pymatching/sparse_blossom/flooder/graph_fill_region.test.cc
        src/pymatching/sparse_blossom/flooder/match.test.cc
        src/pymatching/sparse_blossom/flooder/graph_flooder.test.cc
        src/pymatching/sparse_blossom/matcher/alternating_tree.test.cc
        src/pymatching/sparse_blossom/matcher/mwpm.test.cc
        src/pymatching/sparse_blossom/tracker/flood_check_event.test.cc
        src/pymatching/sparse_blossom/tracker/radix_heap_queue.test.cc
        src/pymatching/sparse_blossom/flooder_matcher_interop/mwpm_event.test.cc
        src/pymatching/sparse_blossom/tracker/queued_event_tracker.test.cc
        src/pymatching/sparse_blossom/tracker/cyclic.test.cc
        src/pymatching/sparse_blossom/diagram/mwpm_diagram.test.cc
        src/pymatching/sparse_blossom/search/search_graph.test.cc
        src/pymatching/sparse_blossom/search/search_flooder.test.cc
        src/pymatching/sparse_blossom/driver/user_graph.test.cc
        )

set(PERF_FILES
        src/pymatching/perf/main.perf.cc
        src/pymatching/perf/util.perf.cc
        src/pymatching/sparse_blossom/driver/mwpm_decoding.perf.cc
        src/pymatching/sparse_blossom/driver/user_graph.perf.cc
        src/pymatching/sparse_blossom/flooder_matcher_interop/varying.perf.cc
        src/pymatching/sparse_blossom/tracker/radix_heap_queue.perf.cc
        )

set(PYTHON_API_FILES
        src/pymatching/sparse_blossom/driver/user_graph.pybind.cc
        src/pymatching/rand/rand_gen.pybind.cc
        src/pymatching/pymatching.pybind.cc
        )

add_executable(pymatching src/main.cc ${SOURCE_FILES_NO_MAIN})

target_compile_options(pymatching PRIVATE ${ARCH_OPT})
if (NOT (MSVC))
    target_link_options(pymatching PRIVATE -pthread ${ARCH_OPT})
endif ()
target_link_libraries(pymatching libstim)
install(TARGETS pymatching RUNTIME DESTINATION bin)

enable_testing()

add_executable(pymatching_tests ${SOURCE_FILES_NO_MAIN} ${TEST_FILES} ${FASTGAP_SOURCE_FILES} ${FASTGAP_TEST_FILES})
# fastgap's debug invariants (CLAUDE.md §5.1) are always on in tests.
target_compile_definitions(pymatching_tests PRIVATE FASTGAP_CHECK)
target_compile_options(pymatching_tests PRIVATE -fsanitize=address -g -fsanitize=undefined -fno-omit-frame-pointer -coverage)
if (NOT (MSVC))
    target_link_options(pymatching_tests PRIVATE -pthread -fsanitize=address -fsanitize=undefined -coverage)
else ()
    target_link_options(pymatching_tests PRIVATE -fsanitize=address -fsanitize=undefined -coverage)
endif ()
target_link_libraries(pymatching_tests GTest::gtest_main GTest::gmock_main libstim)

add_executable(pymatching_perf ${SOURCE_FILES_NO_MAIN} ${PERF_FILES})
target_compile_options(pymatching_perf PRIVATE ${ARCH_OPT})
if (NOT (MSVC))
    target_link_options(pymatching_perf PRIVATE -pthread -O3)
endif ()
target_link_libraries(pymatching_perf libstim)

add_library(libpymatching ${SOURCE_FILES_NO_MAIN})
set_target_properties(libpymatching PROPERTIES PREFIX "")
target_include_directories(libpymatching PUBLIC src)
target_compile_options(libpymatching PRIVATE ${ARCH_OPT} -fPIC)
if(NOT(MSVC))
    target_link_options(libpymatching PRIVATE -pthread -O3)
endif()
target_link_libraries(libpymatching libstim)
install(TARGETS libpymatching LIBRARY DESTINATION)
install(DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/src/" DESTINATION "include" FILES_MATCHING PATTERN "*.h" PATTERN "*.inl")

include(GoogleTest)
gtest_discover_tests(pymatching_tests)

IF(EXISTS ${CMAKE_CURRENT_LIST_DIR}/pybind11/CMakeLists.txt)
    add_subdirectory(pybind11)
    pybind11_add_module(_cpp_pymatching ${PYTHON_API_FILES} ${SOURCE_FILES_NO_MAIN})
    target_link_libraries(_cpp_pymatching PRIVATE libstim)
    target_compile_options(_cpp_pymatching PRIVATE ${ARCH_OPT})

    pybind11_add_module(_cpp_fastgap src/pymatching/fastgap/fastgap.pybind.cc ${SOURCE_FILES_NO_MAIN} ${FASTGAP_SOURCE_FILES})
    target_link_libraries(_cpp_fastgap PRIVATE libstim)
    target_compile_options(_cpp_fastgap PRIVATE ${ARCH_OPT})
    # setup.py builds both extensions from one CMake tree; each lands in its own package directory.
    if (DEFINED PYMATCHING_PY_OUTPUT_DIR)
        set_target_properties(_cpp_pymatching PROPERTIES LIBRARY_OUTPUT_DIRECTORY ${PYMATCHING_PY_OUTPUT_DIR})
    endif ()
    if (DEFINED FASTGAP_PY_OUTPUT_DIR)
        set_target_properties(_cpp_fastgap PROPERTIES LIBRARY_OUTPUT_DIRECTORY ${FASTGAP_PY_OUTPUT_DIR})
    endif ()
else()
    message("WARNING: Skipped the pybind11 module _cpp_pymatching because the `pybind11` git submodule isn't present. To fix, run `git submodule update --init --recursive`")
endif()

# EXAMPLE_VERSION_INFO is defined by setup.py and passed into the C++ code as a
# define (VERSION_INFO) here.
#target_compile_definitions(_cpp_pymatching
#        PRIVATE VERSION_INFO=${PYMATCHING_VERSION_INFO})
