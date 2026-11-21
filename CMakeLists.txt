# polysolver with CMake: the library (shared by default, -DBUILD_SHARED_LIBS=OFF for static), the
# CLI, the test suite (ctest), and an installed package for find_package(polysolver) and pkg-config.
# The Makefile is the development build; both compile src/*.c the same way.
cmake_minimum_required(VERSION 3.21)   # PROJECT_IS_TOP_LEVEL
file(READ "${CMAKE_CURRENT_SOURCE_DIR}/VERSION" PS_VERSION_TEXT)
string(STRIP "${PS_VERSION_TEXT}" PS_VERSION_TEXT)
project(polysolver VERSION ${PS_VERSION_TEXT} LANGUAGES C)

include(GNUInstallDirs)
include(CMakePackageConfigHelpers)
option(BUILD_SHARED_LIBS "Build libpolysolver as a shared library" ON)
option(POLYSOLVER_BUILD_TESTS "Build the test suite" ${PROJECT_IS_TOP_LEVEL})
option(POLYSOLVER_BUILD_CLI "Build the polysolver command-line tool" ${PROJECT_IS_TOP_LEVEL})

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")
find_package(GMP REQUIRED)

# the library sources: every src/*.c (token_list.c and token_stack.c are unused)
file(GLOB PS_SOURCES CONFIGURE_DEPENDS "${CMAKE_CURRENT_SOURCE_DIR}/src/*.c")
list(FILTER PS_SOURCES EXCLUDE REGEX "/token_(list|stack)\\.c$")

# compiled once with the internal headers and the ps_/psi_ rename map; the library, the CLI and
# the tests all use these objects (the CLI and tests call internals the shared library hides)
add_library(polysolver_objects OBJECT ${PS_SOURCES})
target_include_directories(polysolver_objects PUBLIC
    "${CMAKE_CURRENT_SOURCE_DIR}/src" "${CMAKE_CURRENT_SOURCE_DIR}/include")
target_compile_options(polysolver_objects PUBLIC -include "${CMAKE_CURRENT_SOURCE_DIR}/src/namespace.h")
target_link_libraries(polysolver_objects PUBLIC GMP::GMP)
set_target_properties(polysolver_objects PROPERTIES
    POSITION_INDEPENDENT_CODE ON C_VISIBILITY_PRESET hidden C_STANDARD 11)

add_library(polysolver $<TARGET_OBJECTS:polysolver_objects>)
add_library(polysolver::polysolver ALIAS polysolver)
target_include_directories(polysolver PUBLIC
    "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>" "$<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>")
target_link_libraries(polysolver PUBLIC GMP::GMP PRIVATE m)
set_target_properties(polysolver PROPERTIES
    VERSION ${PROJECT_VERSION} SOVERSION ${PROJECT_VERSION_MAJOR} PUBLIC_HEADER include/polysolver.h
    INSTALL_NAME_DIR "${CMAKE_INSTALL_FULL_LIBDIR}")       # macOS: an absolute id, so pkg-config users load it

if(POLYSOLVER_BUILD_CLI)
    add_executable(polysolver_cli cli/main.c)
    target_link_libraries(polysolver_cli PRIVATE polysolver_objects m)
    set_target_properties(polysolver_cli PROPERTIES OUTPUT_NAME polysolver)
    install(TARGETS polysolver_cli RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})
endif()

if(POLYSOLVER_BUILD_TESTS)
    enable_testing()
    add_executable(test_symbolic_solver tests/test_symbolic_solver.c tests/test_properties.c tests/test_subres_ref.c)
    target_link_libraries(test_symbolic_solver PRIVATE polysolver_objects m)
    target_compile_options(test_symbolic_solver PRIVATE -UNDEBUG)   # the suite is assert-based: keep them in Release
    add_test(NAME suite COMMAND test_symbolic_solver WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}")
    add_executable(clash tests/clash.c)                 # only the public header and library
    target_compile_options(clash PRIVATE -UNDEBUG)
    target_link_libraries(clash PRIVATE polysolver)
    add_test(NAME clash COMMAND clash)
endif()

# install: the library, the header, a CMake package and a pkg-config file
install(TARGETS polysolver EXPORT polysolverTargets
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR} ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR} PUBLIC_HEADER DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
set(PS_CMAKE_DIR ${CMAKE_INSTALL_LIBDIR}/cmake/polysolver)
install(EXPORT polysolverTargets NAMESPACE polysolver:: DESTINATION ${PS_CMAKE_DIR})
configure_package_config_file(cmake/polysolverConfig.cmake.in "${CMAKE_CURRENT_BINARY_DIR}/polysolverConfig.cmake"
    INSTALL_DESTINATION ${PS_CMAKE_DIR})
write_basic_package_version_file("${CMAKE_CURRENT_BINARY_DIR}/polysolverConfigVersion.cmake"
    COMPATIBILITY SameMinorVersion)
install(FILES "${CMAKE_CURRENT_BINARY_DIR}/polysolverConfig.cmake" "${CMAKE_CURRENT_BINARY_DIR}/polysolverConfigVersion.cmake"
    cmake/FindGMP.cmake DESTINATION ${PS_CMAKE_DIR})
set(PREFIX "${CMAKE_INSTALL_PREFIX}")
set(VERSION "${PROJECT_VERSION}")
configure_file(polysolver.pc.in "${CMAKE_CURRENT_BINARY_DIR}/polysolver.pc" @ONLY)
install(FILES "${CMAKE_CURRENT_BINARY_DIR}/polysolver.pc" DESTINATION ${CMAKE_INSTALL_LIBDIR}/pkgconfig)
