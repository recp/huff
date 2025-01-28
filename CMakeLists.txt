cmake_minimum_required(VERSION 3.13)
project(huff
  VERSION 0.1.0
  HOMEPAGE_URL https://github.com/recp/huff
  DESCRIPTION "fast & flexible huffman library"
  LANGUAGES C
)

add_library(huff INTERFACE)
target_include_directories(huff
    INTERFACE
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
        $<INSTALL_INTERFACE:include>
)

include(GNUInstallDirs)

install(TARGETS huff
        EXPORT huffTargets
        INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})

install(DIRECTORY include/
        DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})

install(EXPORT huffTargets
        FILE huffTargets.cmake
        NAMESPACE huff::
        DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/huff)

set(CPACK_PROJECT_NAME ${PROJECT_NAME})
set(CPACK_PROJECT_VERSION ${PROJECT_VERSION})

if(NOT CPack_CMake_INCLUDED)
  include(CPack)
endif()
