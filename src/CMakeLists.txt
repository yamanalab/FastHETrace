cmake_minimum_required (VERSION 3.5.1)

project(demo CXX)
set(CMAKE_CXX_STANDARD 17)

set(CMAKE_MODULE_PATH ../cmake)

option(USE_FAST_AUTO "Set to ON to use table-based automorphism in PALISADE" ON)

#find_package(Palisade 1.9.2 EXACT REQUIRED PATHS $ENV{HOME}/lib/CMake/)
if(${WITH_INTEL_HEXL})
  find_package(Palisade 1.9.3 EXACT REQUIRED) # Intel HEXL-enabled
else()
  find_package(Palisade 1.9.2 EXACT REQUIRED)
endif()
if (${Palisade_FOUND})
  message("palisade found")
  message( ${PALISADE_INCLUDE})
  message( ${PALISADE_LIBDIR})
else (${Palisade_FOUND})
  message(FATAL_ERROR "palisade not found")
endif (${Palisade_FOUND})

message( ${PALISADE_INCLUDE})
message( ${PALISADE_INCLUDE}/core )
message( ${PALISADE_INCLUDE}/pke )

if (${WITH_INTEL_HEXL})
  find_package(HEXL 1.2.3 REQUIRED)
endif()

set( CMAKE_EXE_LINKER_FLAGS ${PALISADE_EXE_LINKER_FLAGS} )

include_directories(".")
include_directories("../include")
include_directories( ${PALISADE_INCLUDE} )
include_directories( ${PALISADE_INCLUDE}/core )
include_directories( ${PALISADE_INCLUDE}/pke )

link_directories( ${PALISADE_LIBDIR} )
link_libraries( ${PALISADE_LIBRARIES} )


set( CMAKE_CXX_FLAGS "${PALISADE_CXX_FLAGS} -march=native -O3 -Wall")

if (${USE_FAST_AUTO})
  message(" Fast automorph option ON ")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DUSE_FAST_AUTOMORPH")
endif()

if ("${CMAKE_BUILD_TYPE}" STREQUAL "Debug")
  add_compile_options("-fno-omit-frame-pointer")
  add_compile_options("-fsanitize=address")
  add_link_options("-fno-omit-frame-pointer")
  add_link_options("-fsanitize=address")
endif()

include_directories( ${OPENMP_INCLUDES} )
link_directories( ${OPENMP_LIBRARIES} )
message("CMAKE_CXX_FLAGS_DEBUG = ${CMAKE_CXX_FLAGS_DEBUG}")
message("CMAKE_CXX_FLAGS_RELEASE = ${CMAKE_CXX_FLAGS_RELEASE}")
message("CMAKE_CXX_FLAGS = ${CMAKE_CXX_FLAGS}")

add_executable(bench_unroll bench_unroll.cpp)
add_executable(bench_seq bench_seq.cpp)
