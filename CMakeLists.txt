#==============================================================================
# CMake Script file
# Daniel J. Greenhoe
#==============================================================================

#--------------------------------------
# Build Environment Setup
#--------------------------------------
 cmake_minimum_required( VERSION 3.17 )
 set( CMAKE_VERBOSE_MAKEFILE   $<BOOL:true>                )
 set( CMAKE_C_COMPILER_WORKS   $<BOOL:true>                )
 set( CMAKE_CXX_COMPILER_WORKS $<BOOL:true>                )
 set( CMAKE_SYSTEM_NAME        CMAKE_HOST_SYSTEM_NAME      )
 set( CMAKE_SYSTEM_PROCESSOR   CMAKE_HOST_SYSTEM_PROCESSOR )
 set( CMAKE_CROSSCOMPILING     $<BOOL:false>               )
 set( CMAKE_C_STANDARD         99                          )
 set( CMAKE_CXX_STANDARD       11                          )

#--------------------------------------
# Project
#--------------------------------------
 project( opseq
          VERSION     0
          DESCRIPTION "Symbolic Sequence Processing"
         #HOMEPAGE    https://github.com/dgreenhoe/
          LANGUAGES   CXX
        )
 set( TARGET ${PROJECT_NAME} )

#--------------------------------------
# Directory locations
#--------------------------------------
 set( DIR_GNU C:/CVSSource/Libraries/GNU/mingw-w64/i686-8.1.0-posix-dwarf-rt_v6-rev0/mingw32 )
 set( DIR_BCC C:/CVSSource/Libraries/bcc102 )
 set( DIR_TOOLCHAIN ${DIR_GNU} )
 set( DIR_TOOLCHAIN_INCLUDES  ${DIR_TOOLCHAIN}/include/windows/crtl )
 set( DIR_EIGEN ../../libraries/eigen-3.3.8 )
 set( DIR_GTEST ../../libraries/googletest-1.10.0/googletest )
 set( DIR_GMOCK ../../libraries/googletest-1.10.0/googlemock )

#--------------------------------------
# Toolchain
#--------------------------------------
 set( CMAKE_CXX_COMPILER ${DIR_TOOLCHAIN}/bin/g++.exe )

#--------------------------------------
# Set Flags
#--------------------------------------
 set( FLAGS_COMPILE
      #-v
      -Wall
      -Werror
    )

 set( FLAGS_LINK
      #-v
      -Wall
      -Werror
      #-lgtest -lgtest_main
    )

#--------------------------------------
# Source Files
#ERROR using spin_spintoR1(c): c=0(0x30) is not in the valid domain {A,B,C,D,E,F}
#--------------------------------------
 set( SOURCES_TARGET
      main.cpp
      #tests/TestSuiteMisc.cpp
      #tests/TestSuiteR2.cpp
      #tests/TestSuiteR3.cpp
      #tests/TestSuiteR4.cpp
      #tests/TestSuiteR6.cpp
      #tests/TestSuiteLarc.cpp
      labs/LabOCS.cpp
      #labs/LabDNA.cpp
      #labs/LabLarc.cpp
    )

 set( SOURCES_CMD
      fileplot.cpp
    )

 set( SOURCES_DSP
      ballplot.cpp
      bsplines.cpp
      c1.cpp
      c4.cpp
      c6.cpp
      dft.cpp
      die.cpp
      dna.cpp
      dnan.cpp
      elliptic.cpp
      euclid.cpp
      fairdie.cpp
      fileplot.cpp
      lab2015larc.cpp
      lab2015ssp.cpp
      larc.cpp
      mca.cpp
      r1.cpp
      r1op.cpp
      r2.cpp
      r2op.cpp
      r3.cpp
      r3op.cpp
      r4.cpp
      r4op.cpp
      r6.cpp
      r6op.cpp
      realdie.cpp
      spinner.cpp
      symseq.cpp
    )

# set( SOURCES_TESTS
#      tests/TestSuiteMisc.cpp
#    )

 set( SOURCES_GTEST
      ${DIR_GTEST}/src/gtest-all.cc
      ${DIR_GTEST}/src/gtest-death-test.cc
      ${DIR_GTEST}/src/gtest-filepath.cc
      ${DIR_GTEST}/src/gtest-internal-inl.h
      ${DIR_GTEST}/src/gtest-matchers.cc
      ${DIR_GTEST}/src/gtest-port.cc
      ${DIR_GTEST}/src/gtest-printers.cc
      ${DIR_GTEST}/src/gtest-test-part.cc
      ${DIR_GTEST}/src/gtest-typed-test.cc
      ${DIR_GTEST}/src/gtest.cc
     #${DIR_GTEST}/src/gtest_main.cc
    )

#set( SOURCES_GMOCK
#     ${DIR_GMOCK}/src/gmock-all.cc
#     ${DIR_GMOCK}/src/gmock-cardinalities.cc
#     ${DIR_GMOCK}/src/gmock-internal-utils.cc
#     ${DIR_GMOCK}/src/gmock-matchers.cc
#     ${DIR_GMOCK}/src/gmock-spec-builders.cc
#     ${DIR_GMOCK}/src/gmock.cc
#     ${DIR_GMOCK}/src/gmock_main.cc
#   )

 set( HEADERS_TARGET
      main.h
    )

 set( HEADERS_CMD
      test.h
      testlarc.h
    )

 set( HEADERS_DSP
      ballplot.h
      bsplines.h
      c1.h
      c4.h
      c6.h
      dft.h
      die.h
      dna.h
      dnan.h
      elliptic.h
      euclid.h
      fairdie.h
      fileplot.h
      lab2015larc.h
      lab2015ssp.h
      larc.h
      main.h
      mca.h
      r1.h
      r1op.h
      r2.h
      r2op.h
      r3.h
      r3op.h
      r4.h
      r4op.h
      r6.h
      r6op.h
      realdie.h
      spinner.h
      symseq.h
    )

#--------------------------------------
# Build Libraries
#--------------------------------------
 add_library(                Lib_Cmd   STATIC                                               )
 set_target_properties(      Lib_Cmd   PROPERTIES  PREFIX ""                                )
 target_compile_options(     Lib_Cmd   PRIVATE     ${FLAGS_COMPILE}                         )
 target_sources(             Lib_Cmd   PRIVATE     ${SOURCES_CMD}                           )
 target_include_directories( Lib_Cmd   PRIVATE     ${DIR_TOOLCHAIN_INCLUDES} ${DIR_EIGEN} .. )
                                                                                            
 add_library(                Lib_DSP   STATIC                                               )
 set_target_properties(      Lib_DSP   PROPERTIES PREFIX ""                                 )
 target_compile_options(     Lib_DSP   PRIVATE    ${FLAGS_COMPILE}                          )
 target_sources(             Lib_DSP   PRIVATE    ${SOURCES_DSP}                            )
 target_include_directories( Lib_DSP   PRIVATE    ${DIR_TOOLCHAIN_INCLUDES} ${DIR_GTEST}/include/ ${DIR_GTEST} ${DIR_EIGEN} .. )

#add_library(                Lib_Tests STATIC                                               )
#set_target_properties(      Lib_Tests PROPERTIES PREFIX ""                                 )
#target_compile_options(     Lib_Tests PRIVATE    ${FLAGS_COMPILE}                          )
#target_sources(             Lib_Tests PRIVATE    ${SOURCES_TESTS}                          )
#target_include_directories( Lib_Tests PRIVATE    ${DIR_TOOLCHAIN_INCLUDES} ${PROJECT_SOURCE_DIR} ${DIR_GTEST}/include/ ${DIR_GTEST} .. )

 add_library(                Lib_Gtest STATIC                                               )
 set_target_properties(      Lib_Gtest PROPERTIES PREFIX ""                                 )
 target_compile_options(     Lib_Gtest PRIVATE    ${FLAGS_COMPILE}                          )
 target_sources(             Lib_Gtest PRIVATE    ${SOURCES_GTEST}                          )
 target_include_directories( Lib_Gtest PRIVATE    ${DIR_TOOLCHAIN_INCLUDES} ${DIR_GTEST}/include/ ${DIR_GTEST} .. )

#add_library(                Lib_Gmock STATIC                                               )
#set_target_properties(      Lib_Gmock PROPERTIES PREFIX ""                                 )
#target_compile_options(     Lib_Gmock PRIVATE    ${FLAGS_COMPILE}                          )
#target_sources(             Lib_Gmock PRIVATE    ${SOURCES_GMOCK}                          )
#target_include_directories( Lib_Gmock PRIVATE    ${DIR_TOOLCHAIN_INCLUDES} ${DIR_GTEST}/include/ ${DIR_GMOCK}/include/ ${DIR_GMOCK} .. )

#add_library(                Lib_GtestMainExt STATIC IMPORTED )
#set_target_properties(      Lib_GtestMainExt PROPERTIES C_EXTENSIONS OFF IMPORTED_LOCATION ${PROJECT_SOURCE_DIR}/libgtest_main.a )
#add_library(                Lib_GtestExt STATIC IMPORTED )
#set_target_properties(      Lib_GtestExt PROPERTIES C_EXTENSIONS OFF IMPORTED_LOCATION ${PROJECT_SOURCE_DIR}/libgtest.a )

#--------------------------------------
# Build Executable
#--------------------------------------
 add_executable(             ${TARGET}                                                      )
 target_link_libraries(      ${TARGET} PRIVATE    Lib_Cmd Lib_DSP Lib_Gtest                 )
#target_link_libraries(      ${TARGET} PRIVATE    Lib_Cmd Lib_DSP Lib_Gtest Lib_Tests       )
#target_link_libraries(      ${TARGET} PRIVATE    Lib_Cmd Lib_DSP Lib_GtestMainExt Lib_Tests)
 target_compile_options(     ${TARGET} PRIVATE    ${FLAGS_COMPILE}                          )
 target_link_options(        ${TARGET} PRIVATE    ${FLAGS_LINK}                             )
 target_sources(             ${TARGET} PRIVATE    ${SOURCES_TARGET}                         )
 target_include_directories( ${TARGET} PRIVATE    ${DIR_TOOLCHAIN_INCLUDES} ${PROJECT_SOURCE_DIR} ${DIR_GTEST}/include/ ..   )
 set_target_properties(      ${TARGET} PROPERTIES OUTPUT_NAME ${TARGET} SUFFIX ".exe"       )


#--------------------------------------
# Information
#--------------------------------------
message( STATUS "PROJECT_SOURCE_DIR = ${PROJECT_SOURCE_DIR}" )

