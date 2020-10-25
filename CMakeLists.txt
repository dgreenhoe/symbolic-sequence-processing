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
# Set Directory locations
#--------------------------------------
 set( DIR_GNU C:/CVSSource/Libraries/GNU/mingw-w64/i686-8.1.0-posix-dwarf-rt_v6-rev0/mingw32 )
 set( DIR_BCC C:/CVSSource/Libraries/bcc102 )
 set( DIR_TOOLCHAIN ${DIR_GNU} )
 set( DIR_TOOLCHAIN_INCLUDES  ${DIR_TOOLCHAIN}/include/windows/crtl )
 set( CMAKE_CXX_COMPILER ${DIR_TOOLCHAIN}/bin/g++.exe )

#--------------------------------------
# Set Flags
#--------------------------------------
 set( FLAGS_COMPILE
      #-v
      #-fpermissive
      -Wall
      -Wno-write-strings
      -Wno-return-type
      -Wno-format-overflow
    )

 set( FLAGS_LINK
      #-v
    )

#--------------------------------------
# Source Files
#--------------------------------------
 set( SOURCES_TARGET
      main.cpp
      test.cpp
    )

 set( SOURCES_CMD
      testlarc.cpp
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
# Builds
#--------------------------------------
 add_library(                Lib_Cmd   STATIC                                         )
 set_target_properties(      Lib_Cmd   PROPERTIES  PREFIX ""                          )
 target_compile_options(     Lib_Cmd   PRIVATE     ${FLAGS_COMPILE}                   )
 target_sources(             Lib_Cmd   PRIVATE     ${SOURCES_CMD}                     )
 target_include_directories( Lib_Cmd   PRIVATE     ${DIR_TOOLCHAIN_INCLUDES} ..       )

 add_library(                Lib_DSP   STATIC                                         )
 set_target_properties(      Lib_DSP   PROPERTIES PREFIX ""                           )
 target_compile_options(     Lib_DSP   PRIVATE    ${FLAGS_COMPILE}                    )
 target_sources(             Lib_DSP   PRIVATE    ${SOURCES_DSP}                      )
 target_include_directories( Lib_DSP   PRIVATE    ${DIR_TOOLCHAIN_INCLUDES} ..        )

 add_executable(             ${TARGET}                                                )
 set_target_properties(      ${TARGET} PROPERTIES OUTPUT_NAME ${TARGET} SUFFIX ".exe" )
 target_compile_options(     ${TARGET} PRIVATE    ${FLAGS_COMPILE}                    )
 target_link_options(        ${TARGET} PRIVATE    ${FLAGS_LINK}                       )
 target_link_libraries(      ${TARGET} PRIVATE    Lib_Cmd Lib_DSP                     )
 target_sources(             ${TARGET} PRIVATE    ${SOURCES_TARGET}                   )
 target_include_directories( ${TARGET} PRIVATE    ${DIR_TOOLCHAIN_INCLUDES} ..        )
