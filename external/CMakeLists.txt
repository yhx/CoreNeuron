# Copyright (c) 2017, Blue Brain Project
# All rights reserved.

# Redistribution and use in source and binary forms, with or without modification,
# are permitted provided that the following conditions are met:
# 1. Redistributions of source code must retain the above copyright notice,
#    this list of conditions and the following disclaimer.
# 2. Redistributions in binary form must reproduce the above copyright notice,
#    this list of conditions and the following disclaimer in the documentation
#    and/or other materials provided with the distribution.
# 3. Neither the name of the copyright holder nor the names of its contributors
#    may be used to endorse or promote products derived from this software
#    without specific prior written permission.

# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
# ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
# LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
# CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
# SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
# INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
# CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
# ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF
# THE POSSIBILITY OF SUCH DAMAGE.

include (ExternalProject)
include(FindPkgConfig)

#
# if NO BUNDLE is enable
# find external MOD2C using find_package system
# if not
#  use it as a submodule
if(CORENEURON_NO_BUNDLE)
    message(STATUS "Build configured with external MOD2C")
    find_package(MOD2C REQUIRED)

    # dummy alias for dependency resolution
    add_library(mod2c::mod2c_core INTERFACE IMPORTED GLOBAL)

else()
    message(STATUS "Build configured with embedded MOD2C")

    find_path(MOD2C_SUBMODULE NAMES CMakeLists.txt PATHS "${CMAKE_CURRENT_SOURCE_DIR}/mod2c")
    if (NOT MOD2C_SUBMODULE)
      message (FATAL_ERROR "missing mod2c submodule :run on top directory of your sources (git > 1.8.2 required): \
          git submodule update --init --recursive \
          ")
    endif()

    # take care of frontend / backend for cross-compilation
    set(OLD_CMAKE_C_COMPILER ${CMAKE_C_COMPILER})
    set(OLD_CMAKE_CXX_COMPILER ${CMAKE_CXX_COMPILER})
    if(FRONTEND_C_COMPILER)
        set(CMAKE_C_COMPILER ${FRONTEND_C_COMPILER})
    endif()
    if(FRONTEND_CXX_COMPILER)
        set(CMAKE_CXX_COMPILER ${FRONTEND_C_COMPILER})
    endif()

    # compile MOD2C
    add_subdirectory(mod2c)
    set(MOD2C_BINARY "${CMAKE_CURRENT_BINARY_DIR}/mod2c/src/mod2c_core/mod2c_core" CACHE STRING "")
    set(MOD2C_MODLUNIT "${CMAKE_CURRENT_SOURCE_DIR}/mod2c/share/nrnunits.lib" CACHE STRING "")

    # add mod2c as pre-dependency to make dependency tree for parallel build correct
    add_executable(mod2c::mod2c_core ALIAS mod2c_core)


    # restore     
    set(CMAKE_C_COMPILER "${OLD_CMAKE_C_COMPILER}")
    set(CMAKE_CXX_COMPILER "${OLD_CMAKE_CXX_COMPILER}")


endif()


list(APPEND CORENEURON_PRE_DEPENDENCIES mod2c::mod2c_core)
set(CORENEURON_PRE_DEPENDENCIES "${CORENEURON_PRE_DEPENDENCIES}" PARENT_SCOPE)
