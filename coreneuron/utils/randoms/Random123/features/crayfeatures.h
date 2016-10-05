/*
Copyright (c) 2014 EPFL-BBP, All rights reserved.

THIS SOFTWARE IS PROVIDED BY THE BLUE BRAIN PROJECT "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO,
THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE BLUE BRAIN PROJECT
BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR
BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE
OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN
IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

Note: Minimum/Initial version derived from openclfeatures.h to work
with cray compiler.
*/

#ifndef __crayfeatures_dot_hpp
#define __crayfeatures_dot_hpp

#ifndef R123_STATIC_INLINE
#define R123_STATIC_INLINE static __inline__
#endif

#ifndef R123_FORCE_INLINE
#define R123_FORCE_INLINE(decl) decl
#endif

#ifndef R123_CUDA_DEVICE
#define R123_CUDA_DEVICE
#endif

#ifndef R123_ASSERT
#include <assert.h>
#define R123_ASSERT(x) assert(x)
#endif

#ifndef R123_BUILTIN_EXPECT
#define R123_BUILTIN_EXPECT(expr, likely) expr
#endif

#ifndef R123_USE_WMMINTRIN_H
#define R123_USE_WMMINTRIN_H 0
#endif

#ifndef R123_USE_INTRIN_H
#define R123_USE_INTRIN_H 0
#endif

#ifndef R123_USE_MULHILO32_ASM
#define R123_USE_MULHILO32_ASM 0
#endif

#ifndef R123_USE_MULHILO64_ASM
#define R123_USE_MULHILO64_ASM 0
#endif

#ifndef R123_USE_MULHILO64_MSVC_INTRIN
#define R123_USE_MULHILO64_MSVC_INTRIN 0
#endif

#ifndef R123_USE_MULHILO64_CUDA_INTRIN
#define R123_USE_MULHILO64_CUDA_INTRIN 0
#endif

#ifndef R123_USE_MULHILO64_OPENCL_INTRIN
#define R123_USE_MULHILO64_OPENCL_INTRIN 0
#endif

#ifndef R123_USE_MULHILO64_MULHI_INTRIN
#define R123_USE_MULHILO64_MULHI_INTRIN (defined(__powerpc64__))
#endif

#ifndef R123_MULHILO64_MULHI_INTRIN
#define R123_MULHILO64_MULHI_INTRIN __mulhdu
#endif

#ifndef R123_USE_MULHILO32_MULHI_INTRIN
#define R123_USE_MULHILO32_MULHI_INTRIN 0
#endif

#ifndef R123_MULHILO32_MULHI_INTRIN
#define R123_MULHILO32_MULHI_INTRIN __mulhwu
#endif

#ifndef __STDC_CONSTANT_MACROS
#define __STDC_CONSTANT_MACROS
#endif
#include <stdint.h>
#ifndef UINT64_C
#error UINT64_C not defined.  You must define __STDC_CONSTANT_MACROS before you #include <stdint.h>
#endif

#endif
