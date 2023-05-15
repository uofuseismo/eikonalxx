#include <iostream>
#include <CL/sycl.hpp>
#include <algorithm>
#include <cassert>
#if __has_include(<pstl/execution>)
   #include <pstl/execution>
   #include <pstl/algorithm>
   #define USE_PSTL
#endif
//#if __has_include(<oneapi/dpl/execution>)
//   #include <oneapi/dpl/execution>
//   #include <oneapi/dpl/algorithm>
//   #define ONE_API
//#else
//#endif
#include "eikonalxx/ray/gradientTracer2d.hpp"
#include "eikonalxx/abstractBaseClass/solver2d.hpp"
#include "private/grid.hpp"

