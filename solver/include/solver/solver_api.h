#pragma once

#include "core/macros.h"

#if BUILD_WOS_SOLVER_MODULE
#define WOS_SOLVER_API WOS_EXPORT
#define WOS_SOLVER_EXTERN extern
#else
#define WOS_SOLVER_API WOS_IMPORT
#if defined(_MSC_VER)
#define WOS_SOLVER_EXTERN
#else
#define WOS_SOLVER_EXTERN extern
#endif
#endif
