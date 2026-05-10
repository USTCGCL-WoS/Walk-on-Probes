#pragma once

#include "core/macros.h"

#if BUILD_WOS_CORE_MODULE
#define WOS_CORE_API WOS_EXPORT
#define WOS_CORE_EXTERN extern
#else
#define WOS_CORE_API WOS_IMPORT
#if defined(_MSC_VER)
#define WOS_CORE_EXTERN
#else
#define WOS_CORE_EXTERN extern
#endif
#endif