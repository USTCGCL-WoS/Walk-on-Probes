#pragma once

#include "core/macros.h"

#if BUILD_WOS_SCENE_MODULE
#define WOS_SCENE_API WOS_EXPORT
#define WOS_SCENE_EXTERN extern
#else
#define WOS_SCENE_API WOS_IMPORT
#if defined(_MSC_VER)
#define WOS_SCENE_EXTERN
#else
#define WOS_SCENE_EXTERN extern
#endif
#endif