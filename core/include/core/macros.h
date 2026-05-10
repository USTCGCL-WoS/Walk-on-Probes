#pragma once

#ifndef WOS_NAMESPACE_OPEN_SCOPE
#define WOS_NAMESPACE_OPEN_SCOPE namespace WOS {
#endif
#define WOS_NAMESPACE_CLOSE_SCOPE }

#if defined(_MSC_VER)
#define WOS_EXPORT __declspec(dllexport)
#define WOS_IMPORT __declspec(dllimport)
#else
#define WOS_EXPORT __attribute__((visibility("default")))
#define WOS_IMPORT
#endif