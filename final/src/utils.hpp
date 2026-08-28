#pragma once

#include <objscip/objscip.h>
#include <objscip/objscipdefplugins.h>

#ifndef NDEBUG
#define CALL_CHECK(condition) \
  do {                        \
    assert((condition));      \
  } while (false);
#else
#define CALL_CHECK(condition) \
  do {                        \
    (condition);              \
  } while (false);
#endif
