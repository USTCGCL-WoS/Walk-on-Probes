#include "core/scalar.hpp"

#include <gtest/gtest.h>

using namespace WOS;

TEST(ScalarTest, NaNBehavior)
{
    Scalar<3> nanVal = Scalar<3>::NaN();
    EXPECT_TRUE(nanVal.isNaN());
}