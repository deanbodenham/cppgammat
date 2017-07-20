//first tests

#include<gtest/gtest.h>
#include "../gammatest/src/utils.h"

TEST(testTimesTwo, integerTests){
    EXPECT_EQ(0, timesTwo(0));
    EXPECT_EQ(2, timesTwo(1));
}


TEST(testTimesTwo, doubletests){
    EXPECT_EQ(2.4, timesTwo(1.2));
}

