#include "NCPA/strings.hpp"
#include "NCPA/gtest.hpp"

#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include <gtest/gtest-spi.h>

#include <string>
#include <vector>

using namespace testing;
using namespace NCPA::strings;
using namespace std;

class StringsTest : public ::testing::Test {
    protected:
        void SetUp() override {  // define stuff here
            token_string = "  This is a multi-word string.";
            expected_tokens = { "This", "is", "a", "multi-word", "string."} ;
            mixed = "ThiS TExt Is miXeD CasE";
            lower = "this text is mixed case";
            upper = "THIS TEXT IS MIXED CASE";
        }  // void TearDown() override {}

        // declare stuff here
        std::string token_string, mixed, lower, upper;
        vector<std::string> expected_tokens;
        
};

TEST_F( StringsTest, DeblankRemovesFromFront ) {
    EXPECT_EQ( deblank( "     Test string"), "Test string" );
}

TEST_F( StringsTest, DeblankRemovesFromBack ) {
    EXPECT_EQ( deblank( "Test string        "), "Test string" );
}

TEST_F( StringsTest, DeblankRemovesFromFrontAndBack ) {
    EXPECT_EQ( deblank( "   Test string        "), "Test string" );
}

TEST_F( StringsTest, WhitespaceRedefineableForDeblank) {
    EXPECT_EQ( deblank( "   Test string        ", "."), "   Test string        " );
    EXPECT_EQ( deblank( "......Test string...", "."), "Test string" );
}

TEST_F( StringsTest, StringsSplitOnWhitespace) {
    vector<string> tokens = split( token_string );
    EXPECT_ARRAY_EQ( 5, tokens, expected_tokens );
}

TEST_F( StringsTest, StringsSplitOnCustomCharacters ) {
    vector<string> tokens = split( token_string, " \t\n\r." );
    expected_tokens[ 4 ] = "string";
    EXPECT_ARRAY_EQ( 5, tokens, expected_tokens );
}

TEST_F( StringsTest, TimeAsStringWorksAsExpected ) {
    EXPECT_EQ( time_as_string( 0.0 ), "1970-01-01T00:00:00.000Z");
}

TEST_F( StringsTest, ToLowerWorks ) {
    EXPECT_EQ( to_lower( mixed ), lower );
}

TEST_F( StringsTest, ToUpperWorks ) {
    EXPECT_EQ( to_upper( mixed ), upper );
}