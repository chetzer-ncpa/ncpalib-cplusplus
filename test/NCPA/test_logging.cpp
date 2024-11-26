#include "NCPA/logging.hpp"
#include "NCPA/gtest.hpp"

#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include <gtest/gtest-spi.h>

#include <sstream>
#include <fstream>
#include <string>

using namespace testing;
using namespace std;
using NCPA::logging::logger;
using NCPA::logging::log_level_t;

typedef double test_t;

#define _TEST_EQ_       EXPECT_DOUBLE_EQ
#define _TEST_ARRAY_EQ_ EXPECT_ARRAY_DOUBLE_EQ
#define _TEST_TITLE_    NCPALoggingLibraryTest

class _TEST_TITLE_ : public ::testing::Test {
    protected:
        void SetUp() override {  // define stuff here
        }  // void TearDown() override {}

        // declare stuff here
};

TEST_F( _TEST_TITLE_, LoggingToStdoutWorks ) {
    testing::internal::CaptureStdout();
    logger << log_level_t::ERROR << "The test passed if you can see this";
    std::string output = testing::internal::GetCapturedStdout();
    EXPECT_EQ( output, "The test passed if you can see this" );
}

TEST_F(_TEST_TITLE_, LoggingToStreamWorks) {
    std::ostringstream oss;
    logger.set_output( oss ).set_level( log_level_t::DEBUG );
    logger << log_level_t::DEBUG << "This is a test";
    EXPECT_EQ( oss.str(), "This is a test" );
}

TEST_F( _TEST_TITLE_, ChangingLogLevelsWorks) {
    std::ostringstream oss;
    logger.set_output( oss ).set_level( log_level_t::INFO );
    logger << log_level_t::DEBUG << "This is a test";
    EXPECT_EQ( oss.str(), "" );
    logger << log_level_t::INFO << "This is a test";
    EXPECT_EQ( oss.str(), "This is a test" );
}

TEST_F( _TEST_TITLE_, NewlineWorks ) {
    std::ostringstream oss;
    logger.set_output( oss ).set_level( log_level_t::DEBUG );
    logger << log_level_t::DEBUG << "This is a test" << endl;
    EXPECT_EQ( oss.str(), "This is a test\n" );
}

// TEST_F( _TEST_TITLE_, LoggingToFileWorks ) {
//     logger.set_output( "testlog.log" ).set_level( log_level_t::DEBUG );
//     cout << "Logger set up" << endl;
//     logger << log_level_t::ERROR << "This is a test";
//     cout << "Logger message sent" << endl;
//     logger.flush();
//     cout << "Logger flushed" << endl;
//     ifstream infile( "testlog.log" );
//     string logline;
//     getline( infile, logline );
//     EXPECT_EQ( logline, "This is a test");
// }

