#include "NCPA/gtest.hpp"
#include "NCPA/processing.hpp"

#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include <chrono>
#include <string>
#include <vector>

#include <gtest/gtest-spi.h>

using namespace testing;
using namespace NCPA::processing;
using namespace std;

// class DoublingStep : public ProcessingStep<double,double> {

// };

class NCPAProcessingLibraryTest : public ::testing::Test {
    protected:
        void SetUp() override {  // define stuff here
            now = std::chrono::system_clock::now();

            set_double_wrapper = DataWrapper<double>( doubleval );
            input_ptr = input_ptr_t( new InputPacket( input_id_t::COMMAND, "test" ) );
            data_packet = DataPacket<double>( set_double_wrapper );

            data_packet_with_time = DataPacket<double>( doubleval, "has_time", now, std::chrono::minutes(1) );
        }  // void TearDown() override {}

        // declare stuff here
        double doubleval = -4.2;
        DataWrapper<double> double_wrapper;
        DataWrapper<double> set_double_wrapper;
        input_ptr_t input_ptr;
        DataPacket<double> data_packet, data_packet_with_time;
        time_point_t now;
};

TEST_F( NCPAProcessingLibraryTest, DataWrapperEvaluatesFalseIfUnset ) {
    EXPECT_FALSE( double_wrapper );
}

TEST_F( NCPAProcessingLibraryTest, DataWrapperEvaluatesTrueIfSet ) {
    EXPECT_TRUE( set_double_wrapper );
}

TEST_F( NCPAProcessingLibraryTest, DataWrapperReturnsExpectedContents ) {
    EXPECT_DOUBLE_EQ( set_double_wrapper.get(), doubleval );
}

TEST_F( NCPAProcessingLibraryTest, DataWrapperCanBeSet ) {
    ASSERT_FALSE( double_wrapper );
    double_wrapper.set( doubleval );
    ASSERT_TRUE( double_wrapper );
    EXPECT_DOUBLE_EQ( double_wrapper.get(), doubleval );
}

TEST_F( NCPAProcessingLibraryTest, InputPacketPointerReturnsCorrectID ) {
    EXPECT_EQ( input_ptr->ID(), input_id_t::COMMAND );
}

TEST_F( NCPAProcessingLibraryTest, InputPacketPointerReturnsCorrectTag ) {
    EXPECT_EQ( input_ptr->tag(), "test" );
}

TEST_F( NCPAProcessingLibraryTest, InputPacketTimeIsInPast ) {
    time_point_t now = std::chrono::system_clock::now();
    EXPECT_TRUE( input_ptr->packet_time() < now );
}

TEST_F( NCPAProcessingLibraryTest, DataPacketReturnsCorrectValue ) {
    EXPECT_DOUBLE_EQ( data_packet.get(), doubleval );
}

TEST_F( NCPAProcessingLibraryTest, DataPacketReturnsCorrectTime ) {
    EXPECT_EQ( data_packet_with_time.time(), now );
}

TEST_F( NCPAProcessingLibraryTest, DataPacketReturnsCorrectDuration ) {
    EXPECT_EQ( data_packet.duration(), duration_t::zero() );
    EXPECT_EQ( data_packet_with_time.duration(), std::chrono::seconds(60) );
}