#pragma once

#include "NCPA/processing/declarations.hpp"
#include "NCPA/processing/parameters/declarations.hpp"
#include "NCPA/processing/parameters/Parameter.hpp"
#include "NCPA/processing/parameters/ScalarParameter.hpp"

#include <limits>
#include <memory>
#include <sstream>
#include <string>
#include <type_traits>
#include <vector>

#if HAVE_NLOHMANN_JSON_HPP
#  include "nlohmann/json.hpp"
#endif

namespace NCPA {
    namespace processing {
        class DoubleParameter : public ScalarParameter<double> {
            public:
                DoubleParameter( const std::string& key, double value,
                                 const std::string& description = "" ) :
                    ScalarParameter<double>( key, value, description ) {}

                virtual ~DoubleParameter() {}

                friend void swap(
                    NCPA::processing::DoubleParameter& a,
                    NCPA::processing::DoubleParameter& b ) noexcept {
                    using std::swap;
                    swap(
                        dynamic_cast<
                            NCPA::processing::ScalarParameter<double>&>( a ),
                        dynamic_cast<
                            NCPA::processing::ScalarParameter<double>&>( b ) );
                }

                virtual Parameter& clear() override {
                    this->set( 0.0 );
                    return *this;
                }

                virtual std::unique_ptr<Parameter> clone() const override {
                    return std::unique_ptr<Parameter>(
                        new DoubleParameter( this->key(), this->value() ) );
                }

                virtual parameter_type_t data_type() const override {
                    return parameter_type_t::FLOAT;
                }

                virtual bool asBool() const override {
                    return ( this->value() != 0.0 );
                }

                virtual int asInt() const override {
                    return static_cast<int>(
                        this->value()
                        + std::numeric_limits<double>::epsilon() );
                }

                virtual double asDouble() const override {
                    return this->value();
                }

                virtual std::string asString() const override {
                    std::ostringstream oss;
                    oss << this->value();
                    return oss.str();
                }

                static std::unique_ptr<Parameter> build(
                    const std::string& key, double value,
                    const std::string& description ) {
                    return std::unique_ptr<Parameter>(
                        new DoubleParameter( key, value, description ) );
                }

                static std::unique_ptr<Parameter> build(
                    const std::string& key, double value ) {
                    return std::unique_ptr<Parameter>(
                        new DoubleParameter( key, value ) );
                }

                virtual std::string data_type_str() const override {
                    return "float";
                }
        };
    }  // namespace processing
}  // namespace NCPA
