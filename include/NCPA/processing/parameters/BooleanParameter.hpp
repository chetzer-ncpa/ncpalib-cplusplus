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
        class BooleanParameter : public ScalarParameter<bool> {
            public:
                BooleanParameter( const std::string& key, bool value,
                                  const std::string& description = "" ) :
                    ScalarParameter<bool>( key, value, description ) {}

                virtual ~BooleanParameter() {}

                friend void swap(
                    NCPA::processing::BooleanParameter& a,
                    NCPA::processing::BooleanParameter& b ) noexcept {
                    using std::swap;
                    swap(
                        dynamic_cast<NCPA::processing::ScalarParameter<bool>&>(
                            a ),
                        dynamic_cast<NCPA::processing::ScalarParameter<bool>&>(
                            b ) );
                }

                virtual Parameter& clear() override {
                    this->set( false );
                    return *this;
                }

                virtual std::unique_ptr<Parameter> clone() const override {
                    return std::unique_ptr<Parameter>(
                        new BooleanParameter( this->key(), this->value() ) );
                }

                virtual parameter_type_t data_type() const override {
                    return parameter_type_t::BOOLEAN;
                }

                virtual bool asBool() const override { return this->value(); }

                virtual int asInt() const override {
                    return ( this->value() ? 1 : 0 );
                }

                virtual double asDouble() const override {
                    return ( this->value() ? 1.0 : 0.0 );
                }

                virtual std::string asString() const override {
                    return ( this->value() ? "true" : "false" );
                }

                static std::unique_ptr<Parameter> build(
                    const std::string& key, bool value,
                    const std::string& description ) {
                    return std::unique_ptr<Parameter>(
                        new BooleanParameter( key, value, description ) );
                }

                static std::unique_ptr<Parameter> build(
                    const std::string& key, bool value ) {
                    return std::unique_ptr<Parameter>(
                        new BooleanParameter( key, value ) );
                }

                virtual std::string data_type_str() const override {
                    return "boolean";
                }
        };
    }  // namespace processing
}  // namespace NCPA
