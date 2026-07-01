#pragma once

#include "NCPA/processing/declarations.hpp"
#include "NCPA/processing/parameters/declarations.hpp"
#include "NCPA/processing/parameters/Parameter.hpp"

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
        template<typename T>
        class VectorParameter : public Parameter {
            public:
                VectorParameter( const std::string& key,
                                 const std::vector<T>& value,
                                 const std::string& description = "" ) :
                    Parameter( key, description ), _value { value } {}

                VectorParameter( const std::string& key,
                                 const std::initializer_list<T>& value,
                                 const std::string& description = "" ) :
                    Parameter( key, description ),
                    _value { std::vector<T> { value } } {}

                virtual ~VectorParameter() {}

                friend void swap(
                    NCPA::processing::VectorParameter<T>& a,
                    NCPA::processing::VectorParameter<T>& b ) noexcept {
                    using std::swap;
                    swap( dynamic_cast<NCPA::processing::Parameter&>( a ),
                          dynamic_cast<NCPA::processing::Parameter&>( b ) );
                }

                virtual Parameter& clear() override {
                    this->set( {} );
                    return *this;
                }

                virtual bool is_scalar() const override { return false; }

                virtual const std::vector<T>& value() const { return _value; }

                virtual VectorParameter& set( const std::vector<T>& val ) {
                    _value = val;
                    return *this;
                }

                virtual VectorParameter& set(
                    const std::initializer_list<T>& val ) {
                    _value = std::vector<T> { val };
                    return *this;
                }

#if HAVE_NLOHMANN_JSON_HPP
                virtual void add_value_to_json(
                    nlohmann::json& json ) const override {
                    json[ "value" ] = nlohmann::json( this->value() );
                }
#endif

                virtual bool asBool() const override {
                    return this->value().size() > 0
                             ? this->asBoolVector().at( 0 )
                             : false;
                }

                virtual int asInt() const override {
                    return this->value().size() > 0
                             ? this->asIntVector().at( 0 )
                             : 0;
                }

                virtual double asDouble() const override {
                    return this->value().size() > 0
                             ? this->asDoubleVector().at( 0 )
                             : 0.0;
                }

                virtual std::string asString() const override {
                    return this->value().size() > 0
                             ? this->asStringVector().at( 0 )
                             : "";
                }

                virtual parameter_form_t data_form() const override {
                    return parameter_form_t::ARRAY;
                }

                virtual std::string data_form_str() const override {
                    return "array";
                }

                virtual std::string json_value_str() const override {
                    std::ostringstream oss;
                    oss << "[ ";
                    for (auto it = _value.begin(); it != _value.end(); ++it) {
                        if (it != _value.begin()) {
                            oss << ", ";
                        }
                        oss << *it;
                    }
                    oss << " ]";
                    return oss.str();
                }

            private:
                std::vector<T> _value;
        };

        template<typename T>
        static std::unique_ptr<Parameter> make_parameter(
            const std::string& key, const std::vector<T>& val ) {
            return std::unique_ptr<Parameter>(
                new VectorParameter<T>( key, val ) );
        }

        template<typename T>
        static std::unique_ptr<Parameter> make_parameter(
            const std::string& key, const std::initializer_list<T>& val ) {
            return std::unique_ptr<Parameter>(
                new VectorParameter<T>( key, val ) );
        }
    }  // namespace processing
}  // namespace NCPA
