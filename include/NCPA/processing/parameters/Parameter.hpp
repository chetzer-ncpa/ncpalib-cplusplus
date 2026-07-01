#pragma once

#include "NCPA/processing/declarations.hpp"

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
        class Parameter {
            public:
                Parameter() : _key { "" }, _description { "" } {}

                Parameter( const std::string& key ) :
                    _key { key }, _description { "" } {}

                Parameter( const std::string& key,
                           const std::string& description ) :
                    _key { key }, _description { description } {}

                Parameter( const Parameter& other ) :
                    _key { other._key }, _description { other._description } {}

                Parameter( Parameter&& other ) noexcept {
                    swap( *this, other );
                }

                virtual ~Parameter() {}

                friend void swap( Parameter& a, Parameter& b ) noexcept {
                    using std::swap;
                    swap( a._key, b._key );
                    swap( a._description, b._description );
                }

                virtual const std::string& key() const { return _key; }

                virtual const std::string& description() const {
                    return _description;
                }

                virtual size_t asSizeT() const {
                    int val = this->asInt();
                    if (val < 0) {
                        throw std::out_of_range(
                            "Value is negative, can't retrieve as size_t!" );
                    } else {
                        return static_cast<size_t>( val );
                    }
                }

                virtual std::vector<size_t> asSizeTVector() const {
                    std::vector<int> vi = this->asIntVector();
                    std::vector<size_t> vst( vi.size() );
                    for (size_t i = 0; i < vi.size(); ++i) {
                        if (vi[ i ] < 0) {
                            throw std::out_of_range(
                                "Value is negative, can't retrieve as "
                                "size_t!" );
                        } else {
                            vst[ i ] = static_cast<size_t>( vi[ i ] );
                        }
                    }
                    return vst;
                }

                virtual Parameter& clear()                              = 0;
                virtual std::unique_ptr<Parameter> clone() const        = 0;
                virtual bool is_scalar() const                          = 0;
                virtual parameter_form_t data_form() const              = 0;
                virtual std::string data_form_str() const               = 0;
                virtual parameter_type_t data_type() const              = 0;
                virtual std::string data_type_str() const               = 0;
                virtual std::string json_value_str() const              = 0;
                virtual bool asBool() const                             = 0;
                virtual int asInt() const                               = 0;
                virtual double asDouble() const                         = 0;
                virtual std::string asString() const                    = 0;
                virtual std::vector<bool> asBoolVector() const          = 0;
                virtual std::vector<int> asIntVector() const            = 0;
                virtual std::vector<double> asDoubleVector() const      = 0;
                virtual std::vector<std::string> asStringVector() const = 0;

                template<typename T>
                ScalarParameter<T>& as_scalar() {
                    return dynamic_cast<ScalarParameter<T>&>( *this );
                }

                template<typename T>
                const ScalarParameter<T>& as_scalar() const {
                    return dynamic_cast<const ScalarParameter<T>&>( *this );
                }

                template<typename T>
                VectorParameter<T>& as_vector() {
                    return dynamic_cast<VectorParameter<T>&>( *this );
                }

                template<typename T>
                const VectorParameter<T>& as_vector() const {
                    return dynamic_cast<const VectorParameter<T>&>( *this );
                }

                virtual std::string as_json_str( bool pretty = false,
                                                 int indent  = 1,
                                                 char space  = ' ' ) const {
                    if (pretty) {
                        return this->as_json_str_pretty( indent, space );
                    } else {
                        return this->as_json_str_plain();
                    }
                }

#if HAVE_NLOHMANN_JSON_HPP
                virtual nlohmann::json as_json() const {
                    nlohmann::json json;
                    json[ "comment" ] = this->description();
                    json[ "form" ]    = this->data_form_str();
                    json[ "type" ]    = this->data_type_str();
                    this->add_value_to_json( json );
                    return json;
                }

                virtual std::string as_json_str_plain() const {
                    return this->as_json().dump( -1, 0, true );
                }

                virtual std::string as_json_str_pretty( int indent = 1,
                                                        char space
                                                        = ' ' ) const {
                    return this->as_json().dump( indent, space, true );
                }

                virtual void add_value_to_json( nlohmann::json& json ) const
                    = 0;
#else
                virtual std::string as_json( bool pretty = false,
                                             int indent  = 1,
                                             char space  = '\t' ) const {
                    if (pretty) {
                        return this->as_json_str_pretty( indent, space );
                    } else {
                        return this->as_json_str_plain();
                    }
                }

                virtual std::string as_json_str_plain() const {
                    std::ostringstream json;
                    json << "\"" << this->key() << "\": {"
                         << "\"comment\": \"" << this->description() << "\", "
                         << "\"form\": \"" << this->data_form_str() << "\", "
                         << "\"type\": \"" << this->data_type_str() << "\", "
                         << "\"value\": " << this->json_value_str() << " }";
                    return json.str();
                }

                virtual std::string as_json_str_pretty( int indent = 1,
                                                        char space
                                                        = '\t' ) const {
                    std::ostringstream json;
                    // std::string indent = "";
                    // if (n_indent > 0) {
                    //     std::ostringstream indentstr;
                    //     for (int i = 0; i < n_indent; ++i) {
                    //         indentstr << space;
                    //     }
                    //     indent = indentstr.str();
                    // }
                    int n_indent = ( space == '\t' ? 1 : 4 );
                    std::ostringstream oss;
                    for (int i = 0; i < n_indent; ++i) {
                        oss << space;
                    }
                    std::string indentstr = oss.str();
                    oss.str( "" );
                    for (int i = 0; i < indent; ++i) {
                        oss << indentstr;
                    }
                    std::string baseindent = oss.str();

                    json << baseindent << "\"" << this->key() << "\": {\n"
                         << baseindent << indentstr << "\"comment\": \""
                         << this->description() << "\",\n"
                         << baseindent << indentstr << "\"form\": \""
                         << this->data_form_str() << "\",\n"
                         << baseindent << indentstr << "\"type\": \""
                         << this->data_type_str() << "\",\n"
                         << baseindent << indentstr
                         << "\"value\": " << this->json_value_str() << "\n"
                         << baseindent << "}";
                    return json.str();
                }
#endif
            private:
                std::string _key;
                std::string _description;
        };

        

        

        

        

        

        

        

        

        

        

        

        
    }  // namespace processing
}  // namespace NCPA

// void swap( NCPA::processing::Parameter& a,
//            NCPA::processing::Parameter& b ) noexcept {
//     using std::swap;
//     swap( a._key, b._key );
//     swap( a._description, b._description );
// }

// template<typename T>
// void swap( NCPA::processing::ScalarParameter<T>& a,
//            NCPA::processing::ScalarParameter<T>& b ) noexcept {
//     using std::swap;
//     swap( dynamic_cast<NCPA::processing::Parameter&>( a ),
//           dynamic_cast<NCPA::processing::Parameter&>( b ) );
// }

// template<typename T>
// void swap( NCPA::processing::VectorParameter<T>& a,
//            NCPA::processing::VectorParameter<T>& b ) noexcept {
//     using std::swap;
//     swap( dynamic_cast<NCPA::processing::Parameter&>( a ),
//           dynamic_cast<NCPA::processing::Parameter&>( b ) );
// }

// void swap( NCPA::processing::IntegerParameter& a,
//            NCPA::processing::IntegerParameter& b ) noexcept {
//     using std::swap;
//     swap( dynamic_cast<NCPA::processing::ScalarParameter<int>&>( a ),
//           dynamic_cast<NCPA::processing::ScalarParameter<int>&>( b ) );
// }

// void swap( NCPA::processing::DoubleParameter& a,
//            NCPA::processing::DoubleParameter& b ) noexcept {
//     using std::swap;
//     swap( dynamic_cast<NCPA::processing::ScalarParameter<double>&>( a ),
//           dynamic_cast<NCPA::processing::ScalarParameter<double>&>( b ) );
// }

// void swap( NCPA::processing::StringParameter& a,
//            NCPA::processing::StringParameter& b ) noexcept {
//     using std::swap;
//     swap( dynamic_cast<NCPA::processing::ScalarParameter<std::string>&>( a ),
//           dynamic_cast<NCPA::processing::ScalarParameter<std::string>&>( b ) );
// }

// void swap( NCPA::processing::BooleanParameter& a,
//            NCPA::processing::BooleanParameter& b ) noexcept {
//     using std::swap;
//     swap( dynamic_cast<NCPA::processing::ScalarParameter<bool>&>( a ),
//           dynamic_cast<NCPA::processing::ScalarParameter<bool>&>( b ) );
// }

// void swap( NCPA::processing::IntegerVectorParameter& a,
//            NCPA::processing::IntegerVectorParameter& b ) noexcept {
//     using std::swap;
//     swap( dynamic_cast<NCPA::processing::VectorParameter<int>&>( a ),
//           dynamic_cast<NCPA::processing::VectorParameter<int>&>( b ) );
// }

// void swap( NCPA::processing::DoubleVectorParameter& a,
//            NCPA::processing::DoubleVectorParameter& b ) noexcept {
//     using std::swap;
//     swap( dynamic_cast<NCPA::processing::VectorParameter<double>&>( a ),
//           dynamic_cast<NCPA::processing::VectorParameter<double>&>( b ) );
// }

// void swap( NCPA::processing::StringVectorParameter& a,
//            NCPA::processing::StringVectorParameter& b ) noexcept {
//     using std::swap;
//     swap( dynamic_cast<NCPA::processing::VectorParameter<std::string>&>( a ),
//           dynamic_cast<NCPA::processing::VectorParameter<std::string>&>( b ) );
// }

// void swap( NCPA::processing::BooleanVectorParameter& a,
//            NCPA::processing::BooleanVectorParameter& b ) noexcept {
//     using std::swap;
//     swap( dynamic_cast<NCPA::processing::VectorParameter<bool>&>( a ),
//           dynamic_cast<NCPA::processing::VectorParameter<bool>&>( b ) );
// }
