#pragma once

#include "NCPA/configuration/declarations.hpp"

#include <memory>
#include <string>
#include <vector>

namespace NCPA {
    namespace config {

        class Argument {
            public:
                Argument() : Argument( "", "" ) {}

                Argument( const std::string& tag ) : Argument( tag, "" ) {}

                Argument( const std::string& tag,
                          const std::string helptext ) :
                    _tag { tag }, _help_text { helptext } {}

                virtual ~Argument() {}

                virtual std::string help_text() const { return _help_text; }

                virtual parse_status_t parse( size_t nargs, char **args ) {
                    std::vector<std::string> vargs( nargs );
                    for (size_t i = 0; i < nargs; ++i) {
                        vargs[ i ] = args[ i ];
                    }
                    return this->parse( vargs );
                }

                virtual bool& required() { return _required; }

                virtual const bool& required() const { return _required; }

                virtual Argument& required( bool tf ) {
                    _required = tf;
                    return *this;
                }

                virtual Argument& set_help_text( const std::string& text ) {
                    _help_text = text;
                    return *this;
                }

                virtual Argument& set_tag( const std::string& tag ) {
                    _tag = tag;
                    return *this;
                }

                virtual std::string tag() const { return _tag; }

                virtual void apply()               = 0;
                virtual bool expects_value() const = 0;
                virtual parse_status_t parse( std::vector<std::string>& args )
                    = 0;
                virtual bool was_set() const = 0;

            protected:
                std::string _tag;
                std::string _help_text;
                bool _required;
        };
    }  // namespace config
}  // namespace NCPA
