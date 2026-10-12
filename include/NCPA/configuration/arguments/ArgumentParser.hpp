#pragma once

#include "NCPA/configuration/declarations.hpp"
#include "NCPA/configuration/arguments/Argument.hpp"
#include "NCPA/cloneable.hpp"


namespace NCPA {
    namespace config {
        class ArgumentParser : public CloneBase<ArgumentParser> {
            public:
                ArgumentParser() = default;
                virtual ~ArgumentParser() = default;
                ArgumentParser( const ArgumentParser& other ) = default;
                ArgumentParser( ArgumentParser&& other ) noexcept = default;
                friend void swap( ArgumentParser& a, ArgumentParser& b ) noexcept {}
                virtual ArgumentParseResult parse( ArgumentSet& argset ) const = 0;
        };


    }
}