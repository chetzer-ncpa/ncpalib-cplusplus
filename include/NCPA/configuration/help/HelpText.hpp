#pragma once

#include <iostream>
#include <string>
#include <cstdlib>
#include <stdexcept>

/*
Help text:
    Sections
        Section title
        Section contents    
        Subclasses:
            Command line
            Text
            Parameters & flags
            List
    TextWidth

Contains:
    Configurable
    AppName
    Option[]
        Parameter
        OptionToSetMapping
        HelpText
        Subclasses:
            Flag
            Value
                NumberValue
                TextValue

        
Parameter: 
    value, type, validation tests
*/