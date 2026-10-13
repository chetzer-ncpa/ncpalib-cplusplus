#pragma once

#include <memory>
#include <utility>

namespace NCPA {
    namespace config {

        enum class help_text_marker_t {
            NO_MARKER,
            BREAK,
            NEWLINE,
            SET_HANGING_INDENT,
            SPACE
        };

        struct help_text_marker_options {
                std::string marker_end            = ">";
                std::string marker_contents_regex = R"(([^:]+)(?::(\d+))?)";
                std::string marker_start          = "<";
        };


        typedef std::pair<help_text_marker_t, int> marker_contents_t;

        class HelpTextSection;
        class HelpTextParagraphSection;
        class HelpTextOrganizerSection;
        class HelpTextArgumentSection;
        class HelpTextFormatter;
        struct HelpTextFormatterOptions;
        typedef std::unique_ptr<HelpTextSection> help_text_ptr_t;
        class HelpTextMarkers;

        typedef HelpTextOrganizerSection HelpText;

    }  // namespace config
}  // namespace NCPA
