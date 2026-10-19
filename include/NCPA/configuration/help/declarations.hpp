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

        enum class help_text_box_t { OFF = 0, ON, TITLE_ONLY };

        enum class help_text_align_t { LEFT = 0, CENTER, RIGHT };

        struct help_text_marker_options {
                std::string marker_end            = ">";
                std::string marker_contents_regex = R"(([^:]+)(?::(\d+))?)";
                std::string marker_start          = "<";
        };

        typedef std::pair<help_text_marker_t, int> marker_contents_t;

        class HelpTextSection;
        class HelpTextHeaderSection;
        class HelpTextParagraphSection;
        class HelpTextOrganizerSection;
        class HelpTextArgumentSection;
        class HelpTextFormatter;
        struct HelpTextFormatterOptions;
        typedef std::unique_ptr<HelpTextSection> help_text_ptr_t;
        class HelpTextMarkers;
        class HelpText;
        // typedef HelpTextOrganizerSection HelpText;

    }  // namespace config
}  // namespace NCPA
