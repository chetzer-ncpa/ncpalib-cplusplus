#pragma once



namespace NCPA {
    template<typename state_t>
    class Stateful {
        public:
            Stateful() {}
            virtual ~Stateful() {}
            virtual const state_t& state() const = 0;
    };
}