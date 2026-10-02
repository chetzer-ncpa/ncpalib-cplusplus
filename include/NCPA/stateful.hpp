#pragma once

namespace NCPA {
    template<typename state_type>
    class Stateful {
        public:
            using state_t = state_type;

            Stateful() = default;

            Stateful( const state_type& st ) { _state = st; }

            Stateful( state_type&& st ) { _state = st; }

            virtual ~Stateful()                            = default;
            Stateful( const Stateful<state_type>& other )     = default;
            Stateful( Stateful<state_type>&& other ) noexcept = default;

            virtual const state_type& state() const { return _state; }

            virtual state_type& state() { return _state; }

        protected:
            state_type _state;
    };
}  // namespace NCPA
