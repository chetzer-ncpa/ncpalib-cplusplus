#pragma once

#include "NCPA/atmosphere/abstract/atmospheric_property.hpp"

namespace NCPA {
    namespace atmos {
        namespace abstract {
            class atmospheric_property_scalar : public atmospheric_property {
                public:
                    virtual ~atmospheric_property_scalar() {}

                    virtual size_t dimensions() const override { return 0; }

                    virtual double f() = 0;

                    virtual double get( const std::vector<double>& coords ) override {
                        return this->f();
                    }

                    virtual double get_first_derivative(
                        const std::vector<double>& coords,
                        std::vector<size_t>& rel ) override {
                        return 0.0;
                    }

                    virtual double get_second_derivative(
                        const std::vector<double>& coords,
                        std::vector<size_t>& rel ) override {
                        return 0.0;
                    }

                    virtual double get_first_derivative(
                        const std::vector<double>& coords ) override {
                        return 0.0;
                    }

                    virtual double get_second_derivative(
                        const std::vector<double>& coords ) override {
                        return 0.0;
                    }
            };
        }  // namespace abstract
    }  // namespace atmos
}  // namespace NCPA
