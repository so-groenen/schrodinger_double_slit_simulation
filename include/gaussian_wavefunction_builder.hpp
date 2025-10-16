#ifndef WF_MANAGER_HPP
#define WF_MANAGER_HPP
#include <iostream>
#include "interface_wavefunction_builder.hpp"

#ifndef PI
#include <numbers>
#define PI std::numbers::pi_v<float>
#endif


class GaussianWfBuilder : public IWaveFunctionBuilder
{
    float m_sigma {0.5};
    float m_k     {15*PI};
    float m_Lx{};
    float m_Ly{};
    float m_x0{};
    float m_y0{};
    Eigen::MatrixXcf m_psi{};
public:
    GaussianWfBuilder() = default;
    void set_system_size(float Lx, float Ly) override;
    void set_initial_pos(float x0, float y0) override;
    void set_deviation(float sigma) override;
    auto build_wavefunction(size_t N_y, size_t N_x) const -> Eigen::Reshaped<Eigen::MatrixXcf, -1, 1>  override;  
private:
    void init_wavefunction(Eigen::MatrixXcf& Psi) const;
    void zero_out_boundary(Eigen::MatrixXcf& Psi) const;
};

#endif