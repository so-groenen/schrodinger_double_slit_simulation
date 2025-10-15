#ifndef WF_MANAGER_HPP
#define WF_MANAGER_HPP
#include <iostream>

#ifndef PI
#include <numbers>
#define PI std::numbers::pi_v<float>
#endif
#include <eigen3/Eigen/SparseLU> 

constexpr std::complex<float> imaginary_unit{0, 1}; 


class WaveFunctionBuilder
{
    float m_sigma = 0.5;
    float m_k     = 15*PI;
    float m_Lx{};
    float m_Ly{};
    float m_x0{};
    float m_y0{};
    Eigen::MatrixXcf m_psi{};
public:
    WaveFunctionBuilder(float Lx, float Ly, float x0 = 0, float y0 = 0);
    void set_initial_pos(float x0, float y0);
    void set_deviation(float sigma);
    auto build_wavefunction(size_t N_y, size_t N_x) const -> Eigen::Reshaped<Eigen::MatrixXcf, -1, 1>;  
private:
    void init_wavefunction(Eigen::MatrixXcf& Psi) const;
    void zero_out_boundary(Eigen::MatrixXcf& Psi) const;
};

#endif