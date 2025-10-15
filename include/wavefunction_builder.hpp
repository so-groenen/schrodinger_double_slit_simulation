

#ifndef WF_MANAGER_HPP
#define WF_MANAGER_HPP
#include<iostream>

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
    WaveFunctionBuilder(float Lx, float Ly, float x0 = 0, float y0 = 0)
        :   m_Lx{Lx}, m_Ly{Ly}, m_x0{x0}, m_y0{y0}
    {
    }
    void set_initial_pos(float x0, float y0)
    {
        m_x0 = x0;
        m_y0 = y0;
    }
    void set_deviation(float sigma)
    {
        m_sigma = sigma; 
    }
    auto build_wavefunction(size_t N_y, size_t N_x)
    {
        Eigen::MatrixXcf psi_temp = Eigen::MatrixXcf::Zero(N_y-2, N_x-2); 
        init_wavefunction(psi_temp);
        zero_out_boundary(psi_temp);
        return psi_temp.reshaped<Eigen::ColMajor>();
    }   
private:
    void init_wavefunction(Eigen::MatrixXcf& Psi)
    {
        Eigen::VectorXf x = Eigen::VectorXf::LinSpaced(Psi.cols(), 0, m_Lx);
        Eigen::VectorXf y = Eigen::VectorXf::LinSpaced(Psi.rows(), 0, m_Ly);
        float deltaX{};
        float deltaY{};
        float GaussArg{};
        std::complex<float> planeWave{};

        for (int iy = 0; iy < Psi.rows(); iy++)
        {
            for (int jx = 0; jx < Psi.cols(); jx++)
            {
                deltaX   = (x(jx) - m_x0);
                deltaY   = (y(iy) - m_y0);
                GaussArg = ( pow(deltaX, 2) + pow(deltaY, 2)) / pow(m_sigma, 2);

                planeWave  =  std::exp(imaginary_unit*m_k*deltaX);
                Psi(iy,jx) =  std::exp(-0.5f*GaussArg)*planeWave; //planeWave; // std::exp(-0.5*GaussArg) *
            }
        }
    }
    void zero_out_boundary(Eigen::MatrixXcf& Psi)
    {
        Psi.row(0           ).array() = 0;
        Psi.row(Psi.rows()-1).array() = 0;
        Psi.col(0           ).array() = 0;
        Psi.col(Psi.cols()-1).array() = 0;
    }
};

#endif