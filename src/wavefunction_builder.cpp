#include "wavefunction_builder.hpp"

 
WaveFunctionBuilder::WaveFunctionBuilder(float Lx, float Ly, float x0, float y0)
    :   m_Lx{Lx}, m_Ly{Ly}, m_x0{x0}, m_y0{y0}
{
}
void WaveFunctionBuilder::set_initial_pos(float x0, float y0)
{
    m_x0 = x0;
    m_y0 = y0;
}
void WaveFunctionBuilder::set_deviation(float sigma)
{
    m_sigma = sigma; 
}
auto WaveFunctionBuilder::build_wavefunction(size_t N_y, size_t N_x) const -> Eigen::Reshaped<Eigen::MatrixXcf, -1, 1>
{
    Eigen::MatrixXcf psi_temp = Eigen::MatrixXcf::Zero(N_y-2, N_x-2); 
    init_wavefunction(psi_temp);
    zero_out_boundary(psi_temp);
    return psi_temp.reshaped<Eigen::ColMajor>();
}   
void WaveFunctionBuilder::init_wavefunction(Eigen::MatrixXcf& Psi) const
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
void WaveFunctionBuilder::zero_out_boundary(Eigen::MatrixXcf& Psi) const
{
    Psi.row(0           ).array() = 0;
    Psi.row(Psi.rows()-1).array() = 0;
    Psi.col(0           ).array() = 0;
    Psi.col(Psi.cols()-1).array() = 0;
}
