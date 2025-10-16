#include "gaussian_wavefunction_builder.hpp"

constexpr std::complex<float> imaginary_unit{0, 1}; 

void GaussianWfBuilder::set_initial_pos(float x0, float y0)
{
    m_x0 = x0;
    m_y0 = y0;
}
void GaussianWfBuilder::set_system_size(float Lx, float Ly)
{
    m_Lx = Lx;
    m_Ly = Ly;
}
void GaussianWfBuilder::set_deviation(float sigma)
{
    m_sigma = sigma; 
}
auto GaussianWfBuilder::build_wavefunction(size_t N_y, size_t N_x) const -> Eigen::Reshaped<Eigen::MatrixXcf, -1, 1>
{
    Eigen::MatrixXcf psi_temp = Eigen::MatrixXcf::Zero(N_y-2, N_x-2); 
    init_wavefunction(psi_temp);
    zero_out_boundary(psi_temp);
    return psi_temp.reshaped<Eigen::ColMajor>();
}   
void GaussianWfBuilder::init_wavefunction(Eigen::MatrixXcf& Psi) const
{
    Eigen::VectorXf x = Eigen::VectorXf::LinSpaced(Psi.cols(), 0, m_Lx);
    Eigen::VectorXf y = Eigen::VectorXf::LinSpaced(Psi.rows(), 0, m_Ly);
    float delta_x{};
    float delta_y{};
    float gauss_arg{};
    std::complex<float> plane_wave{};

    for (int iy = 0; iy < Psi.rows(); iy++)
    {
        for (int jx = 0; jx < Psi.cols(); jx++)
        {
            delta_x   = (x(jx) - m_x0);
            delta_y   = (y(iy) - m_y0);
            gauss_arg = ( pow(delta_x, 2) + pow(delta_y, 2)) / pow(m_sigma, 2);

            plane_wave  =  std::exp(imaginary_unit*m_k*delta_x);
            Psi(iy,jx) =  std::exp(-0.5f*gauss_arg)*plane_wave; //plane_wave; // std::exp(-0.5*gauss_arg) *
        }
    }
}
void GaussianWfBuilder::zero_out_boundary(Eigen::MatrixXcf& Psi) const
{
    Psi.row(0           ).array() = 0;
    Psi.row(Psi.rows()-1).array() = 0;
    Psi.col(0           ).array() = 0;
    Psi.col(Psi.cols()-1).array() = 0;
}
