#include "schrodinger_equation.hpp"
  
SchodingerEquation::SchodingerEquation(size_t Nx_, size_t Ny_, const SparseMatrix& sparse_A, SparseMatrix&& sparse_M, Eigen::VectorXcf&& initial_wf)
    : Nx{Nx_}, Ny{Ny_}, m_solver(sparse_A), m_sparse_M{sparse_M}, m_psi{initial_wf}
{
    m_psi_backup = m_psi;
}


float SchodingerEquation::get_wf_modulus(size_t k) const
{
    std::complex<float> z = m_psi(k);
    return sqrt(   pow(z.real(),2) + pow(z.imag(),2)   );
}
void SchodingerEquation::evolve()
{
    m_psi_temp.noalias() = m_sparse_M*m_psi;
    m_psi = m_solver.solve(m_psi_temp);
}
void SchodingerEquation::interact(const Interferometer& double_slit)
{
    double_slit.activate_interaction(m_psi);
}
float SchodingerEquation::get_max_amplitude() const
{
    float vmax{-1};
    float value{};
    for (size_t k = 0; k < (size_t)m_psi.size(); k++)
    {
        value = get_wf_modulus(k);
        if (value > vmax)
        {
            vmax = value;
        }
    }
    return vmax;
}
void SchodingerEquation::reset()
{
    m_psi = m_psi_backup;
}