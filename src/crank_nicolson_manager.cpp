#include "crank_nicolson_manager.hpp"
 
 
CrankNicolsonManager::CrankNicolsonManager(size_t Nx, size_t Ny)
    :   N_x{Nx}, N_y{Ny}
{
    N_center = (N_y-2)*(N_x-2);
}
void CrankNicolsonManager::set_diagonal_elements(std::complex<float> a0, std::complex<float> b0) // define
{
    m_a0 = a0;
    m_b0 = b0;
}
void CrankNicolsonManager::set_off_diag_elements(std::complex<float> rx, std::complex<float> ry)
{
    m_rx = rx;
    m_ry = ry;
}
auto CrankNicolsonManager::get_sparse_matrices() const -> std::tuple<SparseMatrix,SparseMatrix>
{
    SparseMatrix sparse_A, sparse_M;
    {
        Eigen::MatrixXcf A, M;
        init_temp_matrices(A, M);
        set_temp_matrices(A, M);
        sparse_A = A.sparseView(); 
        sparse_M = M.sparseView(); 
    }
    sparse_A.SparseMatrix::makeCompressed();
    sparse_M.SparseMatrix::makeCompressed();
    return std::tuple(sparse_A, sparse_M);
}
void CrankNicolsonManager::init_temp_matrices(Eigen::MatrixXcf& A, Eigen::MatrixXcf& M) const
{
    A   =  Eigen::MatrixXcf::Zero(N_center, N_center);
    M   =  Eigen::MatrixXcf::Zero(N_center, N_center);
}
void CrankNicolsonManager::set_temp_matrices(Eigen::MatrixXcf& A, Eigen::MatrixXcf& M) const
{
    A.diagonal(0) << Eigen::VectorXcf::Constant(N_center, m_a0);
    M.diagonal(0) << Eigen::VectorXcf::Constant(N_center, m_b0);
    
    // Psi_{i,j} = psi(y,x) = {psi(1,1), psi(2,1), psi(3,1), ... psi(Ny-2,1), psi(1,2) etc} // column major ordering
    size_t iy{}, jx{};
    for (size_t k = 0; k < N_center; k++)
    {
        // k = (j-1)(Ny-2) + (i-1)  // row index iy for y axis, column index jx for x axis
        iy = 1 + k%(N_y-2);          
        jx = 1 + k/(N_y-2);

        // Jumping from (x,y) -> (x,y-1)
        if (iy != 1) 
        {
            A(k, k-1) = -m_ry;
            M(k, k-1) = +m_ry;
        }
        // Jumping from (x,y) -> (x,y+1)
        if (iy != N_y-2) //
        {
            A(k, k+1) = -m_ry;
            M(k, k+1) = +m_ry;
        }
        // Jumping from (x,y) -> (x-1,y)
        if (jx != 1)    
        {
            A(k, k - (N_y-2) ) = -m_rx;
            M(k, k - (N_y-2) ) = +m_rx;
        }
        // Jumping from (x,y) -> (x+1,y)
        if (jx != (N_x-2)) 
        {
            A(k, k + (N_y-2)) = -m_rx;
            M(k, k + (N_y-2)) = +m_rx;
        }
    }
}

