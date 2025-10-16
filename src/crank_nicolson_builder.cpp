#include "crank_nicolson_builder.hpp"
 
using Trplt = Eigen::Triplet<std::complex<float>>;


CrankNicolsonBuilder::CrankNicolsonBuilder()
{
}
void CrankNicolsonBuilder::set_num_elements(size_t Nx, size_t Ny)
{
    m_Nx = Nx;
    m_Ny = Ny;
    N_center = (m_Ny-2)*(m_Nx-2);
}
void CrankNicolsonBuilder::set_diagonal_elements(std::complex<float> a0, std::complex<float> b0) // define
{
    m_a0 = a0;
    m_b0 = b0;
}
void CrankNicolsonBuilder::set_off_diag_elements(std::complex<float> rx, std::complex<float> ry)
{
    m_rx = rx;
    m_ry = ry;
}
auto CrankNicolsonBuilder::get_sparse_matrices() const -> std::tuple<SparseMatrix,SparseMatrix>
{
    SparseMatrix sparse_A, sparse_M;
    {
        Eigen::MatrixXcf dense_A, dense_M;
        init_dense_matrices(dense_A, dense_M);
        set_dense_matrices(dense_A, dense_M);
        sparse_A = std::move(dense_A.sparseView()); 
        sparse_M = std::move(dense_M.sparseView()); 
    }
    sparse_A.SparseMatrix::makeCompressed();
    sparse_M.SparseMatrix::makeCompressed();
    return std::tuple(sparse_A, sparse_M);
}
void CrankNicolsonBuilder::init_dense_matrices(Eigen::MatrixXcf& A, Eigen::MatrixXcf& M) const
{
    A   =  Eigen::MatrixXcf::Zero(N_center, N_center);
    M   =  Eigen::MatrixXcf::Zero(N_center, N_center);
}
void CrankNicolsonBuilder::set_dense_matrices(Eigen::MatrixXcf& A, Eigen::MatrixXcf& M) const
{
    A.diagonal(0) << Eigen::VectorXcf::Constant(N_center, m_a0);
    M.diagonal(0) << Eigen::VectorXcf::Constant(N_center, m_b0);

    // Psi_{i,j} = psi(y,x) = {psi(1,1), psi(2,1), psi(3,1), ... psi(Ny-2,1), psi(1,2) etc} // column major ordering
    size_t iy{}, jx{};
    for (size_t k = 0; k < N_center; k++)
    {
        // k = (j-1)(Ny-2) + (i-1)  // row index iy for y axis, column index jx for x axis
        iy = 1 + k%(m_Ny-2);          
        jx = 1 + k/(m_Ny-2);


        // Jumping from (x,y) -> (x,y-1)
        if (iy != 1) 
        {
            A(k, k-1) = -m_ry;
            M(k, k-1) = +m_ry;
        }
        // Jumping from (x,y) -> (x,y+1)
        if (iy != m_Ny-2) //
        {
            A(k, k+1) = -m_ry;
            M(k, k+1) = +m_ry;
        }
        // Jumping from (x,y) -> (x-1,y)
        if (jx != 1)    
        {
            A(k, k - (m_Ny-2) ) = -m_rx;
            M(k, k - (m_Ny-2) ) = +m_rx;
        }
        // Jumping from (x,y) -> (x+1,y)
        if (jx != (m_Nx-2)) 
        {
            A(k, k + (m_Ny-2)) = -m_rx;
            M(k, k + (m_Ny-2)) = +m_rx;
        }
    }
}

