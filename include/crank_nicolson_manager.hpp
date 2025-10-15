

#ifndef CN_MANAGER_HPP
#define CN_MANAGER_HPP

using std::size_t;
#include <cstdint>
#include <eigen3/Eigen/SparseLU>

using SparseMatrix = Eigen::SparseMatrix<std::complex<float>>;

struct CrankNicolsonManager
{
    std::complex<float> a0{};
    std::complex<float> b0{};
    std::complex<float> rx{};
    std::complex<float> ry{};
    size_t N_center{};
    size_t N_x{};
    size_t N_y{};
    CrankNicolsonManager(size_t Nx, size_t Ny)
        :   N_x(Nx), N_y(Ny)
    {
        N_center = (N_y-2)*(N_x-2);
    }
    void set_diagonal_elements(std::complex<float> a0_, std::complex<float> b0_) // define
    {
        a0 = a0_;
        b0 = b0_;
    }
    void set_off_diag_elements(std::complex<float> rx_, std::complex<float> ry_)
    {
        rx = ry_;
        ry = ry_;
    }
    // void set_sparse_matrix(Eigen::SparseMatrix<std::complex<float>>& sparseA, Eigen::SparseMatrix<std::complex<float>>& sparseM)
    // {
    //     Eigen::MatrixXcf A, M;
        
    //     init_temp_matrix(A, M);
    //     set_temp_matrix(A, M);

    //     sparseA = A.sparseView(); 
    //     sparseM = M.sparseView(); 

    //     sparseA.Eigen::SparseMatrix<std::complex<float>>::makeCompressed();
    //     sparseM.Eigen::SparseMatrix<std::complex<float>>::makeCompressed();
    // }
    auto get_sparse_matrices()
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
private:
    void init_temp_matrices(Eigen::MatrixXcf& A, Eigen::MatrixXcf& M)
    {
        A   =  Eigen::MatrixXcf::Zero(N_center, N_center);
        M   =  Eigen::MatrixXcf::Zero(N_center, N_center);
    }
    void set_temp_matrices(Eigen::MatrixXcf& A, Eigen::MatrixXcf& M)
    {
        A.diagonal(0) << Eigen::VectorXcf::Constant(N_center, a0);
        M.diagonal(0) << Eigen::VectorXcf::Constant(N_center, b0);
        
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
                A(k, k-1) = -ry;
                M(k, k-1) = +ry;
            }
            // Jumping from (x,y) -> (x,y+1)
            if (iy != N_y-2) //
            {
                A(k, k+1) = -ry;
                M(k, k+1) = +ry;
            }
            // Jumping from (x,y) -> (x-1,y)
            if (jx != 1)    
            {
                A(k, k - (N_y-2) ) = -rx;
                M(k, k - (N_y-2) ) = +rx;
            }
            // Jumping from (x,y) -> (x+1,y)
            if (jx != (N_x-2)) 
            {
                A(k, k + (N_y-2)) = -rx;
                M(k, k + (N_y-2)) = +rx;
            }
        }
    }
};

#endif