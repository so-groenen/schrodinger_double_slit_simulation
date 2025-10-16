#ifndef IWF_BUILDER_HPP
#define IWF_BUILDER_HPP
#include <iostream>
#include "Eigen/SparseLU" 

class IWaveFunctionBuilder
{
public:
    virtual void set_system_size(float Lx, float Ly) = 0;
    virtual void set_initial_pos(float x0, float y0) = 0;
    virtual void set_deviation(float sigma) = 0;
    virtual auto build_wavefunction(size_t N_y, size_t N_x) const -> Eigen::Reshaped<Eigen::MatrixXcf, -1, 1> = 0;  
};
 

#endif