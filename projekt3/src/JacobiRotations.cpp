#ifndef JACOBIROTATIONS_CPP
#define JACOBIROTATIONS_CPP

#include <vector>
#include <cmath>

namespace JacobiRotations
{
    const double epsilon = 1e-15;
    template <typename T> int sgn(T val) { return (T(0) <= val) - (val < T(0)); }

    inline void findMaxOffdiagonal(std::vector<std::vector<double>>& M, size_t& row, size_t& col, double& max)
    {
        row = 0; col = 1; max = std::abs(M[0][1]);
        for (size_t i = 0; i < M.size(); ++i)
        {
            for (size_t j = i+1; j < M.size(); ++j)
            {
                if (std::abs(M[i][j]) > max)
                {
                    row = i; col = j; max = std::abs(M[i][j]);
                }
            }
        }
    }

    inline void rotate(std::vector<std::vector<double>>& M, size_t k, size_t l, double t)
    {
        double c = 1/std::sqrt(t*t+1);
        double s = c*t;
        double r = (1-c)/s;
        std::vector<double> tmp_k(M.size(), 0);
        std::vector<double> tmp_l(M.size(), 0);

        for (size_t t = 0; t < M.size(); ++t)
        {
            tmp_k[t] = M[k][t] - s*(M[l][t]+r*M[k][t]);
            tmp_l[t] = M[l][t] + s*(M[k][t]-r*M[l][t]);
        }
        tmp_k[k] = M[k][k]-t*M[k][l];
        tmp_k[l] = 0;
        tmp_l[l] = M[l][l]+t*M[k][l];
        tmp_l[k] = 0;
        for (size_t t = 0; t < M.size(); ++t)
        {
            M[k][t] = tmp_k[t];
            M[t][k] = tmp_k[t];
            M[l][t] = tmp_l[t];
            M[t][l] = tmp_l[t];
        }
    }

    inline void diagonalize(std::vector<std::vector<double>>& M)
    {
        std::vector<double> evals(M.size(), 0);
        size_t k,l; double max;
        findMaxOffdiagonal(M, k, l, max);
        while (max > epsilon)
        {
            double beta = (M[l][l]-M[k][k])/(2*M[k][l]);
            double t = sgn(beta)/(abs(beta)+sqrt(beta*beta+1));
            rotate(M,k,l,t);
            findMaxOffdiagonal(M, k, l, max);
        }
    }
};
#endif