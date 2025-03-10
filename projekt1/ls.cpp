#include <vector>
#include <iostream>
#include <stdio.h>

namespace LS
{
    inline void print(int size, std::vector<int> &src)
    {
        for (int i = 0; i < size; ++i)
        {
            for (int j = 0; j < size; ++j) printf("%2i", src[i*size+j]+1);
            std::cout << std::endl;
        }
    }
    
    inline void print_tex(int size, std::vector<int> &src)
    {
        std::cout << "$\\left[ \\begin{array}{";
        for (int i = 0; i < size; ++i) std::cout << "c";
        std::cout << "} ";
        for (int i = 0; i < size; ++i)
        {
            for (int j = 0; j < size; ++j)
            {
                std::cout << src[i*size + j]+1;
                if (j < size - 1) std::cout << " & ";
                //else std::cout << " "
            }
            if (i < size - 1) std::cout << " \\\\ ";
        }
        std::cout << "\\end{array}\\right]$" << std::endl;
    }

    int numTransversals(int size, std::vector<int> &src, int n, std::vector<bool> &allowed_cols, std::vector<bool> &allowed_vals)
    {
        int r = 0;
        for (int i = 0; i < size; ++i)
        {
            if (!allowed_cols[i] || !allowed_vals[src[n*size + i]]) continue;
            if (n < size-1)
            {
                allowed_cols[i] = false;
                allowed_vals[src[n*size+i]] = false;
                r += numTransversals(size, src, n+1, allowed_cols, allowed_vals);
                allowed_cols[i] = true;
                allowed_vals[src[n*size+i]] = true;                
            }
            else r += 1;
        }
        return r;
    }

    bool orthogonalMate(int size, std::vector<int> &src, int n, std::vector<int> &current, std::vector<std::vector<bool> > &allowed_h, std::vector<std::vector<bool> > &allowed_v, std::vector<std::vector<bool>> &allowed_pairs)
    {
       if (n < size)
       {
           current[n] = src[n];
           allowed_h[0][src[n]] = false;
           allowed_v[n][src[n]] = false;
           allowed_pairs[src[n]][src[n]] = false;
           return (size == 1) || orthogonalMate(size, src, n+1, current, allowed_h, allowed_v, allowed_pairs);
       }
       int row = n / size;
       int col = n - row*size;
       for (int i = 0; i < size; ++i)
       {
           if (!allowed_h[row][i] || !allowed_v[col][i] || !allowed_pairs[src[n]][i]) continue;
           allowed_h[row][i] = false;
           allowed_v[col][i] = false;
           allowed_pairs[src[n]][i] = false;
           current[n] = i;
           if (n == size*size-1 || orthogonalMate(size, src, n+1, current, allowed_h, allowed_v, allowed_pairs)) return true;
           allowed_h[row][i] = true;
           allowed_v[col][i] = true;
           allowed_pairs[src[n]][i] = true;
       }
       return false;
    }

    void iterate(int size, int n, std::vector<int> &current, void (*callback)(int, std::vector<int>&))
    {
        if (n < size*size)
        {
            int row = n / size;
            int col = n - row*size;
            
            if (row == 0 || col == 0)
            {
                current[n]=row+col;
                iterate(size, n+1, current, callback);
            }
            else
            {
                std::vector<bool> allowed(size, true);
                for (int c = col-1; c >= 0; --c)
                    allowed[current[n-1-c]] = false;
                for (int id = n-size; id >= 0; id-=size)
                    allowed[current[id]] = false;
                
                #if defined(_OPENMP)
                #pragma omp parallel for
                #endif
                for (int j = 0; j < size; ++j)
                {
                    if (!allowed[j]) continue;
                    std::vector<int> newseq = current;
                    newseq[n] = j;
                    iterate(size, n+1, newseq, callback);      
                }
            }
        }
        else callback(size, current);
    }
}