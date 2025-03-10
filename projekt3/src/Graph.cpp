#ifndef GRAPH_CPP
#define GRAPH_CPP

#include <vector>
#include <numeric>
#include "JacobiRotations.cpp"

struct Graph
{
    std::vector<std::vector<bool>> adjacency;
    bool complete;

    Graph() : complete(false) {};
    Graph(std::vector<std::vector<bool>>& m_adjacency)
    {
        adjacency = std::vector<std::vector<bool>>(m_adjacency.size(), std::vector<bool>(m_adjacency.size(), false));
        complete = true;
        for (int i = 0; i < size(); ++i)
        {
            for (int j = i+1; j < size(); ++j)
            {
                adjacency[i][j] = m_adjacency[i][j];
                adjacency[j][i] = adjacency[i][j];
                if (!adjacency[i][j]) complete = false;
            }
        }        
    }

    std::vector<double> eigenvalues()
    {
        std::vector<std::vector<double>> adjacency_double(size());
        for (int r = 0; r < size(); ++r)
            adjacency_double[r] = std::vector<double>(adjacency[r].begin(), adjacency[r].end());
        JacobiRotations::diagonalize(adjacency_double);

        std::vector<double> evals(adjacency_double.size());
        for (size_t i = 0; i < adjacency_double.size(); ++i)
		    evals[i] = adjacency_double[i][i];

        return evals;
    }

    int size() const { return adjacency.size(); }

    int d(int i)
    {
        return std::accumulate(adjacency[i].begin(), adjacency[i].end(), 0);
    }

    int findMinDegVertex()
    {
        int r = 0; int minDeg = d(0);
        for (int v = 1; v < size(); ++v)
        {
            int d_v = d(v);
            if (d_v < minDeg)
            {
                r = v;
                minDeg = d_v;
            }
        }
        return r;
    }

    void removeVertex(int i)
    {
        adjacency.erase(adjacency.begin() + i);
        for (int k = 0; k < size(); ++k)
        {
            adjacency[k].erase(adjacency[k].begin() + i);
        }
    }
};

#endif