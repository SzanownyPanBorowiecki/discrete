
#ifndef WILF_COLORING_CPP
#define WILF_COLORING_CPP

#include <iostream>
#include <vector>
#include <algorithm>

#include "Graph.cpp"
#include "misc.cpp"

void wilf_coloring(std::vector<std::vector<bool>> adjacency)
{
    Graph g(adjacency);
 
    // find ordering for greedy coloring
    std::vector<int> ordering(g.size());
    std::vector<bool> removed(g.size(), false);
    Graph h = g;
    for (int i = 0; i < g.size(); ++i)
    {
        int v_ind = h.findMinDegVertex();   // v_ind = vertex id relative to h (induced graph of g)
        int v = 0;                          // v     = vertex id relative to g = v_ind'th non-removed vertex of g
        for (int k = 0; v < g.size(); ++v)
        {
            if (!removed[v])
            {
                if (k == v_ind) break;
                k++;
            }
        }

        ordering[g.size()-1-i] = v;
        h.removeVertex(v_ind);
        removed[v] = true;
    }

    // greedy colouring
    std::vector<int> colouring(g.size(), -1);
    for (int v : ordering)
    {
        std::vector<bool> available_colours(g.size(), true);
        
        // for each coloured neighbour of v, make its colour unavailable
        for (int i = 0; i < g.size(); ++i)
        {
            if (g.adjacency[v][i] && colouring[i] != -1) available_colours[colouring[i]] = false;
        }
        // chose the lowest available colour
        int c = 0;
        for (; !available_colours[c]; ++c);
        
        colouring[v] = c;
    }

    // bounds
    std::vector<double> evals = g.eigenvalues();
    std::pair<std::vector<double>::iterator, std::vector<double>::iterator> minmax = std::minmax_element(evals.begin(), evals.end());
    double eval_min = *minmax.first;
    double eval_max = *minmax.second;
    double lower = (1.0 - eval_max/eval_min);
    double upper = 1.0 + eval_max;

    std::cout << "Kolorowanie: " << colouring << std::endl;
    std::cout << "Przedział: [" << lower << ", " << upper << "]" << std::endl;
}

#endif