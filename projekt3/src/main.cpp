#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <bitset>
#include <thread>
#include <chrono>
#include <signal.h>
#include <stdlib.h>
#include <limits.h>
#include <algorithm>

#include "JacobiRotations.cpp"
#include "Graph.cpp"
#include "misc.cpp"
#include "wilf_coloring.cpp"


// spec: https://users.cecs.anu.edu.au/~bdm/data/formats.txt
std::vector<std::vector<bool>> g6strToAdjacencyMatrix(std::string& data)
{
    int size;
    int loc = 0;
    if ((int)data[0] == 126)
    {
        if ((int)data[1] == 126)
        {
            size = (int)data[7]-63 
                + 64*((int)data[6]-63)
                + 4096*((int)data[5]-63) 
                + 262144*((int)data[4]-63)
                + 16777216*((int)data[3]-63)
                + 1073741824*((int)data[2]-63);

            loc = 8;
        }
        else
        {
            size = (int)data[3]-63
                + 64*((int)data[2]-63)
                + 4096*((int)data[1]-63);

            loc = 4;
        }
    }
    else
    {
        size = (int)data[0]-63;
        loc = 1;
    }
    std::vector<std::vector<bool>> adjacency(size, std::vector<bool>(size, 0));
    int row = 0, col = 1;
    for (int i = loc; i < data.size(); ++i)
    {
        std::bitset<CHAR_BIT> bs(data[i]-63);
        for (int j = 0; j < 6 && col < adjacency.size(); ++j)
        {
            adjacency[row][col] = bs[5-j];
            adjacency[col][row] = adjacency[row][col];
            row++;
            if (row == col)
            {
                row = 0; ++col;
            }
        }
    }
    return adjacency;
}

bool kcolouring(Graph& g, int k, std::vector<int>& colouring)
{
    int pos = 0;
    colouring = std::vector<int>(g.size(), 0);

    while (pos >= 0 && pos < g.size())
    {
        std::vector<bool> available_colours(k, true);
        
        for (int i = 0; i < g.size(); ++i)
            if (g.adjacency[pos][i] && colouring[i])
                available_colours[colouring[i]-1] = false;

        int new_colour = colouring[pos]+1;      
        for (; new_colour <= k && !available_colours[new_colour-1]; new_colour++);
        
        if (new_colour == k+1)
        {
            colouring[pos] = 0;
            pos--;
        }
        else
        {
            colouring[pos] = new_colour;
            pos++;
        }
    }
    return (pos >= 0);
}

int binarySearchChromaticNumber(Graph& g, int a, int b)
{
    std::vector<int> colouring;
    while (b-a > 0)
    {
        int m = (a+b)/2;
        if (kcolouring(g,m,colouring)) b = m;
        else a = m+1;
    }
    return a;
}


/* ----------------------------------
    Progress and summary
---------------------------------- */
std::thread t;
int total = 1;
int completed = 0;

double max_diff1 = -1.0;
Graph max_diff1_ex;

double max_diff2 = -1.0;
Graph max_diff2_ex;

int equals_num = 0;
Graph equals_ex;

bool killThread = false;
bool finished = false;

#define PBSTR "========================================"
#define PBWIDTH 40

inline void printStatus()
{
    double percentage = (double)completed/(double)total;
    double val = percentage * 100;
    int lpad = (int) (percentage * PBWIDTH);
    int rpad = PBWIDTH - lpad;
    printf("\r[%.*s>%*s] %6.2f%% (%d/%d)", lpad, PBSTR, rpad, "", val, completed, total);
    fflush(stdout);
}

void statusThread()
{
    while(true)
    {
        if (killThread) return;
        printStatus();
        std::this_thread::sleep_for(std::chrono::milliseconds(50));
    }
}

void printSummary()
{
    std::cout << "\r\e[K" << std::flush;
    std::cout << "Sprawdzono " << completed << " z " << total << " grafów." << std::endl;
    std::cout << std::endl;

    std::cout << "1. Maksymalna różnica między liczbą chromatyczną a dolną granicą (Hoffman): " << max_diff1 << std::endl;
    std::cout << max_diff1_ex.adjacency << std::endl;
    //wilf_coloring(max_diff1_ex.adjacency);
    std::cout << std::endl;
    
    std::cout << "2. Maksymalna różnica między liczbą chromatyczną a górną granicą (Wilff): " << max_diff2 << std::endl;
    std::cout << max_diff2_ex.adjacency << std::endl;
    //wilf_coloring(max_diff2_ex.adjacency);
    std::cout << std::endl;

    std::cout << "3. Liczba grafów gdzie granice są równe: " << equals_num << std::endl;
    std::cout << equals_ex.adjacency << std::endl;
    //wilf_coloring(equals_ex.adjacency);

}

void killStatusThread()
{
    killThread = true;
    t.detach();
    printStatus();
}

void sighandler(int sig)
{
    killStatusThread();
    printSummary();
    exit(0);
}


/* ----------------------------------
    main
---------------------------------- */
int main(int argc, char *argv[])
{
    if (argc < 2){
        std::cout << "Użycie: " << argv[0] << " <nazwa_pliku>" << std::endl;
        return 0;
    }

    std::string filename = argv[1];
    std::cout << "Wczytuję \"" << filename << "\" ... ";

    std::cout.flush();
    std::ifstream ifs(filename, std::ios_base::in);
    if (ifs.fail())
    {
        std::cout << "Błąd!" << std::endl << "Plik nie istnieje?" << std::endl;
        return 0;
    }
    std::vector<std::string> lines;
    for (std::string line; getline(ifs, line); )
        lines.push_back(line);
    std::cout << "OK" << std::endl;
    total = lines.size();
    completed = 0;

    signal(SIGINT, &sighandler);
    t = std::thread(statusThread);
    #pragma omp parallel for
    for (int i = 0; i < lines.size(); ++i)
    {
        std::vector<std::vector<bool>> adjacency = g6strToAdjacencyMatrix(lines[i]);
        Graph G(adjacency);
        
        std::vector<double> evals = G.eigenvalues();
        std::pair<std::vector<double>::iterator, std::vector<double>::iterator> minmax = std::minmax_element(evals.begin(), evals.end());
        double eval_min = *minmax.first;
        double eval_max = *minmax.second;

        double val1 = (1.0 - eval_max/eval_min);
        double val2 = 1.0 + eval_max;

        int chromatic_lowerbound = std::floor(val1);  // Hoffman
        int chromatic_upperbound = std::ceil(val2);   // Wilf
        int chromatic_number = binarySearchChromaticNumber(G, chromatic_lowerbound, chromatic_upperbound);

        double diff1 = (double)chromatic_number - val1;
        double diff2 = (double)val2-chromatic_number;

        #pragma omp critical
        {
            if (diff1 > max_diff1)
            {
                max_diff1 = diff1;
                max_diff1_ex = G;
            }
            if (diff2 > max_diff2)
            {
                max_diff2 = diff2;
                max_diff2_ex = G;
            }
            if (val2 - val1 < 1e-10)
            {
                equals_num++;
                if (equals_ex.complete || equals_ex.size() == 0) equals_ex = G;
            }
        }

        #pragma omp atomic
        completed++;
    }

    finished = true;
    killStatusThread();
    printSummary();

    return 0;
}