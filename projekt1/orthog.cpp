#include <iostream>
#include <vector>
#include <string>
#include <thread>
#include <chrono>
#include <signal.h>
#include <stdlib.h>

#include "ls.cpp"

#define TABLE_HEADER "| %10s | %10s |\n"
#define TABLE_DATA "\r| %10i | %10.0f |"
#define TABLE_SEP "+------------+------------+"

using namespace std;

size_t n;
int id = 0;
vector<int> srcEx;
vector<int> orthogEx;
bool result = false;

auto start = chrono::steady_clock::now();
thread t;
bool killThread = false;

inline void printStatus()
{
    auto seconds = chrono::duration_cast<chrono::seconds>(chrono::steady_clock::now() - start).count();
    float rate;
    if (seconds == 0) rate = id;
    else rate = float(id) / seconds;
    printf(TABLE_DATA, id, rate);
    fflush(stdout);
}

void statusThread()
{
    while(true)
    {
        if (killThread) return;
        printStatus();
        this_thread::sleep_for(chrono::milliseconds(250));
    }
}

void printSummary()
{
    cout << "Przeszukano " << id << " zredukowanych kwadratów łacińskich rzędu " << n << "."  << endl;
    if (result)
    {
        cout << "Znaleziono parę ortogonalnych kwadratów łacińskich:" << endl << endl;
        LS::print(n, srcEx);
        cout << endl;
        LS::print(n, orthogEx);
    }
    else
    {
        cout << "Nie znaleziono pary ortogonalnych kwadratów łacińskich" << endl;
    }
}

void killStatusThread()
{
    killThread = true;
    t.detach();
    printStatus();
    cout << endl << TABLE_SEP << endl;
}

void sighandler(int sig)
{
    killStatusThread();
    cout << "Przerywam!" << endl;
    printSummary();
    exit(0);
}

void finish()
{
    killStatusThread();
    printSummary();
}

inline void process(int size, vector<int> &current)
{
    vector<vector<bool>> allowed_h(size, vector<bool>(size, true));
    vector<vector<bool>> allowed_v(size, vector<bool>(size, true));
    vector<vector<bool>> allowed_pairs(size, vector<bool>(size, true));
    vector<int> tmp(size*size, 0);
    bool orthog = LS::orthogonalMate(size, current, 0, tmp, allowed_h, allowed_v, allowed_pairs);
    #if defined(_OPENMP)
    #pragma omp critical
    #endif
    {
        ++id;
        if (orthog)
        {
            result = true;
            srcEx = current;
            orthogEx = tmp;
            finish();
            exit(0);
        }
    }
}

int main(int argc, char *argv[])
{
    if (argc < 2){
        cout << "Użycie: " << argv[0] << " n" << endl;
        return 0;
    }
    signal(SIGINT, &sighandler);
    n = stoi(argv[1]);
    cout << TABLE_SEP << endl;
    printf(TABLE_HEADER, "Razem", "Kw./s");
    cout << TABLE_SEP << endl;
    t = thread(statusThread);
    vector<int> data(n*n);
    LS::iterate(n, 0, data, process);
    finish();
    return 0;
}