#include <iostream>
#include <vector>
#include <string>
#include <thread>
#include <chrono>
#include <signal.h>
#include <stdlib.h>

#include "ls.cpp"

#define TABLE_HEADER "| %10s | %10s | %10s | %10s |\n"
#define TABLE_DATA "\r| %10i | %10.0f | %10i | %10i |"
#define TABLE_SEP "+------------+------------+------------+------------+"

using namespace std;

size_t n;
int id = 0;
int minTrans = -1;
vector<int> minEx;
int maxTrans = -1;
vector<int> maxEx;
auto start = chrono::steady_clock::now();
thread t;
bool killThread = false;
bool finished = false;

inline void printStatus()
{
    auto seconds = chrono::duration_cast<chrono::seconds>(chrono::steady_clock::now() - start).count();
    float rate;
    if (seconds == 0) rate = id;
    else rate = (float(id) / seconds);
    printf(TABLE_DATA, id, rate, minTrans, maxTrans);
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
    cout << "Przeszukano " << id << " zredukowanych kwadratów łacińskich rzędu " << n << "." << endl << endl;
    cout << "t(" << n << ")" << ((finished || minTrans == 0) ? " = " : " <= ") << minTrans << endl;
    cout << "Kwadrat zawierający " << minTrans << " transwersal: " << endl;
    LS::print(n, minEx);
    cout << endl;
    cout << "T(" << n << ")" << (finished ? " = " : " >= ") << maxTrans << endl;
    cout << "Kwadrat zawierający " << maxTrans << " transwersal: " << endl;
    LS::print(n, maxEx);
    cout << endl;
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

inline void process(int size, vector<int> &current)
{
    vector<bool> allowed_cols(size, true);
    vector<bool> allowed_vals(size, true);
    int numTrans = LS::numTransversals(size, current, 0, allowed_cols, allowed_vals);
    vector<int> tmp(size*size);
    vector<vector<bool>> partition_membership(size, vector<bool>(size, false));
    #if defined(_OPENMP)
    #pragma omp critical
    #endif
    {
        if (numTrans > maxTrans){
            maxTrans = numTrans;
            maxEx = current;
        }
        if (numTrans < minTrans || minTrans == -1)
        {
            minTrans = numTrans;
            minEx = current;
        }
        ++id;
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
    printf(TABLE_HEADER, "Razem", "Kw./s", "t(n)", "T(n)");
    cout << TABLE_SEP << endl;
    t = thread(statusThread);
    vector<int> data(n*n);
    LS::iterate(n, 0, data, process);
    finished = true;
    killStatusThread();
    printSummary();
    return 0;
}