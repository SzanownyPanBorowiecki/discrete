#ifndef UTILS_CPP
#define UTILS_CPP

#include "Poly.cpp"
#include "mod.cpp"

namespace utils
{
    typedef mod<int> Zp;

    template<typename T>
    T pow(T a, int p)
    {
        T r = a;
        for (; p>=2; --p)
        {
            r *= a;
        }
        return r;
    }

    Zp gcd(Zp a, Zp b)
    {
        while (b != Zp(0))
        {
            Zp r = a % b;
            a = b;
            b = r;
        }
        return a;    
    }

    Poly<Zp> gcd(Poly<Zp> a, Poly<Zp> b)
    {
        while (b != Zp(0))
        {
            Poly<Zp> r = a % b;
            a = b;
            b = r;
        }
        a /= a[a.deg()];
        return a;
    }

    bool isIrreducible(Poly<Zp> f)
    {
        int m = f.deg();
        if (m <= 1) return true;

        Poly<Zp> x({0, 1});
        Poly<Zp> u = x;
        for (int i = 1; i <= m/2; ++i)
        {
            u = pow(u, Zp::p);
            u = u % f;
            Poly<Zp> d = gcd(f, u - x);
            if (d != Zp(1)) return false;
        }
        return true;
    }
}
#endif