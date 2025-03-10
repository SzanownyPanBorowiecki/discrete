#ifndef POLY_CPP
#define POLY_CPP

#include <vector>
#include <string>
#include <iostream>

std::string intToSuperscriptString(int n)
{
    char const *superscript[] = { 
        "\u2070", "\u00B9", "\u00B2", "\u00B3", "\u2074", "\u2075", "\u2076", "\u2077", "\u2078", "\u2079"
    };
    std::string r;
    while (n > 0)
    {
        r.insert(0, superscript[n % 10]);
        n /= 10;
    }
    return r;

}

template<typename T>
class Poly
{
private:
    void fixdeg()
    {
        if (c.size() == 0){ c.resize(1, (T)0); return; }
        if (c.size() == 1) return;

        int i = c.size();
        for (; i >= 2 && c[i-1] == T(0); --i);
        c.resize(i);
    }


public:
    std::vector<T> c;

    Poly(){ c.resize(1, (T)0); };
    Poly(std::vector<T> coefs)
    {
        c = coefs;
        fixdeg();
    }
    Poly(T coef){ c.resize(1); c[0] = coef; }

    template<typename S>
    Poly(Poly<S> &f){*this = f;}

    int deg() const { return (c.size() >= 1) ? c.size()-1 : 0; }

    // this = qb + r
    void div(const Poly<T>& b, Poly<T>& q, Poly<T>& r)
    {
        r = *this;
        q.clear();
        while ((r != (T)0) && r.deg() >= b.deg())
        {
            int tmpSize = r.deg() - b.deg();
            std::vector<T> tmpCoefs(tmpSize+1, (T)0);
            T c1 = r[r.deg()];
            T c2 = b[b.deg()];
            tmpCoefs[tmpSize] = c1/c2;
            Poly<T> t(tmpCoefs);
            q += t;
            r -= t*b;
        }
    }

    void clear()
    {
        c.resize(1);
        c[0] = (T)0;
    }
    
 
    virtual operator int(){ return 0; };


    T operator[](int i) const {return (i <= deg() && c.size() >= i+1) ? c[i] : (T)0;}

    Poly<T>& operator   =(const Poly<T>& x){ c = x.c; return *this; }
    Poly<T>& operator   +=(const Poly<T>& x)
    {
        int d = std::max(deg(), x.deg());
        c.resize(d+1, (T)0);
        for (int i = 0; i <= d; ++i)
        {
            c[i]+=x[i];
        }
        fixdeg();
        return *this;
    }
    Poly<T>& operator   -=(const Poly<T>& x)
    {
        int d = std::max(deg(), x.deg());
        c.resize(d+1, (T)0);
        for (int i = 0; i <= d; ++i)
        {
            c[i]-=x[i];
        }
        fixdeg();
        return *this;        
    }
    Poly<T>& operator   *=(const Poly<T>& x)
    {
        int d = deg()+x.deg();
        std::vector<T> tmp(d+1, (T)0);
        for (int i = 0; i <= d; ++i)
        {
            for (int k = 0; k <= i; ++k)
            {
                tmp[i] += this->operator[](k)*x[i-k];
            }
        }
        c = tmp;
        fixdeg();
        return *this;
    }

    Poly<T>& operator   /=(const Poly<T>& x)
    {
        Poly<T> q, r;
        div(x, q, r);
        *this = q;
        fixdeg();
        return *this;

    }

    Poly<T>& operator   %=(const Poly<T>& x)
    {
        Poly<T> q, r;
        div(x, q, r);
        *this = r;
        fixdeg();
        return *this;
    }


    bool operator       ==(const Poly<T>& x) const
    {
        if (deg() != x.deg()) return false;
        for (int i = 0; i <= deg(); ++i) if (c[i]!=x[i]) return false;
        return true;
    }


    bool operator       !=(const Poly<T>& x) const { return !(*this == x); }
        
    bool operator       ==(const T& x){ return (deg() == 0 && c[0] == x); }
    bool operator       !=(const T& x){ return !(*this == x); }    
    
    Poly<T>& operator   +=(T x){ c[0] += x; return *this; }
    Poly<T>& operator   -=(T x){ c[0] -= x; return *this; }
    Poly<T>& operator   *=(T x){ for (int i = 0; i <= deg(); ++i) c[i] *= x; return *this; }
    Poly<T>& operator   /=(T x){ for (int i = 0; i <= deg(); ++i) c[i] /= x; return *this; }
    Poly<T>& operator   %=(T x){ for (int i = 0; i <= deg(); ++i) c[i] %= x; return *this; }

    friend Poly operator    +(Poly<T> x, const Poly<T>& y){ return x += y; }
    friend Poly operator    -(Poly<T> x, const Poly<T>& y){ return x -= y; }
    friend Poly operator    *(Poly<T> x, const Poly<T>& y){ return x *= y; }
    friend Poly operator    /(Poly<T> x, const Poly<T>& y){ return x /= y; }
    friend Poly operator    %(Poly<T> x, const Poly<T>& y){ return x %= y; }

    friend Poly operator    +(Poly<T> x, const T& y){ return x += y; }
    friend Poly operator    -(Poly<T> x, const T& y){ return x -= y; }
    friend Poly operator    *(Poly<T> x, const T& y){ return x *= y; }
    friend Poly operator    /(Poly<T> x, const T& y){ return x /= y; }
    friend Poly operator    %(Poly<T> x, const T& y){ return x %= y; }
    

    friend std::ostream& operator << (std::ostream& output, const Poly<T>& p)
    {
        output << "[";
        for (int i = p.deg(); i >= 1; i--)
        {
            output << p[i] << "x" << ((i > 1) ? intToSuperscriptString(i) : "") << "+";

        }
        output << p[0] << "]";
        return output;
    }
};

#endif