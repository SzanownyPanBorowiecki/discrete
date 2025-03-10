#ifndef MOD_CPP
#define MOD_CPP

#include <iostream>

template<typename T>
class mod{
protected:
    T m(T i)
    {
        return i % p;
    }
 
 
    // ax + by = gcd(a,b)
    T extendedEuclid(T a, T b, T& x, T& y)
    {
        if (a == T(0))
        {
            x = T(0), y = T(1);
            return b;
        }
        T x1(0), y1(0);
        T gcd = extendedEuclid(b % a, a, x1, y1);
        x = y1 - (b / a) * x1;
        y = x1;
        return gcd;    
    }
    
public:
    T val;
    static T p;
   
    mod(){val=T(0);}
    mod(T x){val=m(x);}
    ~mod(){};
    
    mod<T> inv()
    {
        T x(0); T y(0);
        T g = extendedEuclid(val, p, x, y);
        x /= g;
        return mod(x);
    }

    operator int() const { return int(val); }

    mod& operator   =(const mod& x){ val = x.val; return *this; }
    mod& operator   +=(const mod& x){ val = m(val+x.val); return *this; }
    mod& operator   -=(const mod& x){ val = m(val-x.val); return *this; }
    mod& operator   *=(const mod& x){ val = m(val*x.val); return *this; }
    mod& operator   /=(mod x){ val = m(val*x.inv().val); return *this; }

    mod& operator   =(T x){ val = m(x); return *this; }
    mod& operator   +=(T x){ val = m(val+x); return *this;}
    mod& operator   -=(T x){ val = m(val-x); return *this;}
    mod& operator   *=(T x){ val = m(val*x); return *this; }
    mod& operator   /=(T x){ return *this *= mod(x).inv(); }

    friend mod operator +(mod x, const mod& y){ return x += y; }
    friend mod operator -(mod x, const mod& y){ return x -= y; }
    friend mod operator *(mod x, const mod& y){ return x *= y; }
    friend mod operator /(mod x, const mod& y){ return x /= y; }
    
    friend mod operator +(mod x, T y){ return x += y; }
    friend mod operator -(mod x, T y){ return x -= y; }
    friend mod operator *(mod x, T y){ return x *= y; }
    friend mod operator /(mod x, T y){ return x /= y; }

    bool operator   ==(const mod& x) const { return val == x.val; }
    bool operator   !=(const mod& x) const { return val != x.val; }
    

    friend std::ostream& operator << (std::ostream& output, const mod& x)
    {
        return output << x.val;
    }
};
template<>
int mod<int>::m(int i){ int v = i%p; return (v < 0) ? v+p : v; }

template<typename T>
T mod<T>::p;

#endif