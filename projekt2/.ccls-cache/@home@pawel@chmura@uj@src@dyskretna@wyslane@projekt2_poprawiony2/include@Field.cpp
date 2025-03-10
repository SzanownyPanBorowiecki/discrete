#ifndef FIELD_CPP
#define FIELD_CPP

#include <vector>
#include <random>
#include "Poly.cpp"
#include "mod.cpp"
#include "utils.cpp"

class Field
{
    typedef mod<int> Zp;
    typedef mod<Poly<Zp>> F;

private:
    static int n;
    static std::mt19937 rng;
    static std::uniform_int_distribution<> unif;

    static void randomizePoly(Poly<Zp> &p)
    {
        for (int i = 0; i <= p.deg()-1; ++i)
        {
            p.c[i] = Zp(unif(rng));
        }
        p.c[p.deg()] = Zp(1);
    }

    static Poly<Zp> randomIrreduciblePoly(int deg)
    {
        Poly<Zp> ret;
        ret.c.resize(deg+1, Zp(0));
        randomizePoly(ret);
        while (!utils::isIrreducible(ret))
        {
            randomizePoly(ret);
        }
        return ret;
    }


public:
    F val;
    static void init(int prime, int exponent)
    {
        //std::cout << "Field::init(" << prime << ", " << exponent << ")" << std::endl;
        Zp::p = prime;
        std::random_device rd;
        rng = std::mt19937(rd());
        unif = std::uniform_int_distribution<>(0, Zp::p-1);
        F::p = randomIrreduciblePoly(exponent);
        //std::cout << "Irreducible polynomial = " << F::p << std::endl;
        n = utils::pow(prime, exponent);
    }
    static int size(){ return n; }

    Field(){};

    Field(int x)
    {
        std::vector<Zp> coefs;
        while (x > 0)
        {
            coefs.push_back(x % Zp::p);
            x /= Zp::p;
        }
        val = F(Poly<Zp>(coefs));
    }

    Field(Poly<Zp> x)
    {
        val = F(x);
    }

    Field(F x)
    {
        val = x;
    }

    Field inv()
    {
        return Field(val.inv());
    }


    operator int() const
    {
        int r = 0;
        int b = 1;
        for (int i = 0; i <= val.val.deg(); ++i)
        {
            r += int(val.val.c[i])*b;
            b *= Zp::p;
        }
        return r;        
    }

    Field& operator   = (const Field& x){ val = x.val; return *this; }
    Field& operator   +=(const Field& x){ val = val+x.val; return *this; }
    Field& operator   -=(const Field& x){ val = val-x.val; return *this; }
    Field& operator   *=(const Field& x){ val = val*x.val; return *this; }
    Field& operator   /=(Field x){ val = val*x.val.inv().val; return *this; }

    friend Field operator +(Field x, const Field& y){ return x += y; }
    friend Field operator -(Field x, const Field& y){ return x -= y; }
    friend Field operator *(Field x, const Field& y){ return x *= y; }
    friend Field operator /(Field x, const Field& y){ return x /= y; }

    bool operator   ==(const Field& x) const { return val == x.val; }
    bool operator   !=(const Field& x) const { return val != x.val; }

    friend std::ostream& operator << (std::ostream& output, const Field& x)
    {
        return output << x.val;
    }
};

std::mt19937 Field::rng;
std::uniform_int_distribution<> Field::unif;
int Field::n;

#endif