#ifndef MISC_CPP
#define MISC_CPP
#include <iostream>
#include <vector>

template < class T >
std::ostream& operator << (std::ostream& os, const std::vector<T>& v) 
{
    os << "{";
    for (typename std::vector<T>::const_iterator ii = v.begin(); ii != v.end(); ++ii)
    {
        os << *ii << (ii == v.end()-1 ? "" : ",");
    }
    os << "}";
    return os;
}
#endif