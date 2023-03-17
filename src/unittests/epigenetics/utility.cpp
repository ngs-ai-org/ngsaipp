#include "utility.hpp"
#include <cmath>

bool
kmerMapEqual(const ngsai::KmerMap& map1,
             const ngsai::KmerMap& map2)
{   return map1 == map2 ; }


std::vector<std::vector<double>>
operator *(const std::vector<std::vector<double>>& m, double n)
{   std::vector<std::vector<double>> m2(m) ;
    for(auto& v : m2)
    {   for(auto& x : v)
        {   x *= n ; }
    }
    return m2 ;
}


std::vector<std::vector<double>>
operator /(const std::vector<std::vector<double>>& m, double n)
{   return m * (1./ n) ; }


std::vector<std::vector<double>>
log(const std::vector<std::vector<double>>& m)
{   std::vector<std::vector<double>> m2(m) ;
    for(auto& v : m2)
    {   for(auto& x : v)
        {   x = log(x) ; }
    }
    return m2 ;
}