#ifndef NGSAI_IO_UTILITY_HPP
#define NGSAI_IO_UTILITY_HPP

#include <iostream>
#include <vector>
#include <pbcopper/data/CigarOperation.h>


/*!
 * \brief Writes the given vector to the given stream.
 * \param stream the stream to write to.
 * \param v the vector of interest.
 * \returns a reference to the stream after the writing operations.
 */
template<class T>
std::ostream&
operator << (std::ostream& stream, const std::vector<T>& v)
{   for(const auto& x : v)
    {   stream << x << " " ; }
    return stream ;
}


/*!
 * \brief Writes the given pair to the given stream.
 * \param stream the stream to write to.
 * \param p the pair of interest.
 * \returns a reference to the stream after the writing operations.
 */
template<class T1, class T2>
std::ostream&
operator << (std::ostream& stream, const std::pair<T1,T2>& p)
{   stream << p.first << " " << p.second ;
    return stream ;
}

#endif // NGSAI_IO_UTILITY_HPP