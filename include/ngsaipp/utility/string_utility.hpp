#ifndef NGSAI_UTILITY_STRING_UTILITY_HPP
#define NGSAI_UTILITY_STRING_UTILITY_HPP

#include <string>
#include <vector>
#include <sstream>

namespace ngsai
{   
    /*!
     * \brief Splits the string into substrings according to the 
     * given separator.
     * \param string the string to split.
     * \param separator the separator to use.
     * \returns a vector containing the substrings.
     */
    std::vector<std::string> split(const std::string& s, char separator) ;

    /*!
     * \brief Tells whether str ends with the given suffix.
     * \param str the string of interest.
     * \param suffix the suffix of interest.
     * \returns whether str ends with suffix.
     */
    bool endswith(const std::string& str, const std::string& suffix) ;

} // namespace ngsai

#endif // NGSAI_UTILITY_STRING_UTILITY_HPP