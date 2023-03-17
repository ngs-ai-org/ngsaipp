#include <ngsaipp/utility/string_utility.hpp>
#include <vector>
#include <string>


std::vector<std::string>
ngsai::split(const std::string& s, char separator)
{
    std::stringstream ss(s) ;
    std::string segment ;
    std::vector<std::string> substrings ;
    
    while(std::getline(ss, segment, separator))
    {   // if(segment.size() == 0)
        // {   continue ; }
        substrings.push_back(segment);
    }
    return substrings ;
}


bool 
ngsai::endswith(const std::string& str, const std::string& suffix)
{   if (str.length() < suffix.length())
    {   return false; }
    return str.compare(str.length() - suffix.length(), 
                       suffix.length(), suffix) == 0;
}