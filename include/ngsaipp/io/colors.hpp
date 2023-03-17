#ifndef NGSAI_IO_COLORS_HPP
#define NGSAI_IO_COLORS_HPP


#include <string>

/*!
 *\brief This file contains values that are useful for
 * formatting terminal output.
 */
namespace ngsai
{
    std::string CEND      = "\33[0m" ;  //  end formatting
    std::string CBOLD     = "\33[1m" ;  //  start bold
    std::string CITALIC   = "\33[3m" ;  //  start italic
    std::string CURL      = "\33[4m" ;  //  start URL
    std::string CBLINK    = "\33[5m" ;  //  start blink
    std::string CBLINK2   = "\33[6m" ;  //  start blink2
    std::string CSELECTED = "\33[7m" ;  //  start selected

    std::string CBLACK    = "\33[30m" ;   //  start black color
    std::string CRED      = "\33[31m" ;   //  start red color
    std::string CGREEN    = "\33[32m" ;   //  start green color
    std::string CYELLOW   = "\33[33m" ;   //  start yellow color
    std::string CBLUE     = "\33[34m" ;   //  start blue color
    std::string CVIOLET   = "\33[35m" ;   //  start violet color
    std::string CBEIGE    = "\33[36m" ;   //  start beige color
    std::string CWHITE    = "\33[37m" ;   //  start white color

    std::string CBLACKBG  = "\33[40m" ;  // start black background
    std::string CREDBG    = "\33[41m" ;  // start red background
    std::string CGREENBG  = "\33[42m" ;  // start green background
    std::string CYELLOWBG = "\33[43m" ;  // start yellow background
    std::string CBLUEBG   = "\33[44m" ;  // start blue background
    std::string CVIOLETBG = "\33[45m" ;  // start violet background
    std::string CBEIGEBG  = "\33[46m" ;  // start beige background
    std::string CWHITEBG  = "\33[47m" ;  // start white background

    std::string CGREY     = "\33[90m" ;   //  start grey color
    std::string CRED2     = "\33[91m" ;   //  start red color
    std::string CGREEN2   = "\33[92m" ;   //  start green color
    std::string CYELLOW2  = "\33[93m" ;   //  start yellow color
    std::string CBLUE2    = "\33[94m" ;   //  start blue color
    std::string CVIOLET2  = "\33[95m" ;   //  start violet color
    std::string CBEIGE2   = "\33[96m" ;   //  start beige color
    std::string CWHITE2   = "\33[97m" ;   //  start white color

    std::string CGREYBG    = "\33[100m" ;  //  start beige background
    std::string CREDBG2    = "\33[101m" ;  //  start red background
    std::string CGREENBG2  = "\33[102m" ;  //  start green background
    std::string CYELLOWBG2 = "\33[103m" ;  //  start yellow background
    std::string CBLUEBG2   = "\33[104m" ;  //  start blue background
    std::string CVIOLETBG2 = "\33[105m" ;  //  start violet background
    std::string CBEIGEBG2  = "\33[106m" ;  //  start beige background
    std::string CWHITEBG2  = "\33[107m" ;  //  start white background

} // namespace ngsai

#endif // NGSAI_IO_COLORS_HPP