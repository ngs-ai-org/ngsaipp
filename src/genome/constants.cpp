#include <stdexcept>

#include <ngsaipp/genome/constants.hpp>


ngsai::genome::strand ngsai::genome::char_to_strand(char strand)
{   switch(strand)
    {   case '+':
        {   return ngsai::genome::strand::FORWARD ; }
        case '-':
        {   return ngsai::genome::strand::REVERSE ; }
        case '.':
        {   return ngsai::genome::strand::UNORIENTED ; }
        default:
        {   char msg[4096] ;
            sprintf(msg,
                    "invalid strand (%c)", 
                    strand) ;
            throw std::invalid_argument(msg) ;
        }
    }
}


char ngsai::genome::strand_to_char(ngsai::genome::strand strand)
{   switch(strand)
    {   case ngsai::genome::strand::FORWARD:
        {   return '+' ; }
        case ngsai::genome::strand::REVERSE:
        {   return '-' ; }
        case ngsai::genome::strand::UNORIENTED:
        {   return '.' ; }
        default:
        {   char msg[4096] ;
            sprintf(msg,
                    "invalid strand (%c)", 
                    strand) ;
            throw std::invalid_argument(msg) ;
        }
    }
}