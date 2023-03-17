#include <stdexcept>
#include <string>
#include <sstream>    // std::ostringstream

#include <ngsaipp/genome/CpGRegion.hpp>
#include <ngsaipp/genome/constants.hpp>
#include <ngsaipp/genome/GenomeRegion.hpp>



ngsai::genome::CpGRegion::CpGRegion()
    : GenomeRegion()
{ ; }


ngsai::genome::CpGRegion::CpGRegion(const ngsai::genome::GenomeRegion& cpg)
    : GenomeRegion(cpg)
{   this->sanity_check() ; }


ngsai::genome::CpGRegion::CpGRegion(ngsai::genome::GenomeRegion&& cpg)
    : GenomeRegion(cpg)
{   this->sanity_check() ; }


ngsai::genome::CpGRegion::CpGRegion(const std::string& chrom,
                                    size_t start,
                                    size_t end)
    : GenomeRegion(chrom, start, end, ngsai::genome::strand::UNORIENTED)
{   this->sanity_check() ; }


ngsai::genome::CpGRegion::~CpGRegion()
{ ; }


ngsai::genome::GenomeRegion
ngsai::genome::CpGRegion::forwardCoordinates() const
{   return ngsai::genome::GenomeRegion(chrom,
                                       start,
                                       end,
                                       ngsai::genome::strand::FORWARD) ;
}


ngsai::genome::GenomeRegion
ngsai::genome::CpGRegion::reverseCoordinates() const
{   return ngsai::genome::GenomeRegion(chrom,
                                       start,
                                       end,
                                       ngsai::genome::strand::REVERSE) ;
}


void
ngsai::genome::CpGRegion::sanity_check() const
{   
    if(strand != ngsai::genome::strand::UNORIENTED)
    {   char msg[4096] ;
        sprintf(msg, 
                "CpGRegion error : must be unoriented (%c)", 
                ngsai::genome::strand_to_char(strand)) ;
        throw std::invalid_argument(msg) ;
    }
    else if((end - start) != 2)
    {   char msg[4096] ;
        sprintf(msg, 
                "CpGRegion error : forward length must be exactly to 2 (%zu)", 
                end - start) ;
        throw std::invalid_argument(msg) ;
    }
   
}