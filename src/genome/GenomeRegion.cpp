#include <ngsaipp/genome/GenomeRegion.hpp>
#include <ngsaipp/genome/constants.hpp>     // ngsai::genome::strand
#include <ngsaipp/utility/constants.hpp>    // ngsai::bool_ext

#include <stdexcept>
#include <string>
#include <iostream>
#include <sstream>    // std::ostringstream
#include <utility>    // std::move


ngsai::genome::GenomeRegion::GenomeRegion()
    : chrom(""),
      start(0),
      end(0),
      strand(ngsai::genome::strand::UNORIENTED)
{ ; }


ngsai::genome::GenomeRegion::GenomeRegion(const GenomeRegion& other)
    : chrom(other.chrom),
      start(other.start),
      end(other.end),
      strand(other.strand)
{   // check given values
    this->sanity_check() ;
}


ngsai::genome::GenomeRegion::GenomeRegion(GenomeRegion&& other)
    : chrom(std::move(other.chrom)),
      start(std::move(other.start)),
      end(std::move(other.end)),
      strand(std::move(other.strand))
{   // check given values
    this->sanity_check() ;
}


ngsai::genome::GenomeRegion::GenomeRegion(const std::string& chrom,
                                          size_t start,
                                          size_t end,
                                          ngsai::genome::strand strand)
    : chrom(chrom),
      start(start),
      end(end),
      strand(strand)
{   // check given values
    this->sanity_check() ;
}


ngsai::genome::GenomeRegion::~GenomeRegion()
{ ; }


ngsai::genome::GenomeRegion& 
ngsai::genome::GenomeRegion::operator = (
                            const ngsai::genome::GenomeRegion& other)
{   chrom  = other.chrom ;
    start  = other.start ;
    end    = other.end ;
    strand = other.strand ;
    return *this ; 
}


ngsai::genome::GenomeRegion& 
ngsai::genome::GenomeRegion::operator = (ngsai::genome::GenomeRegion&& other)
{   chrom  = std::move(other.chrom) ;
    start  = std::move(other.start) ;
    end    = std::move(other.end) ;
    strand = std::move(other.strand) ;
    return *this ; 
}


bool 
ngsai::genome::GenomeRegion::operator == (
                            const ngsai::genome::GenomeRegion& other) const
{ 
    if((chrom == other.chrom) and
       (start == other.start) and 
       (end == other.end) and
       (strand == other.strand))
    {   return true ; }
    return false ;
}


bool 
ngsai::genome::GenomeRegion::operator != (
                            const ngsai::genome::GenomeRegion& other) const
{ return not ((*this) == other) ; }


bool 
ngsai::genome::GenomeRegion::operator | (const GenomeRegion& other) const
{   if((chrom != other.chrom) or   // on diff chromosomes
       (other.end-1 < start) or    // other upstream this
       (end-1 < other.start))      // other downstream this
    {   return false ; }
    return true ;
}


ngsai::bool_ext 
ngsai::genome::GenomeRegion::operator < (const GenomeRegion& other) const
{   if(chrom != other.chrom)
    {   return ngsai::bool_ext::Undefined ; }
    else if((chrom == other.chrom) and
            (end-1 < other.start))
    { return ngsai::bool_ext::True ; }
    return ngsai::bool_ext::False ;
}


ngsai::bool_ext 
ngsai::genome::GenomeRegion::operator > (const GenomeRegion& other) const
{   if(chrom != other.chrom)
    {   return ngsai::bool_ext::Undefined ; }
    if((chrom == other.chrom) and
       (other.end-1 < start))
    {   return ngsai::bool_ext::True ; }
    return ngsai::bool_ext::False ;

}


size_t
ngsai::genome::GenomeRegion::size() const
{   return end - start ; }


int
ngsai::genome::GenomeRegion::overlap_len(const GenomeRegion& other) const
{   int len = 0 ;
    if((*this) | other)
    {   // this is contained in other or overlap perfectly other
        if(start >= other.start and end <= other.end)
        {   len = this->size() ; }
        // start of this overlaps end other
        else if((other.start < start) and (other.end-1 >= start))
        {   len = other.end - start ; }
        // other contained in this (perect overlap is handled in first case)
        else if(other.start >= start and other.end <= end)
        {   len = other.size() ; }
        // end of this overlaps start of other (only case left)
        else
        {   len = end - other.start ; }
    }
    return len ;
}


std::string 
ngsai::genome::GenomeRegion::toString() const
{
    std::ostringstream oss ;
    oss << chrom << "\t"
        << start << "\t"
        << end << "\t"
        << ngsai::genome::strand_to_char(strand) ;
    return oss.str() ;
}

void 
ngsai::genome::GenomeRegion::sanity_check() const
{
    if(start >= end)
    {   char msg[4096] ;
        sprintf(msg, 
                "GenomeRegion error : start (%zu) must be smaller than end (%zu)", 
                start, end) ;
        throw std::invalid_argument(msg) ;
    }
}

std::ostream& 
ngsai::genome::operator << (std::ostream& stream,
                            const ngsai::genome::GenomeRegion& record)
{
    stream << record.toString() ;
    return stream ;
}