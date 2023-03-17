#include <ngsaipp/io/BedRecord.hpp>

#include <stdexcept>
#include <utility>      // std::move()
#include <string>
#include <iostream>
#include <sstream>      // std::ostringstream


ngsai::BedRecord::BedRecord()
    : GenomeRegion(),
      name(""),
      score(0.)
{ ; }


ngsai::BedRecord::BedRecord(const BedRecord& other)
    : GenomeRegion(other.chrom,
                   other.start,
                   other.end, 
                   other.strand),
      name(other.name),
      score(other.score)
{ ; }


ngsai::BedRecord::BedRecord(BedRecord&& other)
    : GenomeRegion(other),
      name(std::move(other.name)),
      score(std::move(other.score))
{ ; }


ngsai::BedRecord::BedRecord(const std::string& chrom,
                            size_t start,
                            size_t end,
                            const std::string& name,
                            double score,
                            ngsai::genome::strand strand)
    : GenomeRegion(chrom, start, end, strand),
      name(name),
      score(score)
{ ; }


ngsai::BedRecord::~BedRecord()
{ ; }


ngsai::BedRecord& 
ngsai::BedRecord::operator = (const ngsai::BedRecord& other)
{   this->chrom   = other.chrom ;
    this->start   = other.start ;
    this->end     = other.end ;
    this->strand  = other.strand ;
    this->name    = other.name ;
    this->score   = other.score ;
    return *this ; 
}


ngsai::BedRecord& 
ngsai::BedRecord::operator = (ngsai::BedRecord&& other)
{   this->chrom  = std::move(other.chrom) ;
    this->start  = std::move(other.start) ;
    this->end    = std::move(other.end) ;
    this->strand = std::move(other.strand) ;
    this->name   = std::move(other.name) ;
    this->score  = std::move(other.score) ;
    return *this ; 
}


bool 
ngsai::BedRecord::operator == (const ngsai::BedRecord& other) const
{   if(ngsai::genome::GenomeRegion::operator==(other) and
       name == other.name and
       score == other.score)
    {   return true ; }
    return false ;
}


bool 
ngsai::BedRecord::operator != (const ngsai::BedRecord& other) const
{ return not ((*this) == other) ; }


std::string 
ngsai::BedRecord::toString() const
{
    std::ostringstream oss ;
    oss << this->chrom << "\t"
        << this->start << "\t"
        << this->end << "\t"
        << this->name << "\t"
        << this->score << "\t"
        << ngsai::genome::strand_to_char(this->strand) ;
    return oss.str() ;
}


std::ostream& 
ngsai::operator << (std::ostream& stream,
                    const ngsai::BedRecord& record)
{
    stream << record.toString() ;
    return stream ;
}