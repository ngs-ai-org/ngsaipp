#include <ngsaipp/epigenetics/KmerData.hpp>
#include <iostream>
#include <string>
#include <vector>


ngsai::KmerData::KmerData()
    : sequence(""), ipd(), pwd()
{ ; }


ngsai::KmerData::KmerData(const KmerData& other)
    : sequence(other.sequence),
      ipd(other.ipd),
      pwd(other.pwd)
{ 
    if(this->sequence.size() != this->ipd.size())
    {   char msg[4096] ; 
        sprintf(msg, 
                "sequence length (%zu) and IPD length (%zu) differ",
                this->sequence.size(), this->ipd.size()) ;
        throw std::invalid_argument(msg) ;
    }
    if(this->sequence.size() != this->pwd.size())
    {   char msg[4096] ; 
        sprintf(msg, 
                "sequence length (%zu) and PWD length (%zu) differ",
                this->sequence.size(), this->pwd.size()) ;
        throw std::invalid_argument(msg) ;
    }
}


ngsai::KmerData::KmerData(KmerData&& other)
{   if(other.sequence.size() != other.ipd.size())
    {   char msg[4096] ; 
        sprintf(msg, 
                "sequence length (%zu) and IPD length (%zu) differ",
                other.sequence.size(), other.ipd.size()) ;
        throw std::invalid_argument(msg) ;
    }
    if(other.sequence.size() != other.pwd.size())
    {   char msg[4096] ; 
        sprintf(msg, 
                "sequence length (%zu) and PWD length (%zu) differ",
                other.sequence.size(), other.pwd.size()) ;
        throw std::invalid_argument(msg) ;
    }
    this->sequence = std::move(other.sequence) ;
    this->ipd = std::move(other.ipd) ;
    this->pwd = std::move(other.pwd) ;
}


ngsai::KmerData::KmerData(const std::string& sequence)
    : sequence(sequence),
      ipd(sequence.size(), 0),
      pwd(sequence.size(), 0)
{ ; }


ngsai::KmerData::KmerData(const std::string& sequence,
                          const std::vector<uint32_t>& ipd,
                          const std::vector<uint32_t>& pwd)
    : sequence(sequence), ipd(ipd), pwd(pwd)
{   if(this->sequence.size() != this->ipd.size())
    {   char msg[4096] ; 
        sprintf(msg, 
                "sequence length (%zu) and IPD length (%zu) differ",
                this->sequence.size(), this->ipd.size()) ;
        throw std::invalid_argument(msg) ;
    }
    if(this->sequence.size() != this->pwd.size())
    {   char msg[4096] ; 
        sprintf(msg, 
                "sequence length (%zu) and PWD length (%zu) differ",
                this->sequence.size(), this->pwd.size()) ;
        throw std::invalid_argument(msg) ;
    }
}


ngsai::KmerData::~KmerData()
{ ; }


ngsai::KmerData& 
ngsai::KmerData::operator = (const ngsai::KmerData& other)
{   if(other.sequence.size() != other.ipd.size())
    {   char msg[4096] ; 
        sprintf(msg, 
                "sequence length (%zu) and IPD length (%zu) differ",
                other.sequence.size(), other.ipd.size()) ;
        throw std::invalid_argument(msg) ;
    }
    if(other.sequence.size() != other.pwd.size())
    {   char msg[4096] ; 
        sprintf(msg, 
                "sequence length (%zu) and PWD length (%zu) differ",
                other.sequence.size(), other.pwd.size()) ;
        throw std::invalid_argument(msg) ;
    }

    this->sequence = other.sequence ; 
    this->ipd = other.ipd ;
    this->pwd = other.pwd ;
    return *this ;
}


ngsai::KmerData& 
ngsai::KmerData::operator = (const ngsai::KmerData&& other)
{   if(other.sequence.size() != other.ipd.size())
    {   char msg[4096] ; 
        sprintf(msg, 
                "sequence length (%zu) and IPD length (%zu) differ",
                other.sequence.size(), other.ipd.size()) ;
        throw std::invalid_argument(msg) ;
    }
    if(other.sequence.size() != other.pwd.size())
    {   char msg[4096] ; 
        sprintf(msg, 
                "sequence length (%zu) and PWD length (%zu) differ",
                other.sequence.size(), other.pwd.size()) ;
        throw std::invalid_argument(msg) ;
    }
    
    this->sequence = std::move(other.sequence) ; 
    this->ipd = std::move(other.ipd) ;
    this->pwd = std::move(other.pwd) ;
    return *this ;
}


bool
ngsai::KmerData::operator == (const ngsai::KmerData& other) const
{   if((this->sequence     == other.sequence) and
       (this->ipd          == other.ipd) and
       (this->pwd          == other.pwd))
    {   return true ; }
    else
    {   return false ; }
}


bool
ngsai::KmerData::operator != (const ngsai::KmerData& other) const
{   return not (*this == other) ; }


size_t
ngsai::KmerData::size() const
{   return this->sequence.size() ; }


std::ostream& operator << (std::ostream& stream,
                           const ngsai::KmerData& kmer)
{   stream << "(" 
           << kmer.sequence << " "
           << kmer.ipd << " " 
           << kmer.pwd << ")" ;
    return stream ;
}