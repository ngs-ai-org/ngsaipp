#include <ngsaipp/dna/SequenceEnumerator.hpp>
#include <string>
#include <utility>   // std::move()

#include <ngsaipp/dna/dna_utility.hpp>   // ngsai::dna::base_code, ngsai::dna::code_base


ngsai::dna::SequenceEnumerator::SequenceEnumerator()
    : kmer_size(0), kmer()
{ ; }


ngsai::dna::SequenceEnumerator::SequenceEnumerator(size_t kmer_size)
    : kmer_size(kmer_size), kmer()
{ ; }


ngsai::dna::SequenceEnumerator::SequenceEnumerator(
    const ngsai::dna::SequenceEnumerator& other)
    : kmer_size(other.kmer_size), kmer(other.kmer)
{ ; }

ngsai::dna::SequenceEnumerator::SequenceEnumerator(
    const ngsai::dna::SequenceEnumerator&& other)
{   this->kmer_size = std::move(other.kmer_size) ;
    this->kmer      = std::move(other.kmer) ;
}


ngsai::dna::SequenceEnumerator::~SequenceEnumerator()
{ ; }


ngsai::dna::SequenceEnumerator&
ngsai::dna::SequenceEnumerator::operator = (
    const ngsai::dna::SequenceEnumerator& other)
{   this->kmer_size = other.kmer_size ;
    this->kmer = other.kmer ;
    return *this ;
}


ngsai::dna::SequenceEnumerator&
ngsai::dna::SequenceEnumerator::operator = (
    const ngsai::dna::SequenceEnumerator&& other)
{   this->kmer_size = std::move(other.kmer_size) ;
    this->kmer = std::move(other.kmer) ;
    return *this ;
}


void
ngsai::dna::SequenceEnumerator::setSize(size_t kmer_size)
{   this->kmer_size = kmer_size ;
    this->reset() ;
}


bool
ngsai::dna::SequenceEnumerator::getNext(std::string& kmer)
{   
    char smallest_char = ngsai::dna::code_base.at(0) ; // A
    char largest_char = 
        ngsai::dna::code_base.at(ngsai::dna::code_base.size() - 1) ; // T

    // 1st call to method, init the kmer to AA...A
    if(this->kmer.size() == 0)
    {   this->kmer = std::string(this->kmer_size, smallest_char) ;  // A
        kmer = this->kmer ;
        return true ;
    }

    // start from the rightmost side and find the first number less than n
    int p = this->kmer_size - 1;
    while (this->kmer[p] == largest_char)  // T
    {   p--; }
 
    // If all char have max value in kmer then there is no successor, no more
    // kmer
    if (p < 0)
    {   return false ; }
 
    // successor
    this->kmer[p] = 
        ngsai::dna::code_base.at(ngsai::dna::base_code.at(this->kmer[p]) + 1) ;
    for(size_t i=p+1; i<this->kmer_size; i++)
    {   this->kmer[i] = 'A'; }

    kmer = this->kmer ;
    return true ;
}


void
ngsai::dna::SequenceEnumerator::reset()
{   this->kmer = std::string() ; }
