#include  <ngsaipp/epigenetics/KmerMap.hpp>

#include <iostream>
#include <string>
#include <vector>
#include <cmath>
#include <utility>

#include <ngsaipp/epigenetics/KmerData.hpp>
#include <ngsaipp/dna/SequenceEnumerator.hpp>
#include <ngsaipp/dna/dna_utility.hpp>    // ngsai::kmer_hash()


ngsai::KmerMap::KmerMap()
    : m_kmer_n(0),
      m_kmer_size(0),
      m_map()
{ ; }


ngsai::KmerMap::KmerMap(size_t kmer_size)
    : KmerMap()
{ this->resize(kmer_size) ; }


ngsai::KmerMap::KmerMap(const ngsai::KmerMap& other)
    : m_kmer_n(other.m_kmer_n),
      m_kmer_size(other.m_kmer_size),
      m_map(other.m_map)
{ ; }


ngsai::KmerMap::KmerMap(ngsai::KmerMap&& other)
    : m_kmer_n(std::move(other.m_kmer_n)),
      m_kmer_size(std::move(other.m_kmer_size)),
      m_map(std::move(other.m_map))
{ ; }


ngsai::KmerMap::~KmerMap()
{ ; }


ngsai::KmerMap&
ngsai::KmerMap::operator = (const ngsai::KmerMap& other)
{   m_kmer_n     = other.m_kmer_n ;
    m_kmer_size  = other.m_kmer_size ;
    m_map        = other.m_map ;
    return *this ;
}


bool
ngsai::KmerMap::operator == (const ngsai::KmerMap& other) const
{
    if(m_kmer_size != other.m_kmer_size)
    {   return false ; }

    auto iter1 = this->begin() ;
    auto iter2 = other.begin() ;
    while(iter1 != this->end())
    {   if(*iter1 != *iter2)
        {   return false ; }
        iter1++ ;
        iter2++ ;
    }
    return true ;
}


bool
ngsai::KmerMap::operator != (const ngsai::KmerMap& other) const
{   return not ((*this) == other) ; }


std::vector<ngsai::kmerBucket>::iterator 
ngsai::KmerMap::begin()
{   return m_map.begin() ; }


std::vector<ngsai::kmerBucket>::const_iterator 
ngsai::KmerMap::begin() const
{   return m_map.begin() ; }


std::vector<ngsai::kmerBucket>::iterator 
ngsai::KmerMap::end()
{   return m_map.end() ; }


std::vector<ngsai::kmerBucket>::const_iterator 
ngsai::KmerMap::end() const
{   return m_map.end() ; }


size_t
ngsai::KmerMap::size() const
{   return m_kmer_n ; }


size_t
ngsai::KmerMap::capacity() const
{   return m_map.size() ; }


void
ngsai::KmerMap::resize(size_t kmer_size)
{   
    m_kmer_n = 0 ;
    m_kmer_size = kmer_size ;
    m_map = std::vector<ngsai::kmerBucket>(
                pow(ngsai::dna::base_code.size(), kmer_size)) ;

    // initialize the sequences, IPD, PWD and nb of occurences of each KmerData
    ngsai::dna::SequenceEnumerator enumerator(m_kmer_size) ;
    std::string kmer ;
    size_t i=0 ;
    while(enumerator.getNext(kmer))
    {   m_map[i].first = 0 ;
        m_map[i].second= ngsai::KmerData(kmer) ;
        i++ ;
    }
    if(i != m_map.size())
    {   char msg[4096] ;
        sprintf(msg,
                "error, all kmer sequences, IPD and PWD could not all be "
                "initialized (only %zu/%zu)",
                i,
                m_map.size()) ;
    }
}


std::vector<ngsai::kmerBucket>::iterator 
ngsai::KmerMap::find(const ngsai::KmerData& kmer)
{
    size_t hash = ngsai::dna::hash_kmer(kmer.sequence) ;
    // kmer has been inserted
    if(m_map[hash].first)
    {   return m_map.begin() + hash ; }
    // nothing insert here yet -> return end
    else
    {   return m_map.end() ; }
}


std::vector<ngsai::kmerBucket>::const_iterator 
ngsai::KmerMap::find(const ngsai::KmerData& kmer) const 
{
    size_t hash = ngsai::dna::hash_kmer(kmer.sequence) ;
    // kmer has been inserted
    if(m_map[hash].first)
    {   return m_map.begin() + hash ; }
    // nothing insert here yet -> return end
    else
    {   return m_map.end() ; }
}


std::vector<std::vector<ngsai::kmerBucket>::iterator> 
ngsai::KmerMap::find(const std::string& sequence)
{   size_t n_kmer = sequence.size() - m_kmer_size + 1 ;
    std::vector<std::vector<ngsai::kmerBucket>::iterator> iters(n_kmer) ;

    std::vector<size_t> hashes = ngsai::dna::hash_kmer(sequence, 
                                                       m_kmer_size) ;
    for(size_t i=0; i<n_kmer; i++)
    {   // this kmer is in the map
        if(m_map[hashes[i]].first)
        {   iters[i] = m_map.begin() + hashes[i] ; }
        // this kmer is not in the map
        else
        {   iters[i] = m_map.end() ; }
    } 
    return iters ;
}


std::vector<std::vector<ngsai::kmerBucket>::const_iterator> 
ngsai::KmerMap::find(const std::string& sequence) const
{   size_t n_kmer = sequence.size() - m_kmer_size + 1 ;
    std::vector<std::vector<ngsai::kmerBucket>::const_iterator> iters(n_kmer) ;

    std::vector<size_t> hashes = ngsai::dna::hash_kmer(sequence,
                                                       m_kmer_size) ;
    for(size_t i=0; i<n_kmer; i++)
    {   // this kmer is in the map
        if(m_map[hashes[i]].first)
        {   iters[i] = m_map.begin() + hashes[i] ; }
        // this kmer is not in the map
        else
        {   iters[i] = m_map.end() ; }
    } 
    return iters ;
}


void 
ngsai::KmerMap::insert(const ngsai::KmerData& kmer,
                       insert_mode mode)
{   if(kmer.sequence.size() != m_kmer_size)
    {   char msg[4096] ;
        sprintf(msg,
                "kmer sequence length (%zu) and kmer map length (%zu) differ",
                kmer.sequence.size(), m_kmer_size) ;
        throw std::invalid_argument(msg) ;
    }
    if(kmer.ipd.size() != m_kmer_size)
    {   char msg[4096] ;
        sprintf(msg,
                "kmer IPD length (%zu) and kmer map length (%zu) differ",
                kmer.ipd.size(), m_kmer_size) ;
        throw std::invalid_argument(msg) ;
    }
    if(kmer.pwd.size() != m_kmer_size)
    {   char msg[4096] ;
        sprintf(msg,
                "kmer PWD length (%zu) and kmer map length (%zu) differ",
                kmer.pwd.size(), m_kmer_size) ;
        throw std::invalid_argument(msg) ;
    }

    size_t hash = ngsai::dna::hash_kmer(kmer.sequence) ;
    // present
    if(m_map[hash].first)
    {   // replace data
        if(mode == insert_mode::REPLACE)
        {   m_map[hash].first      = 1 ; 
            m_map[hash].second.ipd = kmer.ipd ; 
            m_map[hash].second.pwd = kmer.pwd ;
        }
        // add data
        else
        {   m_map[hash].first++ ;
            for(size_t j=0; j<m_kmer_size; j++)
            {   m_map[hash].second.ipd[j] += kmer.ipd[j] ;
                m_map[hash].second.pwd[j] += kmer.pwd[j] ;
            }
        }
    }
    // not present
    else
    {   m_map[hash].first++ ;
        m_map[hash].second = kmer ;
        m_kmer_n++ ;
    }
}


void 
ngsai::KmerMap::insert(const ngsai::KmerData&& kmer,
                       insert_mode mode)
{   if(kmer.sequence.size() != m_kmer_size)
    {   char msg[4096] ;
        sprintf(msg,
                "kmer sequence length (%zu) and kmer map size (%zu) differ",
                kmer.sequence.size(), m_kmer_size) ;
        throw std::invalid_argument(msg) ;
    }
    if(kmer.ipd.size() != m_kmer_size)
    {   char msg[4096] ;
        sprintf(msg,
                "kmer IPD length (%zu) and kmer map size (%zu) differ",
                kmer.ipd.size(), m_kmer_size) ;
        throw std::invalid_argument(msg) ;
    }
    if(kmer.pwd.size() != m_kmer_size)
    {   char msg[4096] ;
        sprintf(msg,
                "kmer PWD length (%zu) and kmer map size (%zu) differ",
                kmer.pwd.size(), m_kmer_size) ;
        throw std::invalid_argument(msg) ;
    }

    size_t hash = ngsai::dna::hash_kmer(kmer.sequence) ;
    // present
    if(m_map[hash].first)
    {   // replace data
        if(mode == insert_mode::REPLACE)
        {   m_map[hash].first      = 1 ;
            m_map[hash].second.ipd = std::move(kmer.ipd) ;
            m_map[hash].second.pwd = std::move(kmer.pwd) ;
        }
        // add data
        else
        {   m_map[hash].first++ ;
            for(size_t j=0; j<m_kmer_size; j++)
            {   m_map[hash].second.ipd[j] += kmer.ipd[j] ;
                m_map[hash].second.pwd[j] += kmer.pwd[j] ;
            }
        }
    }
    // not present
    else
    {   m_map[hash].first++ ;
        m_map[hash].second = std::move(kmer) ;
        m_kmer_n++ ;
    }
}


void 
ngsai::KmerMap::insert(const std::string& sequence,
                       const std::vector<uint32_t>& ipd,
                       const std::vector<uint32_t>& pwd,
                       insert_mode mode)
{   if(sequence.size() != ipd.size())
    {   char msg[4096] ;
        sprintf(msg,
                "given sequence length (%zu) does not match given IPD length "
                "(%zu)",
                sequence.size(), ipd.size()) ;
        throw std::invalid_argument(msg) ;
    }
    if(sequence.size() != pwd.size())
    {   char msg[4096] ;
        sprintf(msg,
                "given sequence length (%zu) does not match given PWD length "
                "(%zu)",
                sequence.size(), pwd.size()) ;
        throw std::invalid_argument(msg) ;
    }
    if(sequence.size() < m_kmer_size)
    {   char msg[4096] ;
        sprintf(msg,
                "given sequence too short (%zu) for current kmer size (%zu)",
                sequence.size(), m_kmer_size) ;
        throw std::invalid_argument(msg) ;
    }

    // hash kmers
    std::vector<std::string> kmer_sequences ;
    std::vector<size_t> kmer_hashes = ngsai::dna::hash_kmer(sequence, 
                                                            kmer_sequences,
                                                            m_kmer_size) ;

    // insert all kmers
    for(size_t i=0; i<kmer_hashes.size(); i++)
    {   std::string kmer_sequence = kmer_sequences[i] ;
        size_t kmer_hash = kmer_hashes[i] ;

        // present
        if(m_map[kmer_hash].first)
        {   // replace data
            if(mode == insert_mode::REPLACE)
            {   m_map[kmer_hash].first = 1 ;
                m_map[kmer_hash].second.ipd = 
                        std::vector(ipd.begin()+i,
                                    ipd.begin()+i + m_kmer_size) ;
                m_map[kmer_hash].second.pwd = 
                        std::vector(pwd.begin()+i,
                                    pwd.begin()+i + m_kmer_size) ;
            }
            // add data
            else
            {   m_map[kmer_hash].first++ ;
                for(size_t j=0; j<m_kmer_size; j++)
                {   m_map[kmer_hash].second.ipd[j] += ipd[i+j] ;
                    m_map[kmer_hash].second.pwd[j] += pwd[i+j] ;
                }
            }
        }
        // not present
        else
        {   m_map[kmer_hash].first++ ;
            m_map[kmer_hash].second.ipd = 
                    std::vector(ipd.begin()+i,
                                ipd.begin()+i + m_kmer_size) ;
            m_map[kmer_hash].second.pwd = 
                    std::vector(pwd.begin()+i,
                                pwd.begin()+i + m_kmer_size) ;
            m_kmer_n++ ;
        }
    }
}


ngsai::kmerBucket& 
ngsai::KmerMap::at(const ngsai::KmerData& kmer)
{   size_t hash = ngsai::dna::hash_kmer(kmer.sequence) ;
    // nothing insert here yet -> exception
    if(m_map[hash].first == 0)
    {   char msg[4096] ;
        sprintf(msg,
                "no kmer with sequence %s in the map",
                kmer.sequence.c_str()) ;
        throw std::out_of_range(msg) ;

    }
    else
    {   return m_map[hash] ; }
}


const ngsai::kmerBucket& 
ngsai::KmerMap::at(const ngsai::KmerData& kmer) const
{   size_t hash = ngsai::dna::hash_kmer(kmer.sequence) ;
    // kmer has been left empty -> nothing insert here yet -> exception
    if(m_map[hash].first == 0)
    {   char msg[4096] ;
        sprintf(msg,
                "no kmer with sequence %s in the map",
                kmer.sequence.c_str()) ;
        throw std::out_of_range(msg) ;

    }
    else
    {   return m_map[hash] ; }
}


size_t
ngsai::KmerMap::getKmerSize() const
{   return m_kmer_size ; }


std::ostream& operator << (std::ostream& stream,
                           const ngsai::KmerMap& m)
{   for(auto iter=m.begin(); iter!=m.end(); iter++)
    {   if(iter->first)
        {   stream << iter->first << " " 
                   << iter->second << std::endl ;
        }
    }
    return stream ;
}