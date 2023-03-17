#include <cmath>
#include <sstream>
#include <stdexcept>

#include <ngsaipp/dna/dna_utility.hpp>
#include <ngsaipp/dna/constants.hpp>
#include <ngsaipp/io/colors.hpp>


std::string
ngsai::dna::get_reverse_complement(const std::string& seq)
{   
    // reverse
    std::string seq_rc(seq.rbegin(), seq.rend());

    // complement
    for(auto& c:seq_rc)
    {   try
        {    c = ngsai::dna::base_pairing.at(c) ; }
        catch(const std::out_of_range& e)
        {   char msg[8096] ; 
            sprintf(msg,
                    "invalid DNA character (%c) in sequence (%s)",
                    c,
                    seq.c_str()) ;
            throw std::invalid_argument(msg) ;
        }
    }
    
    return seq_rc ;
}

std::string
ngsai::dna::colored_dna_seq(const std::string& seq)
{
    std::stringstream seq_col ;
    for(const char& c : seq)
    {   switch (c)
        {
            case 'A':
                seq_col << ngsai::CGREEN << c ;
                break;
            case 'C':
                seq_col << ngsai::CBLUE << c ;
                break;
            case 'G':
                seq_col << ngsai::CYELLOW << c ;
                break;
            case 'T':
                seq_col << ngsai::CRED << c ;
                break;
             case 'N':
                seq_col << ngsai::CWHITE << c ;
                break;
            default:
                seq_col << ngsai::CEND << c ;
                break ;
        }   
    }
    seq_col << ngsai::CEND ;
    return seq_col.str() ;
}


size_t
ngsai::dna::hash_kmer(const std::string& sequence)
{   return ngsai::dna::hash_kmer(sequence, 0, sequence.size()) ; }


size_t
ngsai::dna::hash_kmer(const std::string& sequence,
                 size_t start,
                 size_t length)
{   
    if(start > sequence.size() - length)
    {   char msg[4096] ;
        sprintf(msg,
               "no kmer of size %zu starts at position %zu in "
               "sequence of length %zu",
               length, 
               start, 
               sequence.size()) ;
        throw std::invalid_argument(msg) ;
    }

    size_t hash = 0 ;
    for(size_t i=start, j=0; i<start+length; i++, j++)
    {   try
        {
            hash += (ngsai::dna::base_code.at(sequence[i]) * 
                    pow(ngsai::dna::dna_alphabet_size, length-j-1)) ;
        }
        catch(const std::out_of_range& e)
        {   char msg[8096] ; 
            sprintf(msg,
                    "invalid DNA character (%c) in sequence (%s)",
                    sequence[i],
                    sequence.c_str()) ;
            throw std::invalid_argument(msg) ;
        }
    }
    return hash ;
}


std::vector<size_t>
ngsai::dna::hash_kmer(const std::string& sequence,
                 size_t kmer_size)
{   if(kmer_size > sequence.size())
    {   char msg[4096] ;
        sprintf(msg,
               "kmer size (%zu) bigger than sequence kmer_size (%zu)",
               kmer_size,  
               sequence.size()) ;
        throw std::invalid_argument(msg) ;
    }

    // reserve space in vector
    std::vector<size_t> hashes(sequence.size() - kmer_size + 1) ;

    // 1st kmer hash
    size_t h = ngsai::dna::hash_kmer(sequence, 0, kmer_size) ;
    hashes[0] = h ;

    // following kmer hashes
    for(size_t start=1; start<sequence.size() - kmer_size + 1; start++)
    {   try
        {
            // remove 1st char weight into hash
            h -= (ngsai::dna::base_code.at(sequence[start-1]) * 
                  pow(ngsai::dna::dna_alphabet_size, kmer_size - 1)) ;
            // update power of all remaining char weight
            h *= ngsai::dna::dna_alphabet_size ;
            // add last char weight to hash
            h += ngsai::dna::base_code.at(sequence[start + kmer_size - 1]) ;
            hashes[start] = h ;
        }
        catch(const std::out_of_range& e)
        {   char msg[8096] ; 
            sprintf(msg,
                    "invalid DNA character (%c or %c) in sequence (%s)",
                    sequence[start-1],
                    sequence[start + kmer_size - 1],
                    sequence.c_str()) ;
            throw std::invalid_argument(msg) ;

        }
    }
    return hashes ;
}


std::vector<size_t>
ngsai::dna::hash_kmer(const std::string& sequence,
                 std::vector<std::string>& kmer_sequences,
                 size_t kmer_size)
{   if(kmer_size > sequence.size())
    {   char msg[4096] ;
        sprintf(msg,
               "kmer size (%zu) bigger than sequence length (%zu)",
               kmer_size,  
               sequence.size()) ;
        throw std::invalid_argument(msg) ;
    }
    
    // reserve space in vectors
    // all sequences are init to 'ZZZ' (nb char = kmer size)
    size_t n = sequence.size() - kmer_size + 1;
    kmer_sequences.resize(n, std::string(kmer_size, 'Z')) ;
    std::vector<size_t> hashes(n) ;

    // 1st kmer hash
    size_t h = ngsai::dna::hash_kmer(sequence, 0, kmer_size) ;
    // 1st kmer sequence
    // NOTE could this be optimized?
    for(size_t i=0; i<kmer_size; i++)
    {   kmer_sequences[0][i] = sequence[i] ; }

    hashes[0] = h ;

    // following kmer hashes
    for(size_t start=1; start<n; start++)
    {   try
        {
            // remove 1st char weight into hash
            h -= (ngsai::dna::base_code.at(sequence[start-1]) * 
                  pow(ngsai::dna::dna_alphabet_size, kmer_size - 1)) ;
            // update power of all remaining char weight
            h *= ngsai::dna::dna_alphabet_size ;
            // add last char weight to hash
            h += ngsai::dna::base_code.at(sequence[start + kmer_size - 1]) ;
            hashes[start] = h ;
        }
        catch(const std::out_of_range& e)
        {   char msg[8096] ; 
            sprintf(msg,
                    "invalid DNA character (%c or %c) in sequence (%s)",
                    sequence[start-1],
                    sequence[start + kmer_size - 1],
                    sequence.c_str()) ;
            throw std::invalid_argument(msg) ;
        }

        // kmer sequence
        // NOTE could this be optimized?
        for(size_t i=0; i<kmer_size; i++)
        {   kmer_sequences[start][i] = sequence[start + i] ; }
    }
    return hashes ;
}


std::vector<size_t>
ngsai::dna::hash_kmer(const std::string& sequence,
                 const std::vector<uint32_t>& ipds,
                 const std::vector<uint32_t>& pwds,
                 std::vector<std::string>& kmer_sequences,
                 std::vector<std::vector<uint32_t>>& kmer_ipds,
                 std::vector<std::vector<uint32_t>>& kmer_pwds,
                 size_t kmer_size)
{   if(kmer_size > sequence.size())
    {   char msg[4096] ;
        sprintf(msg,
               "kmer size (%zu) bigger than sequence length (%zu)",
               kmer_size,  
               sequence.size()) ;
        throw std::invalid_argument(msg) ;
    }
    if(sequence.size() != ipds.size())
    {   char msg[4096] ;
        sprintf(msg,
                "sequence length (%zu) and IPD length (%zu) don't match",
                sequence.size(), ipds.size()) ;
        throw std::invalid_argument(msg) ;
    }
    if(sequence.size() != pwds.size())
    {   char msg[4096] ;
        sprintf(msg,
                "sequence length (%zu) and PWD length (%zu) don't match",
                sequence.size(), pwds.size()) ;
        throw std::invalid_argument(msg) ;
    }
    
    // reserve space in vectors
    // all sequences are init to 'ZZZ' (nb char = kmer size)
    // all IPD/PWD are init to 0
    size_t n = sequence.size() - kmer_size + 1;
    kmer_sequences.resize(n, std::string(kmer_size, 'Z')) ;
    kmer_ipds.resize(n, std::vector<uint32_t>(kmer_size, 0)) ;
    kmer_pwds.resize(n, std::vector<uint32_t>(kmer_size, 0)) ;
    std::vector<size_t> hashes(n) ;

    // 1st kmer hash
    size_t h = ngsai::dna::hash_kmer(sequence, 0, kmer_size) ;
    // 1st kmer sequence/IPD/PWD
    // NOTE could this be optimized?
    for(size_t i=0; i<kmer_size; i++)
    {   kmer_sequences[0][i] = sequence[i] ; 
        kmer_ipds[0][i]      = ipds[i] ;
        kmer_pwds[0][i]      = pwds[i] ; 
    }

    hashes[0] = h ;

    // following kmer hashes
    for(size_t start=1; start<n; start++)
    {   try
        {
            // remove 1st char weight into hash
            h -= (ngsai::dna::base_code.at(sequence[start-1]) * 
                  pow(ngsai::dna::dna_alphabet_size, kmer_size - 1)) ;
            // update power of all remaining char weight
            h *= ngsai::dna::dna_alphabet_size ;
            // add last char weight to hash
            h += ngsai::dna::base_code.at(sequence[start + kmer_size - 1]) ;
            hashes[start] = h ;
        }
        catch(const std::out_of_range& e)
        {   char msg[8096] ; 
            sprintf(msg,
                    "invalid DNA character (%c or %c) in sequence (%s)",
                    sequence[start-1],
                    sequence[start + kmer_size - 1],
                    sequence.c_str()) ;
            throw std::invalid_argument(msg) ;
        }

        // kmer sequence/IPD/PWD
        // NOTE could this be optimized?
        for(size_t i=0; i<kmer_size; i++)
        {   kmer_sequences[start][i] = sequence[start + i] ; 
            kmer_ipds[start][i]      = ipds[start + i] ;
            kmer_pwds[start][i]      = pwds[start + i] ; 
        }
    }
    return hashes ;
}