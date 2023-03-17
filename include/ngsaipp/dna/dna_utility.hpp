#ifndef NGSAI_DNA_DNA_UTILITY_HPP
#define NGSAI_DNA_DNA_UTILITY_HPP

#include <string>
#include <vector>

#include <ngsaipp/dna/constants.hpp>


namespace ngsai
{   
    namespace dna {

        /*!
        * \brief Compute the reverse complement of a DNA sequence.
        * \param seq A DNA sequence of interest.
        * \throws std::invalid_argument if the sequence contains a non-DNA 
        * character.
        * \returns the reverse complement of the given DNA sequence.
        */
        std::string get_reverse_complement(const std::string& seq) ;

        /*!
        * \brief Creates a new string containing a formatted 
        * DNA sequence that appears colored when printed on 
        * a terminal. Only DNA characters (A, C, G, T, N) are 
        * colored.
        * \param seq the DNA sequence of interest.
        * \returns a new string containing the original DNA 
        * sequence plus some formatting characters.
        */
        std::string colored_dna_seq(const std::string& seq) ;

        /*!
        * \brief Computes a hash for the given sequence.
        * The hash represents the index at which the sequence would be located 
        * if all possible sequences of the same length were sorted.
        * This function is designed for fixed length squences. The hashes of 
        * sequences of different lengths should not be compared. 
        * \param sequence the kmer to compute the hash from.
        * \throws std::invalid_argument if the sequence contains an invalid DNA 
        * character.
        * \returns the kmer hash.
        */
        size_t hash_kmer(const std::string& sequence) ;

        /*!
        * \brief Computes a hash for the sub-sequence starting at the given
        * position in the sequence and of the given length.
        * The hash represents the index at which the sub-sequence would be located 
        * if all possible sequences of the same length were sorted.
        * This function is designed for fixed length sequences. The hashes of 
        * sequences of different lengths should not be compared. 
        * \param sequence the sequencee containing the sub-sequence of interest.
        * \param start the 0 based position at which the sub-sequence starts in 
        * the sequence.
        * \param length the length of the sub-sequence.
        * \throws std::invalid_argument if the sequence is too short to contain 
        * the given kmer or if it contains an invalid DNA character.
        * \returns the kmer hash.
        */
        size_t hash_kmer(const std::string& sequence,
                        size_t start,
                        size_t length) ;

        /*!
        * \brief Computes the hash each kmer (at every offset) in the sequence
        * using a Rabin-Karp rolling hash implementation.
        * The hash represents the index at which the kmers would be located 
        * if all possible kmers with this length were sorted.
        * This function is designed for fixed length kmers. The hashes of 
        * kmers of different lengths should not be compared. 
        * \param sequence the sequence of interest.
        * \param kmer_size the size of the kmer to hash from the sequence.
        * \throws std::invalid_argument if the sequence is too short to contain 
        * the given kmer or if it contains an invalid DNA character.
        * \returns a vector containing the hashes of each kmer, starting at each 
        * offset. The i-th hash corresponds to the kmer starting at position i 
        * in the sequence.
        */
        std::vector<size_t> hash_kmer(const std::string& sequence,
                                    size_t kmer_size) ;

        /*!
        * \brief Computes the hash each kmer (at every offset) in the sequence
        * using a Rabin-Karp rolling hash implementation. The kmer hashes are 
        * returned and the kmer sequences are stored in the given vector.
        * The hash represents the index at which the kmers would be located 
        * if all possible kmers with this length were sorted.
        * This function is designed for fixed length kmers. The hashes of 
        * kmers of different lengths should not be compared. 
        * \param sequence the sequence of interest.
        * \param kmer_size the size of the kmer to hash from the sequence.
        * \param kmer_sequences a vector in which each kmer sequence will be 
        * stored. The sequence order will match the hash order.
        * \throws std::invalid_argument if the sequence is too short to contain 
        * the given kmer or if it contains an invalid DNA character.
        * \returns a vector containing the hashes of each kmer, starting at each 
        * offset. The i-th hash corresponds to the kmer starting at position i 
        * in the sequence.
        */
        std::vector<size_t> hash_kmer(const std::string& sequence,
                                    std::vector<std::string>& kmer_sequences,
                                    size_t kmer_size) ;
        
        /*!
        * \brief Does the same as the previous function but also stored the 
        * IPDs/PWDs corresponding to each kmer in the corresponding vectors.
        * \param sequence the sequence of interest.
        * \param ipds the IPD values corresponding to each sequence position.
        * \param pwds the PWD values corresponding to each sequence position.
        * \param kmer_size the size of the kmer to hash from the sequence.
        * \param kmer_sequences a vector in which each kmer sequence will be 
        * stored. The sequence order will match the hash order.
        * \param kmer_ipds a vector in which each kmer IPD values will be 
        * stored. The IPD order will match the hash order.
        * \param kmer_pwds a vector in which each kmer PWD values will be 
        * stored. The PWD order will match the hash order.
        * \throws std::invalid_argument if the sequence is too short to contain 
        * the given kmer, if it contains an invalid DNA character or if the 
        * sequence, IPDs and PWDs length do not match.
        * \returns a vector containing the hashes of each kmer, starting at each 
        * offset. The i-th hash corresponds to the kmer starting at position i 
        * in the sequence.
        */
        std::vector<size_t> hash_kmer(const std::string& sequence,
                                    const std::vector<uint32_t>& ipds,
                                    const std::vector<uint32_t>& pwds,
                                    std::vector<std::string>& kmer_sequences,
                                    std::vector<std::vector<uint32_t>>& kmer_ipds,
                                    std::vector<std::vector<uint32_t>>& kmer_pwds,
                                    size_t kmer_size) ;

    } // namespace dna
}  // namespace ngsai

#endif // NGSAI_DNA_DNA_UTILITY_HPP