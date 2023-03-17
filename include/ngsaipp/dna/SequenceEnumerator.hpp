#ifndef NGSAI_DNA_SEQUENCEENUMERATOR_HPP
#define NGSAI_DNA_SEQUENCEENUMERATOR_HPP

#include <string>


namespace ngsai
{   
    namespace dna {
        /*!
        * This class allows to generate all DNA sequences of a given length 
        * (kmer), in lexicographic order.
        * The DNA base lexicographic order is defined in ngsai::dna::base_code 
        * and in ngsai::dna::code_base in dna/dna_utility.hpp.
        */
        class SequenceEnumerator
        {
            public:
                /*!
                * \brief Constructs an empty object. 
                */
                SequenceEnumerator() ;

                /*!
                * \brief Constructs an object to enumerate all sequences of the 
                * given length.
                * \param kmer_size the size of the sequences. 
                */
                SequenceEnumerator(size_t kmer_size) ;

                /*!
                * \brief Copy constructor.
                * \param other the object to copy from. 
                */
                SequenceEnumerator(const ngsai::dna::SequenceEnumerator& other) ;

                /*!
                * \brief Assignment constructor.
                * \param other the object to copy from. 
                */
                SequenceEnumerator(const ngsai::dna::SequenceEnumerator&& other) ;

                /*!
                * \brief Destructor. 
                */
                ~SequenceEnumerator() ;

                /*!
                * \brief Assignment operator.
                * \param other the object to copy from. 
                * \returns a reference to the current instance.
                */
                ngsai::dna::SequenceEnumerator& 
                operator = (const ngsai::dna::SequenceEnumerator& other) ; 

                /*!
                * \brief Move assignment operator.
                * \param other the object to copy from. 
                * \returns a reference to the current instance.
                */
                ngsai::dna::SequenceEnumerator& 
                operator = (const ngsai::dna::SequenceEnumerator&& other) ; 

                /*!
                * \brief Sets the kmer size and reset the enumerator.
                * \param kmer_size the new kmer size. 
                */
                void
                setSize(size_t kmer_size) ;
                
                /*!
                * \brief Gets the next kmer sequence.
                * \param kmer a reference into which the kmer sequence will 
                * be copied.
                * \returns true if a new kmer sequence was enumerated, false if 
                * the enumeration is done. 
                */
                bool
                getNext(std::string& kmer) ;
                
                /*!
                * \brief Resets the enumerator. The next enumerated kmer will be 
                * the first in lexicographic order. 
                */
                void
                reset() ;    

            protected:
                size_t kmer_size ;
                std::string kmer ;
        } ;
    } // namespace dna
}  // namespace ngsai

#endif // NGSAI_DNA_SEQUENCEENUMERATOR_HPP