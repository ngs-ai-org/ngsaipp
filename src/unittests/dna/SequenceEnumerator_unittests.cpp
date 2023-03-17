#include <gtest/gtest.h>
#include <gmock/gmock.h>
#include <utility>        // std::move()

#include <ngsaipp/dna/SequenceEnumerator.hpp>

/*
 * This file contains tests for the ngsai::dna::SequenceEnumerator class from 
 * src/dna/SequenceEnumerator.cpp
 */


/*!
 * \brief A children class of ngsai::dna::SequenceEnumerator for testing 
 * purposes with getter to the inner data.
 */
class SequenceEnumeratorTest : public ngsai::dna::SequenceEnumerator
{
    public:
        /*!
         * \brief Constructs an empty object. 
         */
        SequenceEnumeratorTest()
            : ngsai::dna::SequenceEnumerator()
        {}

        /*!
         * \brief Constructs an object to enumerate all sequences of the 
         * given length.
         * \param kmer_size the size of the sequences. 
         */
        SequenceEnumeratorTest(size_t kmer_size)
            : ngsai::dna::SequenceEnumerator(kmer_size)
        {}

        /*!
         * \brief Copy constructor.
         * \param other the object to copy from. 
         */
        SequenceEnumeratorTest(const SequenceEnumeratorTest& other)
            : ngsai::dna::SequenceEnumerator(other)
        {}

        /*!
         * \brief Assignment constructor.
         * \param other the object to copy from. 
         */
        SequenceEnumeratorTest(const SequenceEnumeratorTest&& other)
            : ngsai::dna::SequenceEnumerator(other)
        {}

        /*!
         * \brief Destructor. 
         */
        ~SequenceEnumeratorTest()
        {}

        SequenceEnumeratorTest&
        operator = (const SequenceEnumeratorTest& other)
        {   ngsai::dna::SequenceEnumerator::operator = (other) ;
            return *this ;
        }

        SequenceEnumeratorTest&
        operator = (const SequenceEnumeratorTest&& other)
        {   ngsai::dna::SequenceEnumerator::operator = (std::move(other)) ;
            return *this ;
        }

        /*!
         * \brief Returns a reference to the kmer size.
         * \returns a reference to the kmer size.
         */
        size_t& getKmerSize()
        {   return this->kmer_size ; }

        /*!
         * \brief Returns a reference to the current kmer sequence.
         * \returns a reference to the current kmer sequence.
         */
        std::string& getKmer()
        {   return this->kmer ; }
} ;


/*!
 * \brief Constructs a ngsai::dna::SequenceEnumerator and checks that the
 * enumerated sequences match the expectations.
 * \param kmer_size the size of the sequences to enumerate.
 * \param sequences_exp the expected sequences to be enumerated, in 
 * order.
 * \return true if the enumerated sequences match the expectations, false 
 * otherwise.
 */
bool check_sequence_enumeration(size_t kmer_size,
                                const std::vector<std::string>& sequences_exp)
{
    ngsai::dna::SequenceEnumerator e(kmer_size) ;
    std::string sequence ;
    try
    {   for(size_t i=0; e.getNext(sequence); i++)
        {   if(sequence != sequences_exp[i])
            {   return false ; }
        }
    }
    catch(const std::exception& e)
    {
        std::cerr << "something occured during "
                     "SequenceEnumerator::getNext() tests" << std::endl 
                  << e.what() << std::endl ;
        return false ;
    }
    return true ;
}


/*!
 * \brief Checks that the constructors work properly.
 */
TEST(SequenceEnumeratorTest, constructor)
{   
    SequenceEnumeratorTest seq1 ;
    ASSERT_EQ(seq1.getKmerSize(), 0) ;
    ASSERT_EQ(seq1.getKmer(), std::string("")) ;

    std::vector<size_t> kmer_sizes = {1,5,9,13} ;
    for(auto kmer_size : kmer_sizes)
    {   
        SequenceEnumeratorTest seq2(kmer_size) ; 
        ASSERT_EQ(seq2.getKmerSize(), kmer_size) ;
        ASSERT_EQ(seq2.getKmer(), std::string("")) ;

        SequenceEnumeratorTest seq3(seq2) ;
        ASSERT_EQ(seq3.getKmerSize(), kmer_size) ;
        ASSERT_EQ(seq3.getKmer(), std::string("")) ;

        SequenceEnumeratorTest seq4(std::move(seq2)) ;
        ASSERT_EQ(seq4.getKmerSize(), kmer_size) ;
        ASSERT_EQ(seq4.getKmer(), std::string("")) ;
    }
}


/*!
 * \brief Checks that the assignment operators work properly.
 */
TEST(SequenceEnumeratorTest, assignment_operator)
{   size_t kmer_size = 11 ;
    
    SequenceEnumeratorTest s1(kmer_size) ;
    
    SequenceEnumeratorTest s2 = s1 ;
    ASSERT_EQ(s2.getKmerSize(), kmer_size) ;
    ASSERT_EQ(s2.getKmer(),     std::string("")) ;

    SequenceEnumeratorTest s3 = std::move(s1) ;
    ASSERT_EQ(s3.getKmerSize(), kmer_size) ;
    ASSERT_EQ(s3.getKmer(),     std::string("")) ;
}


/*!
 * \brief Checks that the getNext() method works properly.
 */
TEST(SequenceEnumeratorTest, getNext)
{   std::vector<std::string> kmer1 = {"A",
                                      "C",
                                      "G",
                                      "T"} ;
    std::vector<std::string> kmer2 = {"AA",
                                      "AC",
                                      "AG",
                                      "AT",
                                      "CA",
                                      "CC",
                                      "CG",
                                      "CT",
                                      "GA",
                                      "GC",
                                      "GG",
                                      "GT",
                                      "TA",
                                      "TC",
                                      "TG",
                                      "TT"} ;
    std::vector<std::string> kmer3 = {"AAA",
                                      "AAC",
                                      "AAG",
                                      "AAT",
                                      "ACA",
                                      "ACC",
                                      "ACG",
                                      "ACT", 
                                      "AGA",
                                      "AGC",
                                      "AGG",
                                      "AGT", 
                                      "ATA",
                                      "ATC",
                                      "ATG",
                                      "ATT", 
                                      "CAA",
                                      "CAC",
                                      "CAG",
                                      "CAT",
                                      "CCA",
                                      "CCC",
                                      "CCG",
                                      "CCT", 
                                      "CGA",
                                      "CGC",
                                      "CGG",
                                      "CGT", 
                                      "CTA",
                                      "CTC",
                                      "CTG",
                                      "CTT",
                                      "GAA",
                                      "GAC",
                                      "GAG",
                                      "GAT",
                                      "GCA",
                                      "GCC",
                                      "GCG",
                                      "GCT", 
                                      "GGA",
                                      "GGC",
                                      "GGG",
                                      "GGT", 
                                      "GTA",
                                      "GTC",
                                      "GTG",
                                      "GTT",
                                      "TAA",
                                      "TAC",
                                      "TAG",
                                      "TAT",
                                      "TCA",
                                      "TCC",
                                      "TCG",
                                      "TCT", 
                                      "TGA",
                                      "TGC",
                                      "TGG",
                                      "TGT", 
                                      "TTA",
                                      "TTC",
                                      "TTG",
                                      "TTT"} ;

    EXPECT_PRED2(check_sequence_enumeration, 1, kmer1) ;
    EXPECT_PRED2(check_sequence_enumeration, 2, kmer2) ;
    EXPECT_PRED2(check_sequence_enumeration, 3, kmer3) ;
}


/*!
 * \brief Checks that the reset() method works properly.
 */
TEST(SequenceEnumeratorTest, reset)
{   size_t kmer_size = 10 ;
    SequenceEnumeratorTest e(kmer_size) ;
    std::string sequence ;
    e.getNext(sequence) ;
    ASSERT_EQ(sequence, std::string(kmer_size, 'A')) ;

    e.reset() ;
    ASSERT_EQ(e.getKmerSize(), kmer_size) ;
    ASSERT_EQ(e.getKmer(), std::string("")) ;
    e.getNext(sequence) ;
    ASSERT_EQ(sequence, std::string(kmer_size, 'A')) ;
}


/*!
 * \brief Checks that the setSize() method works properly.
 */
TEST(SequenceEnumeratorTest, setSize)
{   size_t kmer_size = 10 ;
    SequenceEnumeratorTest e(kmer_size) ;
    std::string sequence ;
    e.getNext(sequence) ;
    ASSERT_EQ(sequence, std::string(kmer_size, 'A')) ;

    kmer_size = 12 ;
    e.setSize(kmer_size) ;
    ASSERT_EQ(e.getKmer(), std::string("")) ;
    ASSERT_EQ(e.getKmerSize(), kmer_size) ;
    e.getNext(sequence) ;
    ASSERT_EQ(sequence, std::string(kmer_size, 'A')) ;
}