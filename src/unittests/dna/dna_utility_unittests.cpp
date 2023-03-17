#include <gtest/gtest.h>
#include <gmock/gmock.h>
#include <cmath>


#include <ngsaipp/dna/dna_utility.hpp>


using testing::ElementsAre ;
using testing::ContainerEq ;


/*
 * This file contains tests for the functions from 
 * src/dna/dna_utility.cpp
 */

/*!
 * \brief Tests function create proper reverse complements and throw 
 * std::invalid_argument if sequence contains non DNA characters.
 */
TEST(dnaTest, reverse_complement)
{   ASSERT_EQ(ngsai::dna::get_reverse_complement("A"), "T") ;
    ASSERT_EQ(ngsai::dna::get_reverse_complement("C"), "G") ;
    ASSERT_EQ(ngsai::dna::get_reverse_complement("G"), "C") ;
    ASSERT_EQ(ngsai::dna::get_reverse_complement("T"), "A") ;
    ASSERT_EQ(ngsai::dna::get_reverse_complement("N"), "N") ;
    ASSERT_EQ(ngsai::dna::get_reverse_complement("AAA"), "TTT") ;
    ASSERT_EQ(ngsai::dna::get_reverse_complement("CCC"), "GGG") ;
    ASSERT_EQ(ngsai::dna::get_reverse_complement("GGG"), "CCC") ;
    ASSERT_EQ(ngsai::dna::get_reverse_complement("TTT"), "AAA") ;
    ASSERT_EQ(ngsai::dna::get_reverse_complement("NNN"), "NNN") ;
    
    EXPECT_THROW(ngsai::dna::get_reverse_complement("SALUT"), 
                                               std::invalid_argument) ;
}


/*!
 * \brief Tests that the proper hashes are computed and throw 
 * std::invalid_argument is sequence is too short or contains non valid 
 * DNA characters.
 */
TEST(dnaTest, hash_kmer)
{   
    // hash_kmer(const std::string& sequence)
    ASSERT_EQ(ngsai::dna::hash_kmer("AA"), 0) << "hash of AA is not 0" ;
    ASSERT_EQ(ngsai::dna::hash_kmer("AC"), 1) << "hash of AC is not 1" ;
    ASSERT_EQ(ngsai::dna::hash_kmer("AG"), 2) << "hash of AG is not 2" ;
    ASSERT_EQ(ngsai::dna::hash_kmer("AT"), 3) << "hash of AT is not 3" ;
    ASSERT_EQ(ngsai::dna::hash_kmer("CA"), 4) << "hash of CA is not 4" ;
    ASSERT_EQ(ngsai::dna::hash_kmer("CC"), 5) << "hash of CC is not 5" ;
    ASSERT_EQ(ngsai::dna::hash_kmer("CG"), 6) << "hash of CG is not 6" ;
    ASSERT_EQ(ngsai::dna::hash_kmer("CT"), 7) << "hash of CT is not 7" ;
    ASSERT_EQ(ngsai::dna::hash_kmer("GA"), 8) << "hash of GA is not 8" ;
    ASSERT_EQ(ngsai::dna::hash_kmer("GC"), 9) << "hash of GC is not 9" ;
    ASSERT_EQ(ngsai::dna::hash_kmer("GG"), 10) << "hash of GG is not 10" ;
    ASSERT_EQ(ngsai::dna::hash_kmer("GT"), 11) << "hash of GT is not 11" ;
    ASSERT_EQ(ngsai::dna::hash_kmer("TA"), 12) << "hash of TA is not 12" ;
    ASSERT_EQ(ngsai::dna::hash_kmer("TC"), 13) << "hash of TC is not 13" ;
    ASSERT_EQ(ngsai::dna::hash_kmer("TG"), 14) << "hash of TG is not 14" ;
    ASSERT_EQ(ngsai::dna::hash_kmer("TT"), 15) << "hash of TG is not 15" ;
    ASSERT_EQ(ngsai::dna::hash_kmer("AAAAAA"), 0) << "hash of AAAAAA is not 0" ;
    ASSERT_EQ(ngsai::dna::hash_kmer("TTTTTT"), pow(4,6)-1) << 
                                            "hash of TTTTTT is not 4095" ;
    ASSERT_EQ(ngsai::dna::hash_kmer("AAAAAAAAA"), 0) << 
                                            "hash of AAAAAAAAA is not 0" ;
    ASSERT_EQ(ngsai::dna::hash_kmer("TTTTTTTTT"), pow(4,9)-1) << 
                                            "hash of TTTTTTTTT is not 262143" ;
    EXPECT_THROW(ngsai::dna::hash_kmer("salut"), std::invalid_argument) ;

    // some data for next tests
    std::string seq("ACGT") ;
    std::vector<uint32_t> ipd({1,2,3,4}) ;
    std::vector<uint32_t> pwd({5,6,7,8}) ;
    std::vector<size_t> starts = {0,1,2} ;
    std::vector<size_t> hashes_exp = {1, 6, 11} ;
    std::vector<std::string> kmers_exp = {"AC", "CG", "GT"} ;
    std::vector<std::vector<uint32_t>> ipds_exp = {{1,2}, {2,3}, {3,4}} ;
    std::vector<std::vector<uint32_t>> pwds_exp = {{5,6}, {6,7}, {7,8}} ;

    // hash_kmer(const std::string& sequence, size_t start, size_t length)
    for(const auto& start : starts)
    {   ASSERT_EQ(ngsai::dna::hash_kmer(seq, start, 2), hashes_exp[start]) ; }
    EXPECT_THROW(ngsai::dna::hash_kmer(seq, 0, 5),     std::invalid_argument) ;
    EXPECT_THROW(ngsai::dna::hash_kmer(seq, 4, 2),     std::invalid_argument) ;
    EXPECT_THROW(ngsai::dna::hash_kmer("salut", 0, 2), std::invalid_argument) ;

    // hash_kmer(const std::string& sequence, size_t kmer_size)
    auto hashes = ngsai::dna::hash_kmer(seq, 2) ;
    EXPECT_THAT(hashes, ContainerEq(hashes_exp)) ;
    EXPECT_THROW(ngsai::dna::hash_kmer(seq, 5),     std::invalid_argument) ;
    EXPECT_THROW(ngsai::dna::hash_kmer("salut", 2), std::invalid_argument) ;

    // hash_kmer(const std::string& sequence, 
    //           std::vector<std::string>& kmer_sequences, 
    //           size_t kmer_size)
    std::vector<std::string> kmers ;
    hashes = ngsai::dna::hash_kmer(seq, kmers, 2) ;
    EXPECT_THAT(hashes, ContainerEq(hashes_exp)) ;
    EXPECT_THAT(kmers, ContainerEq(kmers_exp)) ;
    EXPECT_THROW(ngsai::dna::hash_kmer(seq, kmers, 5),     std::invalid_argument) ;
    EXPECT_THROW(ngsai::dna::hash_kmer("salut", kmers, 2), std::invalid_argument) ;

    // hash_kmer(const std::string& sequence,
    //           const std::vector<uint32_t>& ipds,
    //           const std::vector<uint32_t>& pwds,
    //           std::vector<std::string>& kmer_sequences,
    //           std::vector<std::vector<uint32_t>>& kmer_ipds,
    //           std::vector<std::vector<uint32_t>>& kmer_pwds,
    //           size_t kmer_size) ;
    std::vector<std::vector<uint32_t>> ipds ;
    std::vector<std::vector<uint32_t>> pwds ; 
    hashes = ngsai::dna::hash_kmer(seq, ipd, pwd, kmers, ipds, pwds, 2) ;
    EXPECT_THAT(hashes, ContainerEq(hashes_exp)) ;
    EXPECT_THAT(kmers, ContainerEq(kmers_exp)) ;
    EXPECT_THAT(ipds, ContainerEq(ipds_exp)) ;
    EXPECT_THAT(pwds, ContainerEq(pwds_exp)) ;
    EXPECT_THROW(ngsai::dna::hash_kmer(seq, ipd, pwd, kmers, ipds, pwds, 5),
                 std::invalid_argument) ;
    EXPECT_THROW(ngsai::dna::hash_kmer("salut", ipd, pwd, kmers, ipds, pwds, 2),
                 std::invalid_argument) ;
    EXPECT_THROW(ngsai::dna::hash_kmer(seq,
                                  std::vector<uint32_t>({1,2,3}),
                                  pwd,
                                  kmers,
                                  ipds,
                                  pwds,
                                  2),
                 std::invalid_argument) ;
    EXPECT_THROW(ngsai::dna::hash_kmer(seq,
                                  ipd,
                                  std::vector<uint32_t>({5,6,7}),
                                  kmers,
                                  ipds,
                                  pwds,
                                  2),
                 std::invalid_argument) ;


}

