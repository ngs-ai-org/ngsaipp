#include <gtest/gtest.h>
#include <gmock/gmock.h>
#include <vector>
#include <string>
#include <cmath>
#include <utility>


#include <ngsaipp/epigenetics/KmerMap.hpp>
#include <ngsaipp/dna/dna_utility.hpp>
#include <ngsaipp/dna/SequenceEnumerator.hpp>

using testing::ContainerEq ;


/*
 * This file contains tests for the ngsai::KmerMap class from 
 * src/epigenetics/KmerMap.cpp
 */

/*!
 * \brief A wrapper around KmerMap to access the inner data for testing
 * purposes.
 */
class KmerMapTest : public ngsai::KmerMap
{
    public:
        KmerMapTest() = delete ;

        KmerMapTest(size_t kmer_size)
            : ngsai::KmerMap(kmer_size)
        { ; }

        ~KmerMapTest()
        { ; }

        std::vector<ngsai::kmerBucket>& getMap() 
        {   return m_map ; }

        const std::vector<ngsai::kmerBucket>& getMapConst() const
        {   return m_map ; }
} ;


/*!
 * \brief Returns true if both maps have strictly identical field values.
 * \param map1 the 1st map to compare.
 * \param map2 the 2nd map to compare.
 * \return true if both maps have the same field values, false otherwise.
 */
bool
kmerMapEqual(const KmerMapTest& map1,
             const KmerMapTest& map2)
{   if((map1.size()              == map2.size()) and
       (map1.capacity()          == map2.capacity()) and
       (map1.getMapConst()       == map2.getMapConst()))
    {   return true ; }
    else
    {   return false ; }
}


/*!
 *\brief Tests that the constructor works properly.
 */
TEST(KmerMapTest, constructor)
{
    std::vector<size_t> kmer_sizes = {1,5,9,11} ;
    for(auto kmer_size : kmer_sizes)
    {   // construct kmer map
        KmerMapTest map(kmer_size) ;
        // construct expected data
        std::vector<ngsai::kmerBucket> map_exp(
                map.capacity(),
                std::make_pair(uint32_t(0),ngsai::KmerData())) ;
        ngsai::dna::SequenceEnumerator seqen(kmer_size) ;
        std::string seq ;
        for(size_t i=0 ; seqen.getNext(seq); i++)
        {   map_exp[i].second.sequence = seq ;
            map_exp[i].second.ipd = std::vector<uint32_t>(kmer_size, 0) ;
            map_exp[i].second.pwd = std::vector<uint32_t>(kmer_size, 0) ;
        }
        // test
        ASSERT_EQ(map.size(), 0) ;
        ASSERT_EQ(map.capacity(), 
                  pow(ngsai::dna::dna_alphabet_size, kmer_size)) ;
        // cannot understand why it does not work
        // EXPECT_THAT(map.getMapConst(), ContainerEq(map_exp)) ;
        ASSERT_EQ(map.getMapConst() == map_exp, true) ;
    }
}


/*!
 * \brief Tests that the copy constructor works properly.
 */
TEST(KmerMapTest, constructor_copy)
{
    std::vector<size_t> kmer_sizes = {5, 7, 9, 11} ;
    for(auto kmer_size : kmer_sizes)
    {   KmerMapTest map1(kmer_size) ;
        KmerMapTest map2(map1) ;
        EXPECT_PRED2(kmerMapEqual, map1, map2) ;
    }

}


/*!
 * \brief Tests that the move constructor works properly.
 */
TEST(KmerMapTest, constructor_move)
{
    std::vector<size_t> kmer_sizes = {5, 7, 9, 11} ;
    for(auto kmer_size : kmer_sizes)
    {   KmerMapTest map1(kmer_size) ;
        KmerMapTest map2(map1) ;
        KmerMapTest map3(std::move(map2)) ;
        EXPECT_PRED2(kmerMapEqual, map1, map3) ;
    }

}


/*!
 *\brief Tests that the assignment operator works properly.
 */
TEST(KmerMapTest, assignment_operator)
{
    std::vector<size_t> kmer_sizes = {5, 7, 9, 11} ;
    for(auto kmer_size : kmer_sizes)
    {   KmerMapTest map1(kmer_size) ;
        KmerMapTest map2 = map1 ;
        EXPECT_PRED2(kmerMapEqual, map1, map2) ;
    }
}

/*!
 *\brief Tests that the begin() method works properly.
 */
TEST(KmerMapTest, begin_method)
{   KmerMapTest map(5) ; 
    const KmerMapTest cmap = map ;
    auto  iter =  map.begin() ;
    auto citer = cmap.begin() ;
    ASSERT_EQ(citer, cmap.getMapConst().begin()) ;
    ASSERT_EQ(iter,   map.getMap().begin()) ;
}


/*!
 *\brief Tests that the end() method works properly.
 */
TEST(KmerMapTest, end_method)
{   KmerMapTest map(5) ; 
    const KmerMapTest cmap = map ;
    auto  iter =  map.end() ;
    auto citer = cmap.end() ;
    ASSERT_EQ(citer, cmap.getMapConst().end()) ;
    ASSERT_EQ(iter,   map.getMap().end()) ;
}


/*!
 *\brief Tests that the insert() and find() methods work properly.
 */
TEST(KmerMapTest, insert_find_methods)
{
    // empty map
    size_t kmer_size = 5 ;
    KmerMapTest map1(kmer_size) ;
    
    std::string seq_acgta("ACGTA") ;
    std::string seq_ggggg("GGGGG") ;
    std::vector<uint32_t> ipd_0({0,0,0,0,0}) ;
    std::vector<uint32_t> ipd_1({1,1,1,1,1}) ;
    std::vector<uint32_t> ipd_2({2,2,2,2,2}) ;

    std::vector<uint32_t> pwd_0({0,0,0,0,0}) ;
    std::vector<uint32_t> pwd_2({2,2,2,2,2}) ;
    std::vector<uint32_t> pwd_4({4,4,4,4,4}) ;

    // a kmer with no IPD/PWD
    ngsai::KmerData acgta0(seq_acgta, ipd_0, pwd_0) ;
    // with same seq but different IPD/PWD
    ngsai::KmerData acgta1(seq_acgta, ipd_1, pwd_2) ;
    // with same seq and 2x IPD/PWD
    ngsai::KmerData acgta2(seq_acgta, ipd_2, pwd_4) ;

    // a kmer with no IPD/PWD
    ngsai::KmerData ggggg0(seq_ggggg, ipd_0, pwd_0) ;
    // with same seq but different IPD/PWD
    ngsai::KmerData ggggg1(seq_ggggg, ipd_1, pwd_2) ;
    // with same seq and 2x IPD/PWD
    ngsai::KmerData ggggg2(seq_ggggg, ipd_2, pwd_4) ;

    // insert new and find KmerData
    ASSERT_NO_THROW(map1.insert(acgta0)) ;
    auto iter1 = map1.find(acgta0) ;
    ASSERT_EQ(map1.size(), 1) ;
    ASSERT_EQ(*iter1 == std::make_pair(uint32_t(1), acgta0), true) ;
    iter1 = map1.find(acgta1) ; // IPD/PWD should not influence find()
    ASSERT_EQ(*iter1 == std::make_pair(uint32_t(1), acgta0), true) ;
    // insert replace and find KmerData
    ASSERT_NO_THROW(map1.insert(acgta1, ngsai::KmerMap::insert_mode::REPLACE)) ;
    ASSERT_EQ(map1.size(), 1) ;
    iter1 = map1.find(acgta1) ;
    ASSERT_EQ(*iter1 == std::make_pair(uint32_t(1), acgta1), true) ;
    iter1 = map1.find(acgta0) ;  // IPD/PWD should not influence find()
    ASSERT_EQ(*iter1 == std::make_pair(uint32_t(1), acgta1), true) ;
    // insert sum and find KmerData
    ASSERT_NO_THROW(map1.insert(acgta1, ngsai::KmerMap::insert_mode::SUM)) ;
    ASSERT_EQ(map1.size(), 1) ;
    iter1 = map1.find(acgta1) ;
    ASSERT_EQ(*iter1 == std::make_pair(uint32_t(2), 
                                       ngsai::KmerData(acgta2.sequence,
                                                       acgta2.ipd,
                                                       acgta2.pwd)),
              true) ;
    iter1 = map1.find(acgta0) ;  // IPD/PWD should not influence find()
    ASSERT_EQ(*iter1 == std::make_pair(uint32_t(2), 
                                       ngsai::KmerData(acgta2.sequence,
                                                       acgta2.ipd,
                                                       acgta2.pwd)),
              true) ;
    
    // mv insert new and find KmerData
    ASSERT_NO_THROW(map1.insert(ngsai::KmerData(ggggg0),
                                ngsai::KmerMap::insert_mode::REPLACE)) ;
    ASSERT_EQ(map1.size(), 2) ;
    auto iter2 = map1.find(ggggg0) ;
    ASSERT_EQ(*iter2 == std::make_pair(uint32_t(1), ggggg0), true) ;
    iter2 = map1.find(ggggg1) ;  // IPD/PWD should not influence find()
    ASSERT_EQ(*iter2 == std::make_pair(uint32_t(1), ggggg0), true) ;
    // mv insert replace and find KmerData
    ASSERT_NO_THROW(map1.insert(ngsai::KmerData(ggggg1), 
                                ngsai::KmerMap::insert_mode::REPLACE)) ;
    ASSERT_EQ(map1.size(), 2) ;
    iter2 = map1.find(ggggg1) ;
    ASSERT_EQ(*iter2 == std::make_pair(uint32_t(1), ggggg1), true) ;
    iter2 = map1.find(ggggg0) ;   // IPD/PWD should not influence find()
    ASSERT_EQ(*iter2 == std::make_pair(uint32_t(1), ggggg1), true) ;
    
    // mv insert sum and find KmerData
    ASSERT_NO_THROW(map1.insert(ngsai::KmerData(ggggg1), 
                                ngsai::KmerMap::insert_mode::SUM)) ;
    ASSERT_EQ(map1.size(), 2) ;
    iter2 = map1.find(ggggg1) ;
    ASSERT_EQ(*iter2 == std::make_pair(uint32_t(2), ggggg2), true) ;
    iter2 = map1.find(ggggg0) ;   // IPD/PWD should not influence find()
    ASSERT_EQ(*iter2 == std::make_pair(uint32_t(2), ggggg2), true) ;
    

    // insert replace from sequence with IPD/PWD and find KmerData
    std::string seq("CCCCCG") ; // should insert CCCCC and CCCCG with 0 IPD/PWD
    ngsai::KmerData ccccc1("CCCCC", ipd_1, pwd_2) ;
    ngsai::KmerData ccccg1("CCCCG", ipd_1, pwd_2) ;
    ngsai::KmerData ccccc2("CCCCC", ipd_2, pwd_4) ;
    ngsai::KmerData ccccg2("CCCCG", ipd_2, pwd_4) ;
    ASSERT_NO_THROW(map1.insert(seq,
                                std::vector<uint32_t>(6, 1),
                                std::vector<uint32_t>(6, 2),
                                ngsai::KmerMap::insert_mode::REPLACE)) ;
    ASSERT_EQ(map1.size(), 4) ;
    auto iter5 = map1.find(ccccc1) ;
    auto iter6 = map1.find(ccccg1) ;
    ASSERT_EQ(*iter5 == std::make_pair(uint32_t(1), ccccc1), true) ;
    ASSERT_EQ(*iter6 == std::make_pair(uint32_t(1), ccccg1), true) ;

   
    // insert sum from sequence with IPD/PWD and find KmerData
    ASSERT_NO_THROW(map1.insert(seq,
                                std::vector<uint32_t>(6, 1),
                                std::vector<uint32_t>(6, 2),
                                ngsai::KmerMap::insert_mode::SUM)) ;
    ASSERT_EQ(map1.size(), 4) ;
    auto iter7 = map1.find(ccccc1) ;
    auto iter8 = map1.find(ccccg1) ;
    ASSERT_EQ(*iter7 == std::make_pair(uint32_t(2), ccccc2), true) ;
    ASSERT_EQ(*iter8 == std::make_pair(uint32_t(2), ccccg2), true) ;
    
    // find something not in map
    ngsai::KmerData aaaaa0("AAAAA", ipd_0, pwd_0) ;
    auto  iter9 = map1.find(aaaaa0) ;
    ASSERT_EQ(iter9 == map1.end(), true) ;
    
    // find from sequence
    ngsai::KmerMap map2(4) ;
    std::string acgtagccc   = "ACGTAGCCC" ;
    std::string acgtagcccaa = "ACGTAGCCCAA" ;
    std::vector<uint32_t> ipds = { 0, 1, 2, 3, 4, 5, 6, 7, 8} ;
    std::vector<uint32_t> pwds = {10,11,12,13,14,15,16,17,18} ;
    map2.insert(acgtagccc, ipds, pwds, ngsai::KmerMap::insert_mode::REPLACE) ;
    const ngsai::KmerMap map2c = map2 ;
    // find with exact same sequence as insertion was made
    auto   iters1 =  map2.find(acgtagccc) ;
    auto  citers1 = map2c.find(acgtagccc) ;
    // create vector of expected iterators
    std::vector<ngsai::kmerBucket> kmers_exp_acgtagccc = {
        std::make_pair(1, 
                       ngsai::KmerData("ACGT",
                                       std::vector<uint32_t>({ 0, 1, 2, 3}),
                                       std::vector<uint32_t>({10,11,12,13}))),
        std::make_pair(1, 
                       ngsai::KmerData("CGTA",
                                       std::vector<uint32_t>({ 1, 2, 3, 4}),
                                       std::vector<uint32_t>({11,12,13,14}))),
        std::make_pair(1,
                       ngsai::KmerData("GTAG",
                                       std::vector<uint32_t>({ 2, 3, 4, 5}),
                                       std::vector<uint32_t>({12,13,14,15}))),
        std::make_pair(1,
                       ngsai::KmerData("TAGC",
                                       std::vector<uint32_t>({ 3, 4, 5, 6}),
                                       std::vector<uint32_t>({13,14,15,16}))),
        std::make_pair(1,
                       ngsai::KmerData("AGCC",
                                       std::vector<uint32_t>({ 4, 5, 6, 7}),
                                       std::vector<uint32_t>({14,15,16,17}))),
        std::make_pair(1,
                       ngsai::KmerData("GCCC",
                                       std::vector<uint32_t>({ 5, 6, 7, 8}),
                                       std::vector<uint32_t>({15,16,17,18})))
    } ;
    ASSERT_EQ(iters1.size(), kmers_exp_acgtagccc.size()) ;
    ASSERT_EQ(*(iters1[0])  ==  kmers_exp_acgtagccc[0], true) ;
    ASSERT_EQ(*(iters1[1])  ==  kmers_exp_acgtagccc[1], true) ;
    ASSERT_EQ(*(iters1[2])  ==  kmers_exp_acgtagccc[2], true) ;
    ASSERT_EQ(*(iters1[3])  ==  kmers_exp_acgtagccc[3], true) ;
    ASSERT_EQ(*(iters1[4])  ==  kmers_exp_acgtagccc[4], true) ;
    ASSERT_EQ(*(iters1[5])  ==  kmers_exp_acgtagccc[5], true) ;
    ASSERT_EQ(citers1.size(), kmers_exp_acgtagccc.size()) ;
    ASSERT_EQ(*(citers1[0]) ==  kmers_exp_acgtagccc[0], true) ;
    ASSERT_EQ(*(citers1[1]) ==  kmers_exp_acgtagccc[1], true) ;
    ASSERT_EQ(*(citers1[2]) ==  kmers_exp_acgtagccc[2], true) ;
    ASSERT_EQ(*(citers1[3]) ==  kmers_exp_acgtagccc[3], true) ;
    ASSERT_EQ(*(citers1[4]) ==  kmers_exp_acgtagccc[4], true) ;
    ASSERT_EQ(*(citers1[5]) ==  kmers_exp_acgtagccc[5], true) ;
    
    // find with sequence extended compared to insertion
     iters1 =  map2.find(acgtagcccaa) ;
    citers1 = map2c.find(acgtagcccaa) ;
    ASSERT_EQ(iters1.size(), kmers_exp_acgtagccc.size() + 2) ;
    ASSERT_EQ(*(iters1[0]) ==  kmers_exp_acgtagccc[0], true) ;
    ASSERT_EQ(*(iters1[1]) ==  kmers_exp_acgtagccc[1], true) ;
    ASSERT_EQ(*(iters1[2]) ==  kmers_exp_acgtagccc[2], true) ;
    ASSERT_EQ(*(iters1[3]) ==  kmers_exp_acgtagccc[3], true) ;
    ASSERT_EQ(*(iters1[4]) ==  kmers_exp_acgtagccc[4], true) ;
    ASSERT_EQ(*(iters1[5]) ==  kmers_exp_acgtagccc[5], true) ;
    ASSERT_EQ(iters1[6] == map2.end(), true) ;
    ASSERT_EQ(iters1[7] == map2.end(), true) ;
    ASSERT_EQ(citers1.size(), kmers_exp_acgtagccc.size() + 2) ;
    ASSERT_EQ(*(citers1[0]) ==  kmers_exp_acgtagccc[0], true) ;
    ASSERT_EQ(*(citers1[1]) ==  kmers_exp_acgtagccc[1], true) ;
    ASSERT_EQ(*(citers1[2]) ==  kmers_exp_acgtagccc[2], true) ;
    ASSERT_EQ(*(citers1[3]) ==  kmers_exp_acgtagccc[3], true) ;
    ASSERT_EQ(*(citers1[4]) ==  kmers_exp_acgtagccc[4], true) ;
    ASSERT_EQ(*(citers1[5]) ==  kmers_exp_acgtagccc[5], true) ;
    ASSERT_EQ(citers1[6] == map2c.end(), true) ;
    ASSERT_EQ(citers1[7] == map2c.end(), true) ;
}


/*!
 *\brief Tests that the at() method works properly.
 */
TEST(KmerMapTest, at_method)
{
    // empty map
    size_t kmer_size = 5 ;
    KmerMapTest map(kmer_size) ;
    
    ngsai::KmerData kmer1("ACGTA",
                          std::vector<uint32_t>({1,1,1,1,1}),
                          std::vector<uint32_t>({2,2,2,2,2})) ;
    ngsai::KmerData kmer2("CCCCC",
                          std::vector<uint32_t>({3,3,3,3,3}),
                          std::vector<uint32_t>({4,4,4,4,4})) ;
    ngsai::KmerData kmer3("GGGGG",
                          std::vector<uint32_t>({5,5,5,5,5}),
                          std::vector<uint32_t>({6,6,6,6,6})) ;
    ngsai::KmerData kmer4("CGCGC",
                          std::vector<uint32_t>({5,5,5,5,5}),
                          std::vector<uint32_t>({6,6,6,6,6})) ;
    
    // nothing in map
    EXPECT_THROW(auto  ref = map.at(kmer1), std::out_of_range) ;
    EXPECT_THROW(auto  ref = map.at(kmer2), std::out_of_range) ;
    EXPECT_THROW(auto  ref = map.at(kmer3), std::out_of_range) ;
    EXPECT_THROW(auto  ref = map.at(kmer4), std::out_of_range) ;

    // insert
    map.insert(kmer1) ;
    map.insert(kmer2) ;
    map.insert(kmer3) ;

    // const map
    const KmerMapTest cmap = map ;

    // retrieve in map
    auto cref1 = cmap.at(kmer1) ;
    auto cref2 = cmap.at(kmer2) ;
    auto cref3 = cmap.at(kmer3) ;
    auto  ref1 =  map.at(kmer1) ;
    auto  ref2 =  map.at(kmer2) ;
    auto  ref3 =  map.at(kmer3) ;
    ASSERT_EQ(cref1 == std::make_pair(uint32_t(1), kmer1), true) ;
    ASSERT_EQ(cref2 == std::make_pair(uint32_t(1), kmer2), true) ;
    ASSERT_EQ(cref3 == std::make_pair(uint32_t(1), kmer3), true) ;
    ASSERT_EQ(ref1  == std::make_pair(uint32_t(1), kmer1), true) ;
    ASSERT_EQ(ref2  == std::make_pair(uint32_t(1), kmer2), true) ;
    ASSERT_EQ(ref3  == std::make_pair(uint32_t(1), kmer3), true) ;

    // still not in map
    EXPECT_THROW(auto   refx =  map.at(kmer4), std::out_of_range) ;
    EXPECT_THROW(auto  crefx = cmap.at(kmer4), std::out_of_range) ;
}
