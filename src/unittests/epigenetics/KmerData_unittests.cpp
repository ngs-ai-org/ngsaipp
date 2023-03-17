#include <gtest/gtest.h>
#include <gmock/gmock.h>
#include <utility>         // std::move


#include <ngsaipp/epigenetics/KmerData.hpp>


using testing::ElementsAre ;
using testing::ContainerEq ;


/*
 * This file contains tests for the ngsai::KmerData class from 
 * src/epigenetics/KmerData.cpp
 */


/*!
 * \brief Tests that the constructor works properly and throw the expected
 * std::invalid_argument exceptions if given invalid parameters.
 */
TEST(KmerDataTest, constructor)
{   // default constructor
    ngsai::KmerData kmer1 ;
    ASSERT_EQ(kmer1.sequence, std::string()) ;
    ASSERT_EQ(kmer1.ipd, std::vector<uint32_t>()) ;
    ASSERT_EQ(kmer1.pwd, std::vector<uint32_t>()) ;

    // arg constructor
    std::string sequence = "ACGTA" ;
    std::vector<uint32_t> ipd({0,1,2,3,4}) ;
    std::vector<uint32_t> pwd({5,6,7,8,9}) ;
    ngsai::KmerData kmer2(sequence, ipd, pwd) ;
    EXPECT_THAT(kmer2.sequence, ContainerEq(sequence));
    EXPECT_THAT(kmer2.ipd, ContainerEq(ipd));
    EXPECT_THAT(kmer2.pwd, ContainerEq(pwd));
    EXPECT_THROW(ngsai::KmerData("ACGT",
                                 std::vector<uint32_t>({0,1,2,3,4}),
                                 std::vector<uint32_t>({0,1,2,3,4})),
                 std::invalid_argument) ;
    EXPECT_THROW(ngsai::KmerData("ACGTA",
                                 std::vector<uint32_t>({0,1,2,3}),
                                 std::vector<uint32_t>({0,1,2,3,4})),
                 std::invalid_argument) ;
    EXPECT_THROW(ngsai::KmerData("ACGTA",
                                 std::vector<uint32_t>({0,1,2,3,4}),
                                 std::vector<uint32_t>({0,1,2,3})),
                 std::invalid_argument) ;

    // copy constructor
    ngsai::KmerData kmer3(kmer2) ;
    EXPECT_THAT(kmer3.sequence, ContainerEq(sequence));
    EXPECT_THAT(kmer3.ipd,      ContainerEq(ipd));
    EXPECT_THAT(kmer3.pwd,      ContainerEq(pwd));
    kmer2.sequence = "ACGT" ;
    EXPECT_THROW(ngsai::KmerData{kmer2}, std::invalid_argument) ;
    kmer2.sequence = "ACGTA" ;
    kmer2.ipd.pop_back() ;
    EXPECT_THROW(ngsai::KmerData{kmer2}, std::invalid_argument) ;
    kmer2.ipd.push_back(4) ;
    kmer2.pwd.pop_back() ;
    EXPECT_THROW(ngsai::KmerData{kmer2}, std::invalid_argument) ;
    kmer2.pwd.push_back(4) ;

    // move constructor
    ngsai::KmerData kmer4(ngsai::KmerData(sequence, ipd, pwd)) ;
    EXPECT_THAT(kmer4.sequence, ContainerEq(sequence));
    EXPECT_THAT(kmer4.ipd,      ContainerEq(ipd));
    EXPECT_THAT(kmer4.pwd,      ContainerEq(pwd));
    ngsai::KmerData kmer_tmp(sequence, ipd, pwd) ;
    kmer_tmp.sequence = "ACGT" ;
    EXPECT_THROW(ngsai::KmerData{std::move(kmer_tmp)}, std::invalid_argument) ;
    kmer_tmp = ngsai::KmerData(sequence, ipd, pwd) ;
    kmer_tmp.ipd.pop_back() ;
    EXPECT_THROW(ngsai::KmerData{std::move(kmer_tmp)}, std::invalid_argument) ;
    kmer_tmp = ngsai::KmerData(sequence, ipd, pwd) ;
    kmer_tmp.pwd.pop_back() ;
    EXPECT_THROW(ngsai::KmerData{std::move(kmer_tmp)}, std::invalid_argument) ;
}


/*!
 * \brief Tests that the assignment operator works properly and throw the 
 * expected std::invalid_argument exception if given invalid parameters.
 */
TEST(KmerDataTest, assignment_operator)
{
    std::string sequence = "ACGTA" ;
    std::vector<uint32_t> ipd({0,1,2,3,4}) ;
    std::vector<uint32_t> pwd({5,6,7,8,9}) ;
    ngsai::KmerData kmer1(sequence, ipd, pwd) ;

    ngsai::KmerData kmer2 = kmer1 ;

    EXPECT_THAT(kmer2.sequence, ContainerEq(sequence));
    EXPECT_THAT(kmer2.ipd,      ContainerEq(ipd));
    EXPECT_THAT(kmer2.pwd,      ContainerEq(pwd));
    kmer2.sequence = "ACGT" ;
    EXPECT_THROW(kmer1 = kmer2, std::invalid_argument) ;
    kmer2.sequence = "ACGTA" ;
    kmer2.ipd.pop_back() ;
    EXPECT_THROW(kmer1 = kmer2, std::invalid_argument) ;
    kmer2.ipd.push_back(4) ;
    kmer2.pwd.pop_back() ;
    EXPECT_THROW(kmer1 = kmer2, std::invalid_argument) ;
    kmer2.pwd.push_back(4) ;
}


/*!
 * \brief Tests that the move assignment operatorr works properly and throw the 
 * expected std::invalid_argument exception if given invalid parameters.
 */
TEST(KmerDataTest, move_assignment_operator)
{
    std::string sequence = "ACGTA" ;
    std::vector<uint32_t> ipd({0,1,2,3,4}) ;
    std::vector<uint32_t> pwd({5,6,7,8,9}) ;
    ngsai::KmerData kmer_tmp(sequence, ipd, pwd) ;

    ngsai::KmerData kmer1 = std::move(kmer_tmp) ;
    EXPECT_THAT(kmer1.sequence, ContainerEq(sequence));
    EXPECT_THAT(kmer1.ipd,      ContainerEq(ipd));
    EXPECT_THAT(kmer1.pwd,      ContainerEq(pwd));

    // to do the throw tests
    kmer_tmp = ngsai::KmerData(sequence, ipd, pwd) ;
    kmer_tmp.sequence = "ACGT" ;
    EXPECT_THROW(kmer1 = std::move(kmer_tmp), std::invalid_argument) ;

    kmer_tmp = ngsai::KmerData(sequence, ipd, pwd) ;
    kmer_tmp.ipd.pop_back() ;
    EXPECT_THROW(kmer1 = std::move(kmer_tmp), std::invalid_argument) ;

    kmer_tmp = ngsai::KmerData(sequence, ipd, pwd) ;
    kmer_tmp.pwd.pop_back() ;
    EXPECT_THROW(kmer1 = std::move(kmer_tmp), std::invalid_argument) ;
}


/*!
 * \brief Tests that the equality operatorr works properly.
 */
TEST(KmerDataTest, equality_operator)
{

    ngsai::KmerData kmer1("ACGTA",
                          std::vector<uint32_t>({1,1,1,1,1}),
                          std::vector<uint32_t>({2,2,2,2,2})) ;
    ngsai::KmerData kmer2("ACGTA",
                          std::vector<uint32_t>({1,1,1,1,1}),
                          std::vector<uint32_t>({2,2,2,2,2})) ;
    ASSERT_EQ(kmer1 == kmer1, true) 
        << "failed to assess that kmers are identical";
    ASSERT_EQ(kmer2 == kmer2, true) 
        << "failed to assess that kmers are identical";
    ASSERT_EQ(kmer1 == kmer2, true) 
        << "failed to assess that kmers are identical";
    
    
    ngsai::KmerData kmer3("GGGGG",
                          std::vector<uint32_t>({1,1,1,1,1}),
                          std::vector<uint32_t>({2,2,2,2,2})) ;
    ASSERT_EQ(kmer1 == kmer3, false) 
        << "failed to assess that kmers are different";
    
    ngsai::KmerData kmer4("ACGTA",
                          std::vector<uint32_t>({0,0,0,0,0}),
                          std::vector<uint32_t>({2,2,2,2,2})) ;
    ASSERT_EQ(kmer1 == kmer4, false) 
        << "failed to assess that kmers are different";

    ngsai::KmerData kmer5("ACGTA",
                          std::vector<uint32_t>({1,1,1,1,1}),
                          std::vector<uint32_t>({0,0,0,0,0})) ;
    ASSERT_EQ(kmer1 == kmer4, false) 
        << "failed to assess that kmers are different";
}


/*!
 * \brief Tests that the inequality operatorr works properly.
 */
TEST(KmerDataTest, inequality_operator)
{
    
    ngsai::KmerData kmer1("ACGTA",
                          std::vector<uint32_t>({1,1,1,1,1}),
                          std::vector<uint32_t>({2,2,2,2,2})) ;
    ngsai::KmerData kmer2("ACGTA",
                          std::vector<uint32_t>({1,1,1,1,1}),
                          std::vector<uint32_t>({2,2,2,2,2})) ;
    ASSERT_EQ(kmer1 != kmer1, false) 
        << "failed to assess that kmers are identical";
    ASSERT_EQ(kmer2 != kmer2, false) 
        << "failed to assess that kmers are identical";
    ASSERT_EQ(kmer1 != kmer2, false) 
        << "failed to assess that kmers are identical";
    
    
    ngsai::KmerData kmer3("GGGGG",
                          std::vector<uint32_t>({1,1,1,1,1}),
                          std::vector<uint32_t>({2,2,2,2,2})) ;
    ASSERT_EQ(kmer1 != kmer3, true) 
        << "failed to assess that kmers are different";
    
    ngsai::KmerData kmer4("ACGTA",
                          std::vector<uint32_t>({0,0,0,0,0}),
                          std::vector<uint32_t>({2,2,2,2,2})) ;
    ASSERT_EQ(kmer1 != kmer4, true) 
        << "failed to assess that kmers are different";

    ngsai::KmerData kmer5("ACGTA",
                          std::vector<uint32_t>({1,1,1,1,1}),
                          std::vector<uint32_t>({0,0,0,0,0})) ;
    ASSERT_EQ(kmer1 != kmer4, true) 
        << "failed to assess that kmers are different";
}


/*!
 * \brief Tests that the size() method works properly.
 */
TEST(KmerDataTest, size)
{
    std::string sequence = "ACGTA" ;
    std::vector<uint32_t> ipd({0,1,2,3,4}) ;
    std::vector<uint32_t> pwd({5,6,7,8,9}) ;
    ngsai::KmerData kmer1(sequence, ipd, pwd) ;
    ASSERT_EQ(kmer1.size(), sequence.size());
}