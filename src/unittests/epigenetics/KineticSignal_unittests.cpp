#include <gtest/gtest.h>
#include <gmock/gmock.h>
#include <string>
#include <vector>
#include <stdexcept>

#include <ngsaipp/epigenetics/KineticSignal.hpp>


using testing::ElementsAre ;
using testing::ContainerEq ;


/*
 * This file ontains tests for the ngsai::KineticSignal class from 
 * src/epigenetics/KineticSignal.cpp
 */


/*!
 * \brief Tests that the default constructor works properly.
 */
TEST(KineticSignalTest, constructor_default)
{
    ngsai::KineticSignal kinetic ;

    ASSERT_EQ(kinetic.size(), 0) ;
    ASSERT_THAT(kinetic.getSequenceFw(), ContainerEq(std::string(""))) ;
    ASSERT_THAT(kinetic.getSequenceRv(), ContainerEq(std::string(""))) ;
    ASSERT_THAT(kinetic.getIPDFw(),      ContainerEq(std::vector<double>())) ;
    ASSERT_THAT(kinetic.getIPDRv(),      ContainerEq(std::vector<double>())) ;
    ASSERT_THAT(kinetic.getPWDFw(),      ContainerEq(std::vector<double>())) ;
    ASSERT_THAT(kinetic.getPWDRv(),      ContainerEq(std::vector<double>())) ;   
}


/*!
 * \brief Tests that the argument constructor works properly.
 */
TEST(KineticSignalTest, constructor_argument)
{
    std::string seq_0 ;
    std::string seq_fw("ACGT") ;
    std::string seq_rv("CGCG") ;
    std::string seq_err1 ;
    std::string seq_err2("ACGTACGT") ;
    std::vector<double> ipd_fw({1.,2.,3.,4.}) ;
    std::vector<double> ipd_rv({5.,6.,7.,8.}) ;
    std::vector<double> ipd_err1 ;
    std::vector<double> ipd_err2({1.,2.,3.,4.,5.,6.,7.,8.}) ;
    std::vector<double> pwd_fw({1.,2.,3.,4.}) ;
    std::vector<double> pwd_rv({5.,6.,7.,8.}) ;
    std::vector<double> pwd_err1 ;
    std::vector<double> pwd_err2({1.,2.,3.,4.,5.,6.,7.,8.}) ;

    // must work
    ngsai::KineticSignal kinetic(seq_fw, seq_rv, ipd_fw, ipd_rv, pwd_fw, pwd_rv) ;
    ASSERT_EQ(kinetic.size(), seq_fw.size()) ;
    ASSERT_THAT(kinetic.getSequenceFw(), ContainerEq(seq_fw)) ;
    ASSERT_THAT(kinetic.getSequenceRv(), ContainerEq(seq_rv)) ;
    ASSERT_THAT(kinetic.getIPDFw(),      ContainerEq(ipd_fw)) ;
    ASSERT_THAT(kinetic.getIPDRv(),      ContainerEq(ipd_rv)) ;
    ASSERT_THAT(kinetic.getPWDFw(),      ContainerEq(pwd_fw)) ;
    ASSERT_THAT(kinetic.getPWDRv(),      ContainerEq(pwd_rv)) ;

    // seq fw has diff length -> must fail
    EXPECT_THROW(ngsai::KineticSignal kinetic(seq_err1,
                                           seq_rv,
                                           ipd_fw,
                                           ipd_rv,
                                           pwd_fw,
                                           pwd_rv),
                 std::invalid_argument) ;
    EXPECT_THROW(ngsai::KineticSignal kinetic(seq_err2,
                                           seq_rv,
                                           ipd_fw,
                                           ipd_rv,
                                           pwd_fw,
                                           pwd_rv),
                 std::invalid_argument) ;
    
    // seq rv has diff length -> must fail
    EXPECT_THROW(ngsai::KineticSignal kinetic(seq_fw,
                                           seq_err1,
                                           ipd_fw,
                                           ipd_rv,
                                           pwd_fw,
                                           pwd_rv),
                 std::invalid_argument) ;
    EXPECT_THROW(ngsai::KineticSignal kinetic(seq_fw,
                                           seq_err2,
                                           ipd_fw,
                                           ipd_rv,
                                           pwd_fw,
                                           pwd_rv),
                 std::invalid_argument) ;
    
    // ipd fw has diff length -> must fail
    EXPECT_THROW(ngsai::KineticSignal kinetic(seq_fw,
                                           seq_rv,
                                           ipd_err1,
                                           ipd_rv,
                                           pwd_fw,
                                           pwd_rv),
                 std::invalid_argument) ;
    EXPECT_THROW(ngsai::KineticSignal kinetic(seq_fw,
                                           seq_rv,
                                           ipd_err2,
                                           ipd_rv,
                                           pwd_fw,
                                           pwd_rv),
                 std::invalid_argument) ;
    
    // ipd rv has diff length -> must fail
    EXPECT_THROW(ngsai::KineticSignal kinetic(seq_fw,
                                           seq_rv,
                                           ipd_fw,
                                           ipd_err1,
                                           pwd_fw,
                                           pwd_rv),
                 std::invalid_argument) ;
    EXPECT_THROW(ngsai::KineticSignal kinetic(seq_fw,
                                           seq_rv,
                                           ipd_fw,
                                           ipd_err2,
                                           pwd_fw,
                                           pwd_rv),
                 std::invalid_argument) ;
    
    // pwd fw has diff length -> must fail
    EXPECT_THROW(ngsai::KineticSignal kinetic(seq_fw,
                                           seq_rv,
                                           ipd_fw,
                                           ipd_rv,
                                           pwd_err1,
                                           pwd_rv),
                 std::invalid_argument) ;
    EXPECT_THROW(ngsai::KineticSignal kinetic(seq_fw,
                                           seq_rv,
                                           ipd_fw,
                                           ipd_rv,
                                           pwd_err2,
                                           pwd_rv),
                 std::invalid_argument) ;
    
    // pwd rv has diff length -> must fail
    EXPECT_THROW(ngsai::KineticSignal kinetic(seq_fw,
                                           seq_rv,
                                           ipd_fw,
                                           ipd_rv,
                                           pwd_fw,
                                           pwd_err1),
                 std::invalid_argument) ;
    EXPECT_THROW(ngsai::KineticSignal kinetic(seq_fw,
                                           seq_rv,
                                           ipd_fw,
                                           ipd_rv,
                                           pwd_fw,
                                           pwd_err2),
                 std::invalid_argument) ;
}


/*!
 * \brief Tests that the copy constructor works properly.
 */
TEST(KineticSignalTest, constructor_copy)
{
    std::string seq_0 ;
    std::string seq_fw("ACGT") ;
    std::string seq_rv("CGCG") ;
    std::string seq_err1 ;
    std::string seq_err2("ACGTACGT") ;
    std::vector<double> ipd_fw({1.,2.,3.,4.}) ;
    std::vector<double> ipd_rv({5.,6.,7.,8.}) ;
    std::vector<double> ipd_err1 ;
    std::vector<double> ipd_err2({1.,2.,3.,4.,5.,6.,7.,8.}) ;
    std::vector<double> pwd_fw({1.,2.,3.,4.}) ;
    std::vector<double> pwd_rv({5.,6.,7.,8.}) ;
    std::vector<double> pwd_err1 ;
    std::vector<double> pwd_err2({1.,2.,3.,4.,5.,6.,7.,8.}) ;

    // must work
    ngsai::KineticSignal kinetic_1(seq_fw, seq_rv, 
                                ipd_fw, ipd_rv,
                                pwd_fw, pwd_rv) ;
    ngsai::KineticSignal kinetic_2(kinetic_1) ;
    ASSERT_EQ(kinetic_2.size(), kinetic_1.size()) ;
    ASSERT_THAT(kinetic_2.getSequenceFw(),
                ContainerEq(kinetic_1.getSequenceFw())) ;
    ASSERT_THAT(kinetic_2.getSequenceRv(),
                ContainerEq(kinetic_1.getSequenceRv())) ;
    ASSERT_THAT(kinetic_2.getIPDFw(),
                ContainerEq(kinetic_1.getIPDFw())) ;
    ASSERT_THAT(kinetic_2.getIPDRv(),
                ContainerEq(kinetic_1.getIPDRv())) ;
    ASSERT_THAT(kinetic_2.getPWDFw(),
                ContainerEq(kinetic_1.getPWDFw())) ;
    ASSERT_THAT(kinetic_2.getPWDRv(),
                ContainerEq(kinetic_1.getPWDRv())) ;
}


/*!
 * \brief Tests that the move constructor works properly.
 */
TEST(KineticSignalTest, constructor_move)
{
    std::string seq_0 ;
    std::string seq_fw("ACGT") ;
    std::string seq_rv("CGCG") ;
    std::string seq_err1 ;
    std::string seq_err2("ACGTACGT") ;
    std::vector<double> ipd_fw({1.,2.,3.,4.}) ;
    std::vector<double> ipd_rv({5.,6.,7.,8.}) ;
    std::vector<double> ipd_err1 ;
    std::vector<double> ipd_err2({1.,2.,3.,4.,5.,6.,7.,8.}) ;
    std::vector<double> pwd_fw({1.,2.,3.,4.}) ;
    std::vector<double> pwd_rv({5.,6.,7.,8.}) ;
    std::vector<double> pwd_err1 ;
    std::vector<double> pwd_err2({1.,2.,3.,4.,5.,6.,7.,8.}) ;

    // must work
    ngsai::KineticSignal kinetic_1(seq_fw, seq_rv, 
                                ipd_fw, ipd_rv,
                                pwd_fw, pwd_rv) ;
    ngsai::KineticSignal kinetic_tmp(kinetic_1) ;
    ngsai::KineticSignal kinetic_2(std::move(kinetic_tmp)) ;
    ASSERT_EQ(kinetic_2.size(), kinetic_1.size()) ;
    ASSERT_THAT(kinetic_2.getSequenceFw(),
                ContainerEq(kinetic_1.getSequenceFw())) ;
    ASSERT_THAT(kinetic_2.getSequenceRv(),
                ContainerEq(kinetic_1.getSequenceRv())) ;
    ASSERT_THAT(kinetic_2.getIPDFw(),
                ContainerEq(kinetic_1.getIPDFw())) ;
    ASSERT_THAT(kinetic_2.getIPDRv(),
                ContainerEq(kinetic_1.getIPDRv())) ;
    ASSERT_THAT(kinetic_2.getPWDFw(),
                ContainerEq(kinetic_1.getPWDFw())) ;
    ASSERT_THAT(kinetic_2.getPWDRv(),
                ContainerEq(kinetic_1.getPWDRv())) ;
}


/*!
 * \brief Constructs an empty instance and populate the sequence fw first, 
 * then the other fields and checks that it behaves properly.
 */
TEST(KineticSignalTest, sequence_fw)
{   
    std::string seq_0 ;
    std::string seq_1("ACGT") ;
    std::string seq_2("CGCG") ;
    std::string seq_3("ACGTACGTACGT") ;

    std::vector<double> ipd_0 ;
    std::vector<double> ipd_1({1.,2.,3.,4.}) ;
    std::vector<double> ipd_2({5.,6.,7.,8.}) ;
    std::vector<double> ipd_3({1.,2.,3.,4.,5.,6.,7.,8.}) ;

    std::vector<double> pwd_0 ;
    std::vector<double> pwd_1({1.,2.,3.,4.}) ;
    std::vector<double> pwd_2({5.,6.,7.,8.}) ;
    std::vector<double> pwd_3({1.,2.,3.,4.,5.,6.,7.,8.}) ;

    // init using setSequenceFw by reference
    ngsai::KineticSignal kinetic_1 ;
    kinetic_1.setSequenceFw(seq_1) ; // sets size from sequence fw
    kinetic_1.setSequenceRv(seq_2) ;
    kinetic_1.setIPDFw(ipd_1) ;
    kinetic_1.setIPDRv(ipd_2) ;
    kinetic_1.setPWDFw(pwd_1) ;
    kinetic_1.setPWDRv(pwd_2) ;
    ASSERT_EQ(kinetic_1.size(), seq_1.size()) ;
    ASSERT_THAT(kinetic_1.getSequenceFw(), ContainerEq(seq_1)) ;
    ASSERT_THAT(kinetic_1.getSequenceRv(), ContainerEq(seq_2)) ;
    ASSERT_THAT(kinetic_1.getIPDFw(), ContainerEq(ipd_1)) ;
    ASSERT_THAT(kinetic_1.getIPDRv(), ContainerEq(ipd_2)) ;
    ASSERT_THAT(kinetic_1.getPWDFw(), ContainerEq(pwd_1)) ;
    ASSERT_THAT(kinetic_1.getPWDRv(), ContainerEq(pwd_2)) ;
    // modify field values with same length values -> must work
    kinetic_1.setSequenceFw(seq_2) ;
    kinetic_1.setSequenceRv(seq_1) ;
    kinetic_1.setIPDFw(ipd_2) ;
    kinetic_1.setIPDRv(ipd_1) ;
    kinetic_1.setPWDFw(pwd_2) ;
    kinetic_1.setPWDRv(pwd_1) ;
    ASSERT_EQ(kinetic_1.size(), seq_1.size()) ;
    ASSERT_THAT(kinetic_1.getSequenceFw(), ContainerEq(seq_2)) ;
    ASSERT_THAT(kinetic_1.getSequenceRv(), ContainerEq(seq_1)) ;
    ASSERT_THAT(kinetic_1.getIPDFw(), ContainerEq(ipd_2)) ;
    ASSERT_THAT(kinetic_1.getIPDRv(), ContainerEq(ipd_1)) ;
    ASSERT_THAT(kinetic_1.getPWDFw(), ContainerEq(pwd_2)) ;
    ASSERT_THAT(kinetic_1.getPWDRv(), ContainerEq(pwd_1)) ;     
    // modify field values with diff length values ->  must fail
    EXPECT_THROW(kinetic_1.setSequenceFw(seq_0), std::invalid_argument) ;
    EXPECT_THROW(kinetic_1.setSequenceFw(seq_3), std::invalid_argument) ;
    EXPECT_THROW(kinetic_1.setSequenceRv(seq_0), std::invalid_argument) ;
    EXPECT_THROW(kinetic_1.setSequenceRv(seq_3), std::invalid_argument) ;
    EXPECT_THROW(kinetic_1.setIPDFw(ipd_0),      std::invalid_argument) ;
    EXPECT_THROW(kinetic_1.setIPDFw(ipd_3),      std::invalid_argument) ;
    EXPECT_THROW(kinetic_1.setIPDRv(ipd_0),      std::invalid_argument) ;
    EXPECT_THROW(kinetic_1.setIPDRv(ipd_3),      std::invalid_argument) ;
    EXPECT_THROW(kinetic_1.setPWDFw(pwd_0),      std::invalid_argument) ;
    EXPECT_THROW(kinetic_1.setPWDFw(pwd_3),      std::invalid_argument) ;
    EXPECT_THROW(kinetic_1.setPWDRv(pwd_0),      std::invalid_argument) ;
    EXPECT_THROW(kinetic_1.setPWDRv(pwd_3),      std::invalid_argument) ;


    // init using setSequenceFw by rvalue
    ngsai::KineticSignal kinetic_2 ;
    kinetic_2.setSequenceFw(std::string(seq_1)) ; // sets size from sequence fw
    kinetic_2.setSequenceRv(seq_2) ;
    kinetic_2.setIPDFw(ipd_1) ;
    kinetic_2.setIPDRv(ipd_2) ;
    kinetic_2.setPWDFw(pwd_1) ;
    kinetic_2.setPWDRv(pwd_2) ;
    ASSERT_EQ(kinetic_2.size(), seq_1.size()) ;
    ASSERT_THAT(kinetic_2.getSequenceFw(), ContainerEq(seq_1)) ;
    ASSERT_THAT(kinetic_2.getSequenceRv(), ContainerEq(seq_2)) ;
    ASSERT_THAT(kinetic_2.getIPDFw(), ContainerEq(ipd_1)) ;
    ASSERT_THAT(kinetic_2.getIPDRv(), ContainerEq(ipd_2)) ;
    ASSERT_THAT(kinetic_2.getPWDFw(), ContainerEq(pwd_1)) ;
    ASSERT_THAT(kinetic_2.getPWDRv(), ContainerEq(pwd_2)) ;
    // modify field values with same length values -> must work
    kinetic_2.setSequenceFw(seq_2) ;
    kinetic_2.setSequenceRv(seq_1) ;
    kinetic_2.setIPDFw(ipd_2) ;
    kinetic_2.setIPDRv(ipd_1) ;
    kinetic_2.setPWDFw(pwd_2) ;
    kinetic_2.setPWDRv(pwd_1) ;
    ASSERT_EQ(kinetic_2.size(), seq_1.size()) ;
    ASSERT_THAT(kinetic_2.getSequenceFw(), ContainerEq(seq_2)) ;
    ASSERT_THAT(kinetic_2.getSequenceRv(), ContainerEq(seq_1)) ;
    ASSERT_THAT(kinetic_2.getIPDFw(), ContainerEq(ipd_2)) ;
    ASSERT_THAT(kinetic_2.getIPDRv(), ContainerEq(ipd_1)) ;
    ASSERT_THAT(kinetic_2.getPWDFw(), ContainerEq(pwd_2)) ;
    ASSERT_THAT(kinetic_2.getPWDRv(), ContainerEq(pwd_1)) ;     
    // modify field values with diff length values ->  must fail
    EXPECT_THROW(kinetic_2.setSequenceFw(seq_0), std::invalid_argument) ;
    EXPECT_THROW(kinetic_2.setSequenceFw(seq_3), std::invalid_argument) ;
    EXPECT_THROW(kinetic_2.setSequenceRv(seq_0), std::invalid_argument) ;
    EXPECT_THROW(kinetic_2.setSequenceRv(seq_3), std::invalid_argument) ;
    EXPECT_THROW(kinetic_2.setIPDFw(ipd_0),      std::invalid_argument) ;
    EXPECT_THROW(kinetic_2.setIPDFw(ipd_3),      std::invalid_argument) ;
    EXPECT_THROW(kinetic_2.setIPDRv(ipd_0),      std::invalid_argument) ;
    EXPECT_THROW(kinetic_2.setIPDRv(ipd_3),      std::invalid_argument) ;
    EXPECT_THROW(kinetic_2.setPWDFw(pwd_0),      std::invalid_argument) ;
    EXPECT_THROW(kinetic_2.setPWDFw(pwd_3),      std::invalid_argument) ;
    EXPECT_THROW(kinetic_2.setPWDRv(pwd_0),      std::invalid_argument) ;
    EXPECT_THROW(kinetic_2.setPWDRv(pwd_3),      std::invalid_argument) ;
}


/*!
 * \brief Constructs an empty instance and populate the sequence rv first, 
 * then the other fields and checks that it behaves properly.
 */
TEST(KineticSignalTest, sequence_rv)
{   
    std::string seq_0 ;
    std::string seq_1("ACGT") ;
    std::string seq_2("CGCG") ;
    std::string seq_3("ACGTACGTACGT") ;

    std::vector<double> ipd_0 ;
    std::vector<double> ipd_1({1.,2.,3.,4.}) ;
    std::vector<double> ipd_2({5.,6.,7.,8.}) ;
    std::vector<double> ipd_3({1.,2.,3.,4.,5.,6.,7.,8.}) ;

    std::vector<double> pwd_0 ;
    std::vector<double> pwd_1({1.,2.,3.,4.}) ;
    std::vector<double> pwd_2({5.,6.,7.,8.}) ;
    std::vector<double> pwd_3({1.,2.,3.,4.,5.,6.,7.,8.}) ;

    // init using setSequenceRv by reference
    ngsai::KineticSignal kinetic_1 ;
    kinetic_1.setSequenceRv(seq_2) ; // sets size from sequence rv
    kinetic_1.setSequenceFw(seq_1) ;
    kinetic_1.setIPDFw(ipd_1) ;
    kinetic_1.setIPDRv(ipd_2) ;
    kinetic_1.setPWDFw(pwd_1) ;
    kinetic_1.setPWDRv(pwd_2) ;
    ASSERT_EQ(kinetic_1.size(), seq_1.size()) ;
    ASSERT_THAT(kinetic_1.getSequenceFw(), ContainerEq(seq_1)) ;
    ASSERT_THAT(kinetic_1.getSequenceRv(), ContainerEq(seq_2)) ;
    ASSERT_THAT(kinetic_1.getIPDFw(), ContainerEq(ipd_1)) ;
    ASSERT_THAT(kinetic_1.getIPDRv(), ContainerEq(ipd_2)) ;
    ASSERT_THAT(kinetic_1.getPWDFw(), ContainerEq(pwd_1)) ;
    ASSERT_THAT(kinetic_1.getPWDRv(), ContainerEq(pwd_2)) ;
    // modify field values with same length values -> must work
    kinetic_1.setSequenceFw(seq_2) ;
    kinetic_1.setSequenceRv(seq_1) ;
    kinetic_1.setIPDFw(ipd_2) ;
    kinetic_1.setIPDRv(ipd_1) ;
    kinetic_1.setPWDFw(pwd_2) ;
    kinetic_1.setPWDRv(pwd_1) ;
    ASSERT_EQ(kinetic_1.size(), seq_1.size()) ;
    ASSERT_THAT(kinetic_1.getSequenceFw(), ContainerEq(seq_2)) ;
    ASSERT_THAT(kinetic_1.getSequenceRv(), ContainerEq(seq_1)) ;
    ASSERT_THAT(kinetic_1.getIPDFw(), ContainerEq(ipd_2)) ;
    ASSERT_THAT(kinetic_1.getIPDRv(), ContainerEq(ipd_1)) ;
    ASSERT_THAT(kinetic_1.getPWDFw(), ContainerEq(pwd_2)) ;
    ASSERT_THAT(kinetic_1.getPWDRv(), ContainerEq(pwd_1)) ;     
    // modify field values with diff length values ->  must fail
    EXPECT_THROW(kinetic_1.setSequenceFw(seq_0), std::invalid_argument) ;
    EXPECT_THROW(kinetic_1.setSequenceFw(seq_3), std::invalid_argument) ;
    EXPECT_THROW(kinetic_1.setSequenceRv(seq_0), std::invalid_argument) ;
    EXPECT_THROW(kinetic_1.setSequenceRv(seq_3), std::invalid_argument) ;
    EXPECT_THROW(kinetic_1.setIPDFw(ipd_0),      std::invalid_argument) ;
    EXPECT_THROW(kinetic_1.setIPDFw(ipd_3),      std::invalid_argument) ;
    EXPECT_THROW(kinetic_1.setIPDRv(ipd_0),      std::invalid_argument) ;
    EXPECT_THROW(kinetic_1.setIPDRv(ipd_3),      std::invalid_argument) ;
    EXPECT_THROW(kinetic_1.setPWDFw(pwd_0),      std::invalid_argument) ;
    EXPECT_THROW(kinetic_1.setPWDFw(pwd_3),      std::invalid_argument) ;
    EXPECT_THROW(kinetic_1.setPWDRv(pwd_0),      std::invalid_argument) ;
    EXPECT_THROW(kinetic_1.setPWDRv(pwd_3),      std::invalid_argument) ;


    // init using setSequenceRv by rvalue
    ngsai::KineticSignal kinetic_2 ;
    kinetic_2.setSequenceRv(std::string(seq_2)) ; // sets size from sequence fw
    kinetic_2.setSequenceFw(seq_1) ;
    kinetic_2.setIPDFw(ipd_1) ;
    kinetic_2.setIPDRv(ipd_2) ;
    kinetic_2.setPWDFw(pwd_1) ;
    kinetic_2.setPWDRv(pwd_2) ;
    ASSERT_EQ(kinetic_2.size(), seq_1.size()) ;
    ASSERT_THAT(kinetic_2.getSequenceFw(), ContainerEq(seq_1)) ;
    ASSERT_THAT(kinetic_2.getSequenceRv(), ContainerEq(seq_2)) ;
    ASSERT_THAT(kinetic_2.getIPDFw(), ContainerEq(ipd_1)) ;
    ASSERT_THAT(kinetic_2.getIPDRv(), ContainerEq(ipd_2)) ;
    ASSERT_THAT(kinetic_2.getPWDFw(), ContainerEq(pwd_1)) ;
    ASSERT_THAT(kinetic_2.getPWDRv(), ContainerEq(pwd_2)) ;
    // modify field values with same length values -> must work
    kinetic_2.setSequenceFw(seq_2) ;
    kinetic_2.setSequenceRv(seq_1) ;
    kinetic_2.setIPDFw(ipd_2) ;
    kinetic_2.setIPDRv(ipd_1) ;
    kinetic_2.setPWDFw(pwd_2) ;
    kinetic_2.setPWDRv(pwd_1) ;
    ASSERT_EQ(kinetic_2.size(), seq_1.size()) ;
    ASSERT_THAT(kinetic_2.getSequenceFw(), ContainerEq(seq_2)) ;
    ASSERT_THAT(kinetic_2.getSequenceRv(), ContainerEq(seq_1)) ;
    ASSERT_THAT(kinetic_2.getIPDFw(), ContainerEq(ipd_2)) ;
    ASSERT_THAT(kinetic_2.getIPDRv(), ContainerEq(ipd_1)) ;
    ASSERT_THAT(kinetic_2.getPWDFw(), ContainerEq(pwd_2)) ;
    ASSERT_THAT(kinetic_2.getPWDRv(), ContainerEq(pwd_1)) ;     
    // modify field values with diff length values ->  must fail
    EXPECT_THROW(kinetic_2.setSequenceFw(seq_0), std::invalid_argument) ;
    EXPECT_THROW(kinetic_2.setSequenceFw(seq_3), std::invalid_argument) ;
    EXPECT_THROW(kinetic_2.setSequenceRv(seq_0), std::invalid_argument) ;
    EXPECT_THROW(kinetic_2.setSequenceRv(seq_3), std::invalid_argument) ;
    EXPECT_THROW(kinetic_2.setIPDFw(ipd_0),      std::invalid_argument) ;
    EXPECT_THROW(kinetic_2.setIPDFw(ipd_3),      std::invalid_argument) ;
    EXPECT_THROW(kinetic_2.setIPDRv(ipd_0),      std::invalid_argument) ;
    EXPECT_THROW(kinetic_2.setIPDRv(ipd_3),      std::invalid_argument) ;
    EXPECT_THROW(kinetic_2.setPWDFw(pwd_0),      std::invalid_argument) ;
    EXPECT_THROW(kinetic_2.setPWDFw(pwd_3),      std::invalid_argument) ;
    EXPECT_THROW(kinetic_2.setPWDRv(pwd_0),      std::invalid_argument) ;
    EXPECT_THROW(kinetic_2.setPWDRv(pwd_3),      std::invalid_argument) ;
}


/*!
 * \brief Constructs an empty instance and populate the ipd fw first, 
 * then the other fields and checks that it behaves properly.
 */
TEST(KineticSignalTest, ipd_fw)
{   
    std::string seq_0 ;
    std::string seq_1("ACGT") ;
    std::string seq_2("CGCG") ;
    std::string seq_3("ACGTACGTACGT") ;

    std::vector<double> ipd_0 ;
    std::vector<double> ipd_1({1.,2.,3.,4.}) ;
    std::vector<double> ipd_2({5.,6.,7.,8.}) ;
    std::vector<double> ipd_3({1.,2.,3.,4.,5.,6.,7.,8.}) ;

    std::vector<double> pwd_0 ;
    std::vector<double> pwd_1({1.,2.,3.,4.}) ;
    std::vector<double> pwd_2({5.,6.,7.,8.}) ;
    std::vector<double> pwd_3({1.,2.,3.,4.,5.,6.,7.,8.}) ;

    // init using setIPDFw by reference
    ngsai::KineticSignal kinetic_1 ;
    kinetic_1.setIPDFw(ipd_1) ; // sets size from ipd fw
    kinetic_1.setSequenceFw(seq_1) ;
    kinetic_1.setSequenceRv(seq_2) ;
    kinetic_1.setIPDRv(ipd_2) ;
    kinetic_1.setPWDFw(pwd_1) ;
    kinetic_1.setPWDRv(pwd_2) ;
    ASSERT_EQ(kinetic_1.size(), seq_1.size()) ;
    ASSERT_THAT(kinetic_1.getSequenceFw(), ContainerEq(seq_1)) ;
    ASSERT_THAT(kinetic_1.getSequenceRv(), ContainerEq(seq_2)) ;
    ASSERT_THAT(kinetic_1.getIPDFw(), ContainerEq(ipd_1)) ;
    ASSERT_THAT(kinetic_1.getIPDRv(), ContainerEq(ipd_2)) ;
    ASSERT_THAT(kinetic_1.getPWDFw(), ContainerEq(pwd_1)) ;
    ASSERT_THAT(kinetic_1.getPWDRv(), ContainerEq(pwd_2)) ;
    // modify field values with same length values -> must work
    kinetic_1.setSequenceFw(seq_2) ;
    kinetic_1.setSequenceRv(seq_1) ;
    kinetic_1.setIPDFw(ipd_2) ;
    kinetic_1.setIPDRv(ipd_1) ;
    kinetic_1.setPWDFw(pwd_2) ;
    kinetic_1.setPWDRv(pwd_1) ;
    ASSERT_EQ(kinetic_1.size(), seq_1.size()) ;
    ASSERT_THAT(kinetic_1.getSequenceFw(), ContainerEq(seq_2)) ;
    ASSERT_THAT(kinetic_1.getSequenceRv(), ContainerEq(seq_1)) ;
    ASSERT_THAT(kinetic_1.getIPDFw(), ContainerEq(ipd_2)) ;
    ASSERT_THAT(kinetic_1.getIPDRv(), ContainerEq(ipd_1)) ;
    ASSERT_THAT(kinetic_1.getPWDFw(), ContainerEq(pwd_2)) ;
    ASSERT_THAT(kinetic_1.getPWDRv(), ContainerEq(pwd_1)) ;     
    // modify field values with diff length values ->  must fail
    EXPECT_THROW(kinetic_1.setSequenceFw(seq_0), std::invalid_argument) ;
    EXPECT_THROW(kinetic_1.setSequenceFw(seq_3), std::invalid_argument) ;
    EXPECT_THROW(kinetic_1.setSequenceRv(seq_0), std::invalid_argument) ;
    EXPECT_THROW(kinetic_1.setSequenceRv(seq_3), std::invalid_argument) ;
    EXPECT_THROW(kinetic_1.setIPDFw(ipd_0),      std::invalid_argument) ;
    EXPECT_THROW(kinetic_1.setIPDFw(ipd_3),      std::invalid_argument) ;
    EXPECT_THROW(kinetic_1.setIPDRv(ipd_0),      std::invalid_argument) ;
    EXPECT_THROW(kinetic_1.setIPDRv(ipd_3),      std::invalid_argument) ;
    EXPECT_THROW(kinetic_1.setPWDFw(pwd_0),      std::invalid_argument) ;
    EXPECT_THROW(kinetic_1.setPWDFw(pwd_3),      std::invalid_argument) ;
    EXPECT_THROW(kinetic_1.setPWDRv(pwd_0),      std::invalid_argument) ;
    EXPECT_THROW(kinetic_1.setPWDRv(pwd_3),      std::invalid_argument) ;


    // init using setIPDFw by rvalue
    ngsai::KineticSignal kinetic_2 ;
    kinetic_2.setIPDFw(std::vector<double>(ipd_1)) ; // sets size from ipd fw
    kinetic_2.setSequenceFw(seq_1) ;
    kinetic_2.setSequenceRv(seq_2) ;
    kinetic_2.setIPDRv(ipd_2) ;
    kinetic_2.setPWDFw(pwd_1) ;
    kinetic_2.setPWDRv(pwd_2) ;
    ASSERT_EQ(kinetic_2.size(), seq_1.size()) ;
    ASSERT_THAT(kinetic_2.getSequenceFw(), ContainerEq(seq_1)) ;
    ASSERT_THAT(kinetic_2.getSequenceRv(), ContainerEq(seq_2)) ;
    ASSERT_THAT(kinetic_2.getIPDFw(), ContainerEq(ipd_1)) ;
    ASSERT_THAT(kinetic_2.getIPDRv(), ContainerEq(ipd_2)) ;
    ASSERT_THAT(kinetic_2.getPWDFw(), ContainerEq(pwd_1)) ;
    ASSERT_THAT(kinetic_2.getPWDRv(), ContainerEq(pwd_2)) ;
    // modify field values with same length values -> must work
    kinetic_2.setSequenceFw(seq_2) ;
    kinetic_2.setSequenceRv(seq_1) ;
    kinetic_2.setIPDFw(ipd_2) ;
    kinetic_2.setIPDRv(ipd_1) ;
    kinetic_2.setPWDFw(pwd_2) ;
    kinetic_2.setPWDRv(pwd_1) ;
    ASSERT_EQ(kinetic_2.size(), seq_1.size()) ;
    ASSERT_THAT(kinetic_2.getSequenceFw(), ContainerEq(seq_2)) ;
    ASSERT_THAT(kinetic_2.getSequenceRv(), ContainerEq(seq_1)) ;
    ASSERT_THAT(kinetic_2.getIPDFw(), ContainerEq(ipd_2)) ;
    ASSERT_THAT(kinetic_2.getIPDRv(), ContainerEq(ipd_1)) ;
    ASSERT_THAT(kinetic_2.getPWDFw(), ContainerEq(pwd_2)) ;
    ASSERT_THAT(kinetic_2.getPWDRv(), ContainerEq(pwd_1)) ;     
    // modify field values with diff length values ->  must fail
    EXPECT_THROW(kinetic_2.setSequenceFw(seq_0), std::invalid_argument) ;
    EXPECT_THROW(kinetic_2.setSequenceFw(seq_3), std::invalid_argument) ;
    EXPECT_THROW(kinetic_2.setSequenceRv(seq_0), std::invalid_argument) ;
    EXPECT_THROW(kinetic_2.setSequenceRv(seq_3), std::invalid_argument) ;
    EXPECT_THROW(kinetic_2.setIPDFw(ipd_0),      std::invalid_argument) ;
    EXPECT_THROW(kinetic_2.setIPDFw(ipd_3),      std::invalid_argument) ;
    EXPECT_THROW(kinetic_2.setIPDRv(ipd_0),      std::invalid_argument) ;
    EXPECT_THROW(kinetic_2.setIPDRv(ipd_3),      std::invalid_argument) ;
    EXPECT_THROW(kinetic_2.setPWDFw(pwd_0),      std::invalid_argument) ;
    EXPECT_THROW(kinetic_2.setPWDFw(pwd_3),      std::invalid_argument) ;
    EXPECT_THROW(kinetic_2.setPWDRv(pwd_0),      std::invalid_argument) ;
    EXPECT_THROW(kinetic_2.setPWDRv(pwd_3),      std::invalid_argument) ;
}


/*!
 * \brief Constructs an empty instance and populate the ipd rv first, 
 * then the other fields and checks that it behaves properly.
 */
TEST(KineticSignalTest, ipd_rv)
{   
    std::string seq_0 ;
    std::string seq_1("ACGT") ;
    std::string seq_2("CGCG") ;
    std::string seq_3("ACGTACGTACGT") ;

    std::vector<double> ipd_0 ;
    std::vector<double> ipd_1({1.,2.,3.,4.}) ;
    std::vector<double> ipd_2({5.,6.,7.,8.}) ;
    std::vector<double> ipd_3({1.,2.,3.,4.,5.,6.,7.,8.}) ;

    std::vector<double> pwd_0 ;
    std::vector<double> pwd_1({1.,2.,3.,4.}) ;
    std::vector<double> pwd_2({5.,6.,7.,8.}) ;
    std::vector<double> pwd_3({1.,2.,3.,4.,5.,6.,7.,8.}) ;

    // init using setIPDRv by reference
    ngsai::KineticSignal kinetic_1 ;
    kinetic_1.setIPDRv(ipd_2) ; // sets size from ipd rv
    kinetic_1.setSequenceFw(seq_1) ;
    kinetic_1.setSequenceRv(seq_2) ;
    kinetic_1.setIPDFw(ipd_1) ;
    kinetic_1.setPWDFw(pwd_1) ;
    kinetic_1.setPWDRv(pwd_2) ;
    ASSERT_EQ(kinetic_1.size(), seq_1.size()) ;
    ASSERT_THAT(kinetic_1.getSequenceFw(), ContainerEq(seq_1)) ;
    ASSERT_THAT(kinetic_1.getSequenceRv(), ContainerEq(seq_2)) ;
    ASSERT_THAT(kinetic_1.getIPDFw(), ContainerEq(ipd_1)) ;
    ASSERT_THAT(kinetic_1.getIPDRv(), ContainerEq(ipd_2)) ;
    ASSERT_THAT(kinetic_1.getPWDFw(), ContainerEq(pwd_1)) ;
    ASSERT_THAT(kinetic_1.getPWDRv(), ContainerEq(pwd_2)) ;
    // modify field values with same length values -> must work
    kinetic_1.setSequenceFw(seq_2) ;
    kinetic_1.setSequenceRv(seq_1) ;
    kinetic_1.setIPDFw(ipd_2) ;
    kinetic_1.setIPDRv(ipd_1) ;
    kinetic_1.setPWDFw(pwd_2) ;
    kinetic_1.setPWDRv(pwd_1) ;
    ASSERT_EQ(kinetic_1.size(), seq_1.size()) ;
    ASSERT_THAT(kinetic_1.getSequenceFw(), ContainerEq(seq_2)) ;
    ASSERT_THAT(kinetic_1.getSequenceRv(), ContainerEq(seq_1)) ;
    ASSERT_THAT(kinetic_1.getIPDFw(), ContainerEq(ipd_2)) ;
    ASSERT_THAT(kinetic_1.getIPDRv(), ContainerEq(ipd_1)) ;
    ASSERT_THAT(kinetic_1.getPWDFw(), ContainerEq(pwd_2)) ;
    ASSERT_THAT(kinetic_1.getPWDRv(), ContainerEq(pwd_1)) ;     
    // modify field values with diff length values ->  must fail
    EXPECT_THROW(kinetic_1.setSequenceFw(seq_0), std::invalid_argument) ;
    EXPECT_THROW(kinetic_1.setSequenceFw(seq_3), std::invalid_argument) ;
    EXPECT_THROW(kinetic_1.setSequenceRv(seq_0), std::invalid_argument) ;
    EXPECT_THROW(kinetic_1.setSequenceRv(seq_3), std::invalid_argument) ;
    EXPECT_THROW(kinetic_1.setIPDFw(ipd_0),      std::invalid_argument) ;
    EXPECT_THROW(kinetic_1.setIPDFw(ipd_3),      std::invalid_argument) ;
    EXPECT_THROW(kinetic_1.setIPDRv(ipd_0),      std::invalid_argument) ;
    EXPECT_THROW(kinetic_1.setIPDRv(ipd_3),      std::invalid_argument) ;
    EXPECT_THROW(kinetic_1.setPWDFw(pwd_0),      std::invalid_argument) ;
    EXPECT_THROW(kinetic_1.setPWDFw(pwd_3),      std::invalid_argument) ;
    EXPECT_THROW(kinetic_1.setPWDRv(pwd_0),      std::invalid_argument) ;
    EXPECT_THROW(kinetic_1.setPWDRv(pwd_3),      std::invalid_argument) ;


    // init using setIPDRv by rvalue
    ngsai::KineticSignal kinetic_2 ;
    kinetic_2.setIPDRv(std::vector<double>(ipd_2)) ; // sets size from ipd rv
    kinetic_2.setSequenceFw(seq_1) ;
    kinetic_2.setSequenceRv(seq_2) ;
    kinetic_2.setIPDFw(ipd_1) ;
    kinetic_2.setPWDFw(pwd_1) ;
    kinetic_2.setPWDRv(pwd_2) ;
    ASSERT_EQ(kinetic_2.size(), seq_1.size()) ;
    ASSERT_THAT(kinetic_2.getSequenceFw(), ContainerEq(seq_1)) ;
    ASSERT_THAT(kinetic_2.getSequenceRv(), ContainerEq(seq_2)) ;
    ASSERT_THAT(kinetic_2.getIPDFw(), ContainerEq(ipd_1)) ;
    ASSERT_THAT(kinetic_2.getIPDRv(), ContainerEq(ipd_2)) ;
    ASSERT_THAT(kinetic_2.getPWDFw(), ContainerEq(pwd_1)) ;
    ASSERT_THAT(kinetic_2.getPWDRv(), ContainerEq(pwd_2)) ;
    // modify field values with same length values -> must work
    kinetic_2.setSequenceFw(seq_2) ;
    kinetic_2.setSequenceRv(seq_1) ;
    kinetic_2.setIPDFw(ipd_2) ;
    kinetic_2.setIPDRv(ipd_1) ;
    kinetic_2.setPWDFw(pwd_2) ;
    kinetic_2.setPWDRv(pwd_1) ;
    ASSERT_EQ(kinetic_2.size(), seq_1.size()) ;
    ASSERT_THAT(kinetic_2.getSequenceFw(), ContainerEq(seq_2)) ;
    ASSERT_THAT(kinetic_2.getSequenceRv(), ContainerEq(seq_1)) ;
    ASSERT_THAT(kinetic_2.getIPDFw(), ContainerEq(ipd_2)) ;
    ASSERT_THAT(kinetic_2.getIPDRv(), ContainerEq(ipd_1)) ;
    ASSERT_THAT(kinetic_2.getPWDFw(), ContainerEq(pwd_2)) ;
    ASSERT_THAT(kinetic_2.getPWDRv(), ContainerEq(pwd_1)) ;     
    // modify field values with diff length values ->  must fail
    EXPECT_THROW(kinetic_2.setSequenceFw(seq_0), std::invalid_argument) ;
    EXPECT_THROW(kinetic_2.setSequenceFw(seq_3), std::invalid_argument) ;
    EXPECT_THROW(kinetic_2.setSequenceRv(seq_0), std::invalid_argument) ;
    EXPECT_THROW(kinetic_2.setSequenceRv(seq_3), std::invalid_argument) ;
    EXPECT_THROW(kinetic_2.setIPDFw(ipd_0),      std::invalid_argument) ;
    EXPECT_THROW(kinetic_2.setIPDFw(ipd_3),      std::invalid_argument) ;
    EXPECT_THROW(kinetic_2.setIPDRv(ipd_0),      std::invalid_argument) ;
    EXPECT_THROW(kinetic_2.setIPDRv(ipd_3),      std::invalid_argument) ;
    EXPECT_THROW(kinetic_2.setPWDFw(pwd_0),      std::invalid_argument) ;
    EXPECT_THROW(kinetic_2.setPWDFw(pwd_3),      std::invalid_argument) ;
    EXPECT_THROW(kinetic_2.setPWDRv(pwd_0),      std::invalid_argument) ;
    EXPECT_THROW(kinetic_2.setPWDRv(pwd_3),      std::invalid_argument) ;
}


/*!
 * \brief Constructs an empty instance and populate the pwd fw first, 
 * then the other fields and checks that it behaves properly.
 */
TEST(KineticSignalTest, pwd_fw)
{   
    std::string seq_0 ;
    std::string seq_1("ACGT") ;
    std::string seq_2("CGCG") ;
    std::string seq_3("ACGTACGTACGT") ;

    std::vector<double> ipd_0 ;
    std::vector<double> ipd_1({1.,2.,3.,4.}) ;
    std::vector<double> ipd_2({5.,6.,7.,8.}) ;
    std::vector<double> ipd_3({1.,2.,3.,4.,5.,6.,7.,8.}) ;

    std::vector<double> pwd_0 ;
    std::vector<double> pwd_1({1.,2.,3.,4.}) ;
    std::vector<double> pwd_2({5.,6.,7.,8.}) ;
    std::vector<double> pwd_3({1.,2.,3.,4.,5.,6.,7.,8.}) ;

    // init using setPWDFw by reference
    ngsai::KineticSignal kinetic_1 ;
    kinetic_1.setPWDFw(pwd_1) ; // sets size from pwd fw
    kinetic_1.setSequenceFw(seq_1) ;
    kinetic_1.setSequenceRv(seq_2) ;
    kinetic_1.setIPDFw(ipd_1) ;
    kinetic_1.setIPDRv(ipd_2) ;
    kinetic_1.setPWDRv(pwd_2) ;
    ASSERT_EQ(kinetic_1.size(), seq_1.size()) ;
    ASSERT_THAT(kinetic_1.getSequenceFw(), ContainerEq(seq_1)) ;
    ASSERT_THAT(kinetic_1.getSequenceRv(), ContainerEq(seq_2)) ;
    ASSERT_THAT(kinetic_1.getIPDFw(), ContainerEq(ipd_1)) ;
    ASSERT_THAT(kinetic_1.getIPDRv(), ContainerEq(ipd_2)) ;
    ASSERT_THAT(kinetic_1.getPWDFw(), ContainerEq(pwd_1)) ;
    ASSERT_THAT(kinetic_1.getPWDRv(), ContainerEq(pwd_2)) ;
    // modify field values with same length values -> must work
    kinetic_1.setSequenceFw(seq_2) ;
    kinetic_1.setSequenceRv(seq_1) ;
    kinetic_1.setIPDFw(ipd_2) ;
    kinetic_1.setIPDRv(ipd_1) ;
    kinetic_1.setPWDFw(pwd_2) ;
    kinetic_1.setPWDRv(pwd_1) ;
    ASSERT_EQ(kinetic_1.size(), seq_1.size()) ;
    ASSERT_THAT(kinetic_1.getSequenceFw(), ContainerEq(seq_2)) ;
    ASSERT_THAT(kinetic_1.getSequenceRv(), ContainerEq(seq_1)) ;
    ASSERT_THAT(kinetic_1.getIPDFw(), ContainerEq(ipd_2)) ;
    ASSERT_THAT(kinetic_1.getIPDRv(), ContainerEq(ipd_1)) ;
    ASSERT_THAT(kinetic_1.getPWDFw(), ContainerEq(pwd_2)) ;
    ASSERT_THAT(kinetic_1.getPWDRv(), ContainerEq(pwd_1)) ;     
    // modify field values with diff length values ->  must fail
    EXPECT_THROW(kinetic_1.setSequenceFw(seq_0), std::invalid_argument) ;
    EXPECT_THROW(kinetic_1.setSequenceFw(seq_3), std::invalid_argument) ;
    EXPECT_THROW(kinetic_1.setSequenceRv(seq_0), std::invalid_argument) ;
    EXPECT_THROW(kinetic_1.setSequenceRv(seq_3), std::invalid_argument) ;
    EXPECT_THROW(kinetic_1.setIPDFw(ipd_0),      std::invalid_argument) ;
    EXPECT_THROW(kinetic_1.setIPDFw(ipd_3),      std::invalid_argument) ;
    EXPECT_THROW(kinetic_1.setIPDRv(ipd_0),      std::invalid_argument) ;
    EXPECT_THROW(kinetic_1.setIPDRv(ipd_3),      std::invalid_argument) ;
    EXPECT_THROW(kinetic_1.setPWDFw(pwd_0),      std::invalid_argument) ;
    EXPECT_THROW(kinetic_1.setPWDFw(pwd_3),      std::invalid_argument) ;
    EXPECT_THROW(kinetic_1.setPWDRv(pwd_0),      std::invalid_argument) ;
    EXPECT_THROW(kinetic_1.setPWDRv(pwd_3),      std::invalid_argument) ;


    // init using setPWDFw by rvalue
    ngsai::KineticSignal kinetic_2 ;
    kinetic_2.setPWDFw(std::vector<double>(pwd_1)) ; // sets size from pwd fw
    kinetic_2.setSequenceFw(seq_1) ;
    kinetic_2.setSequenceRv(seq_2) ;
    kinetic_2.setIPDFw(ipd_1) ;
    kinetic_2.setIPDRv(ipd_2) ;
    kinetic_2.setPWDRv(pwd_2) ;
    ASSERT_EQ(kinetic_2.size(), seq_1.size()) ;
    ASSERT_THAT(kinetic_2.getSequenceFw(), ContainerEq(seq_1)) ;
    ASSERT_THAT(kinetic_2.getSequenceRv(), ContainerEq(seq_2)) ;
    ASSERT_THAT(kinetic_2.getIPDFw(), ContainerEq(ipd_1)) ;
    ASSERT_THAT(kinetic_2.getIPDRv(), ContainerEq(ipd_2)) ;
    ASSERT_THAT(kinetic_2.getPWDFw(), ContainerEq(pwd_1)) ;
    ASSERT_THAT(kinetic_2.getPWDRv(), ContainerEq(pwd_2)) ;
    // modify field values with same length values -> must work
    kinetic_2.setSequenceFw(seq_2) ;
    kinetic_2.setSequenceRv(seq_1) ;
    kinetic_2.setIPDFw(ipd_2) ;
    kinetic_2.setIPDRv(ipd_1) ;
    kinetic_2.setPWDFw(pwd_2) ;
    kinetic_2.setPWDRv(pwd_1) ;
    ASSERT_EQ(kinetic_2.size(), seq_1.size()) ;
    ASSERT_THAT(kinetic_2.getSequenceFw(), ContainerEq(seq_2)) ;
    ASSERT_THAT(kinetic_2.getSequenceRv(), ContainerEq(seq_1)) ;
    ASSERT_THAT(kinetic_2.getIPDFw(), ContainerEq(ipd_2)) ;
    ASSERT_THAT(kinetic_2.getIPDRv(), ContainerEq(ipd_1)) ;
    ASSERT_THAT(kinetic_2.getPWDFw(), ContainerEq(pwd_2)) ;
    ASSERT_THAT(kinetic_2.getPWDRv(), ContainerEq(pwd_1)) ;     
    // modify field values with diff length values ->  must fail
    EXPECT_THROW(kinetic_2.setSequenceFw(seq_0), std::invalid_argument) ;
    EXPECT_THROW(kinetic_2.setSequenceFw(seq_3), std::invalid_argument) ;
    EXPECT_THROW(kinetic_2.setSequenceRv(seq_0), std::invalid_argument) ;
    EXPECT_THROW(kinetic_2.setSequenceRv(seq_3), std::invalid_argument) ;
    EXPECT_THROW(kinetic_2.setIPDFw(ipd_0),      std::invalid_argument) ;
    EXPECT_THROW(kinetic_2.setIPDFw(ipd_3),      std::invalid_argument) ;
    EXPECT_THROW(kinetic_2.setIPDRv(ipd_0),      std::invalid_argument) ;
    EXPECT_THROW(kinetic_2.setIPDRv(ipd_3),      std::invalid_argument) ;
    EXPECT_THROW(kinetic_2.setPWDFw(pwd_0),      std::invalid_argument) ;
    EXPECT_THROW(kinetic_2.setPWDFw(pwd_3),      std::invalid_argument) ;
    EXPECT_THROW(kinetic_2.setPWDRv(pwd_0),      std::invalid_argument) ;
    EXPECT_THROW(kinetic_2.setPWDRv(pwd_3),      std::invalid_argument) ;
}


/*!
 * \brief Constructs an empty instance and populate the pwd rv first, 
 * then the other fields and checks that it behaves properly.
 */
TEST(KineticSignalTest, pwd_rv)
{   
    std::string seq_0 ;
    std::string seq_1("ACGT") ;
    std::string seq_2("CGCG") ;
    std::string seq_3("ACGTACGTACGT") ;

    std::vector<double> ipd_0 ;
    std::vector<double> ipd_1({1.,2.,3.,4.}) ;
    std::vector<double> ipd_2({5.,6.,7.,8.}) ;
    std::vector<double> ipd_3({1.,2.,3.,4.,5.,6.,7.,8.}) ;

    std::vector<double> pwd_0 ;
    std::vector<double> pwd_1({1.,2.,3.,4.}) ;
    std::vector<double> pwd_2({5.,6.,7.,8.}) ;
    std::vector<double> pwd_3({1.,2.,3.,4.,5.,6.,7.,8.}) ;

    // init using setPWDRv by reference
    ngsai::KineticSignal kinetic_1 ;
    kinetic_1.setPWDRv(pwd_2) ; // sets size from pwd rv
    kinetic_1.setSequenceFw(seq_1) ;
    kinetic_1.setSequenceRv(seq_2) ;
    kinetic_1.setIPDFw(ipd_1) ;
    kinetic_1.setIPDRv(ipd_2) ;
    kinetic_1.setPWDFw(pwd_1) ;
    ASSERT_EQ(kinetic_1.size(), seq_1.size()) ;
    ASSERT_THAT(kinetic_1.getSequenceFw(), ContainerEq(seq_1)) ;
    ASSERT_THAT(kinetic_1.getSequenceRv(), ContainerEq(seq_2)) ;
    ASSERT_THAT(kinetic_1.getIPDFw(), ContainerEq(ipd_1)) ;
    ASSERT_THAT(kinetic_1.getIPDRv(), ContainerEq(ipd_2)) ;
    ASSERT_THAT(kinetic_1.getPWDFw(), ContainerEq(pwd_1)) ;
    ASSERT_THAT(kinetic_1.getPWDRv(), ContainerEq(pwd_2)) ;
    // modify field values with same length values -> must work
    kinetic_1.setSequenceFw(seq_2) ;
    kinetic_1.setSequenceRv(seq_1) ;
    kinetic_1.setIPDFw(ipd_2) ;
    kinetic_1.setIPDRv(ipd_1) ;
    kinetic_1.setPWDFw(pwd_2) ;
    kinetic_1.setPWDRv(pwd_1) ;
    ASSERT_EQ(kinetic_1.size(), seq_1.size()) ;
    ASSERT_THAT(kinetic_1.getSequenceFw(), ContainerEq(seq_2)) ;
    ASSERT_THAT(kinetic_1.getSequenceRv(), ContainerEq(seq_1)) ;
    ASSERT_THAT(kinetic_1.getIPDFw(), ContainerEq(ipd_2)) ;
    ASSERT_THAT(kinetic_1.getIPDRv(), ContainerEq(ipd_1)) ;
    ASSERT_THAT(kinetic_1.getPWDFw(), ContainerEq(pwd_2)) ;
    ASSERT_THAT(kinetic_1.getPWDRv(), ContainerEq(pwd_1)) ;     
    // modify field values with diff length values ->  must fail
    EXPECT_THROW(kinetic_1.setSequenceFw(seq_0), std::invalid_argument) ;
    EXPECT_THROW(kinetic_1.setSequenceFw(seq_3), std::invalid_argument) ;
    EXPECT_THROW(kinetic_1.setSequenceRv(seq_0), std::invalid_argument) ;
    EXPECT_THROW(kinetic_1.setSequenceRv(seq_3), std::invalid_argument) ;
    EXPECT_THROW(kinetic_1.setIPDFw(ipd_0),      std::invalid_argument) ;
    EXPECT_THROW(kinetic_1.setIPDFw(ipd_3),      std::invalid_argument) ;
    EXPECT_THROW(kinetic_1.setIPDRv(ipd_0),      std::invalid_argument) ;
    EXPECT_THROW(kinetic_1.setIPDRv(ipd_3),      std::invalid_argument) ;
    EXPECT_THROW(kinetic_1.setPWDFw(pwd_0),      std::invalid_argument) ;
    EXPECT_THROW(kinetic_1.setPWDFw(pwd_3),      std::invalid_argument) ;
    EXPECT_THROW(kinetic_1.setPWDRv(pwd_0),      std::invalid_argument) ;
    EXPECT_THROW(kinetic_1.setPWDRv(pwd_3),      std::invalid_argument) ;


    // init using setPWDRv by rvalue
    ngsai::KineticSignal kinetic_2 ;
    kinetic_2.setPWDRv(std::vector<double>(pwd_2)) ; // sets size from pwd rv
    kinetic_2.setSequenceFw(seq_1) ;
    kinetic_2.setSequenceRv(seq_2) ;
    kinetic_2.setIPDFw(ipd_1) ;
    kinetic_2.setIPDRv(ipd_2) ;
    kinetic_2.setPWDFw(pwd_1) ;
    ASSERT_EQ(kinetic_2.size(), seq_1.size()) ;
    ASSERT_THAT(kinetic_2.getSequenceFw(), ContainerEq(seq_1)) ;
    ASSERT_THAT(kinetic_2.getSequenceRv(), ContainerEq(seq_2)) ;
    ASSERT_THAT(kinetic_2.getIPDFw(), ContainerEq(ipd_1)) ;
    ASSERT_THAT(kinetic_2.getIPDRv(), ContainerEq(ipd_2)) ;
    ASSERT_THAT(kinetic_2.getPWDFw(), ContainerEq(pwd_1)) ;
    ASSERT_THAT(kinetic_2.getPWDRv(), ContainerEq(pwd_2)) ;
    // modify field values with same length values -> must work
    kinetic_2.setSequenceFw(seq_2) ;
    kinetic_2.setSequenceRv(seq_1) ;
    kinetic_2.setIPDFw(ipd_2) ;
    kinetic_2.setIPDRv(ipd_1) ;
    kinetic_2.setPWDFw(pwd_2) ;
    kinetic_2.setPWDRv(pwd_1) ;
    ASSERT_EQ(kinetic_2.size(), seq_1.size()) ;
    ASSERT_THAT(kinetic_2.getSequenceFw(), ContainerEq(seq_2)) ;
    ASSERT_THAT(kinetic_2.getSequenceRv(), ContainerEq(seq_1)) ;
    ASSERT_THAT(kinetic_2.getIPDFw(), ContainerEq(ipd_2)) ;
    ASSERT_THAT(kinetic_2.getIPDRv(), ContainerEq(ipd_1)) ;
    ASSERT_THAT(kinetic_2.getPWDFw(), ContainerEq(pwd_2)) ;
    ASSERT_THAT(kinetic_2.getPWDRv(), ContainerEq(pwd_1)) ;     
    // modify field values with diff length values ->  must fail
    EXPECT_THROW(kinetic_2.setSequenceFw(seq_0), std::invalid_argument) ;
    EXPECT_THROW(kinetic_2.setSequenceFw(seq_3), std::invalid_argument) ;
    EXPECT_THROW(kinetic_2.setSequenceRv(seq_0), std::invalid_argument) ;
    EXPECT_THROW(kinetic_2.setSequenceRv(seq_3), std::invalid_argument) ;
    EXPECT_THROW(kinetic_2.setIPDFw(ipd_0),      std::invalid_argument) ;
    EXPECT_THROW(kinetic_2.setIPDFw(ipd_3),      std::invalid_argument) ;
    EXPECT_THROW(kinetic_2.setIPDRv(ipd_0),      std::invalid_argument) ;
    EXPECT_THROW(kinetic_2.setIPDRv(ipd_3),      std::invalid_argument) ;
    EXPECT_THROW(kinetic_2.setPWDFw(pwd_0),      std::invalid_argument) ;
    EXPECT_THROW(kinetic_2.setPWDFw(pwd_3),      std::invalid_argument) ;
    EXPECT_THROW(kinetic_2.setPWDRv(pwd_0),      std::invalid_argument) ;
    EXPECT_THROW(kinetic_2.setPWDRv(pwd_3),      std::invalid_argument) ;
}


/*!
 * \brief Checks that the getIPDMean method works properly.
 */
TEST(KineticSignalTest, getIPDMean)
{
    std::vector<double> ipd_1({1.,2.,3.,4.}) ;
    std::vector<double> ipd_2({5.,6.,7.,8.}) ;

    // starts by seting fw strand
    ngsai::KineticSignal kinetics_1 ;
    // no IPD value yet
    ASSERT_THAT(kinetics_1.getIPDMean(), ContainerEq(std::vector<double>())) ;
    // set fw only
    kinetics_1.setIPDFw(ipd_1) ;
    ASSERT_THAT(kinetics_1.getIPDMean(), ContainerEq(ipd_1)) ;
    // set rv also
    kinetics_1.setIPDRv(ipd_2) ;
    ASSERT_THAT(kinetics_1.getIPDMean(),
                ContainerEq(std::vector<double>({3.,4.,5.,6.}))) ;
    // change fw
    kinetics_1.setIPDFw(ipd_2) ;
    ASSERT_THAT(kinetics_1.getIPDMean(), ContainerEq(ipd_2)) ;
    // change rv
    kinetics_1.setIPDFw(ipd_1) ;
    ASSERT_THAT(kinetics_1.getIPDMean(),
                ContainerEq(std::vector<double>({3.,4.,5.,6.}))) ;

    
    // starts by seting fw strand
    ngsai::KineticSignal kinetics_2 ;
    
    // set rv only
    kinetics_2.setIPDRv(ipd_2) ;
    ASSERT_THAT(kinetics_2.getIPDMean(), ContainerEq(ipd_2)) ;
    // set rv also
    kinetics_2.setIPDFw(ipd_1) ;
    ASSERT_THAT(kinetics_2.getIPDMean(),
                ContainerEq(std::vector<double>({3.,4.,5.,6.}))) ;
    // change rv
    kinetics_2.setIPDRv(ipd_1) ;
    ASSERT_THAT(kinetics_2.getIPDMean(), ContainerEq(ipd_1)) ;
    // change fw
    kinetics_2.setIPDFw(ipd_2) ;
    ASSERT_THAT(kinetics_1.getIPDMean(),
                ContainerEq(std::vector<double>({3.,4.,5.,6.}))) ;

}


/*!
 * \brief Checks that the getPWDMean method works properly.
 */
TEST(KineticSignalTest, getPWDMean)
{
    std::vector<double> pwd_1({1.,2.,3.,4.}) ;
    std::vector<double> pwd_2({5.,6.,7.,8.}) ;

    // starts by seting fw strand
    ngsai::KineticSignal kinetics_1 ;
    // no PWD value yet
    ASSERT_THAT(kinetics_1.getPWDMean(), ContainerEq(std::vector<double>())) ;
    // set fw only
    kinetics_1.setPWDFw(pwd_1) ;
    ASSERT_THAT(kinetics_1.getPWDMean(), ContainerEq(pwd_1)) ;
    // set rv also
    kinetics_1.setPWDRv(pwd_2) ;
    ASSERT_THAT(kinetics_1.getPWDMean(),
                ContainerEq(std::vector<double>({3.,4.,5.,6.}))) ;
    // change fw
    kinetics_1.setPWDFw(pwd_2) ;
    ASSERT_THAT(kinetics_1.getPWDMean(), ContainerEq(pwd_2)) ;
    // change rv
    kinetics_1.setPWDFw(pwd_1) ;
    ASSERT_THAT(kinetics_1.getPWDMean(),
                ContainerEq(std::vector<double>({3.,4.,5.,6.}))) ;

    
    // starts by seting fw strand
    ngsai::KineticSignal kinetics_2 ;
    
    // set rv only
    kinetics_2.setPWDRv(pwd_2) ;
    ASSERT_THAT(kinetics_2.getPWDMean(), ContainerEq(pwd_2)) ;
    // set rv also
    kinetics_2.setPWDFw(pwd_1) ;
    ASSERT_THAT(kinetics_2.getPWDMean(),
                ContainerEq(std::vector<double>({3.,4.,5.,6.}))) ;
    // change rv
    kinetics_2.setPWDRv(pwd_1) ;
    ASSERT_THAT(kinetics_2.getPWDMean(), ContainerEq(pwd_1)) ;
    // change fw
    kinetics_2.setPWDFw(pwd_2) ;
    ASSERT_THAT(kinetics_1.getPWDMean(),
                ContainerEq(std::vector<double>({3.,4.,5.,6.}))) ;

}


/*!
 * \brief Checks that the isFwComplete method works properly.
 */
TEST(KineticSignalTest, isFwComplete)
{   
    std::string s("AAAAA") ;
    std::vector<double> k({1., 2., 3., 4., 5.}) ;
    
    ngsai::KineticSignal s_fw ;

    // setting IPD then PWD then seq
    s_fw = ngsai::KineticSignal() ;
    ASSERT_EQ(s_fw.isFwComplete(), false) ;
    s_fw.setIPDFw(k)      ; ASSERT_EQ(s_fw.isFwComplete(), false) ;
    s_fw.setPWDFw(k)      ; ASSERT_EQ(s_fw.isFwComplete(), false) ;
    s_fw.setSequenceFw(s) ; ASSERT_EQ(s_fw.isFwComplete(), true) ;

    // setting IPD then seq then PWD
    s_fw = ngsai::KineticSignal() ;
    ASSERT_EQ(s_fw.isFwComplete(), false) ;
    s_fw.setIPDFw(k)      ; ASSERT_EQ(s_fw.isFwComplete(), false) ;
    s_fw.setSequenceFw(s) ; ASSERT_EQ(s_fw.isFwComplete(), false) ;
    s_fw.setPWDFw(k)      ; ASSERT_EQ(s_fw.isFwComplete(), true) ;

    // setting PWD then IPD then seq
    s_fw = ngsai::KineticSignal() ;
    ASSERT_EQ(s_fw.isFwComplete(), false) ;
    s_fw.setPWDFw(k)      ; ASSERT_EQ(s_fw.isFwComplete(), false) ;
    s_fw.setIPDFw(k)      ; ASSERT_EQ(s_fw.isFwComplete(), false) ;
    s_fw.setSequenceFw(s) ; ASSERT_EQ(s_fw.isFwComplete(), true) ;

    // setting PWD then seq then IPD
    s_fw = ngsai::KineticSignal() ;
    ASSERT_EQ(s_fw.isFwComplete(), false) ;
    s_fw.setPWDFw(k)      ; ASSERT_EQ(s_fw.isFwComplete(), false) ;
    s_fw.setSequenceFw(s) ; ASSERT_EQ(s_fw.isFwComplete(), false) ;
    s_fw.setIPDFw(k)      ; ASSERT_EQ(s_fw.isFwComplete(), true) ;

    // setting seq then IPD then PWD
    s_fw = ngsai::KineticSignal() ;
    ASSERT_EQ(s_fw.isFwComplete(), false) ;
    s_fw.setSequenceFw(s) ; ASSERT_EQ(s_fw.isFwComplete(), false) ;
    s_fw.setIPDFw(k)      ; ASSERT_EQ(s_fw.isFwComplete(), false) ;
    s_fw.setPWDFw(k)      ; ASSERT_EQ(s_fw.isFwComplete(), true) ;

    // setting seq then PWD then IPD
    s_fw = ngsai::KineticSignal() ;
    ASSERT_EQ(s_fw.isFwComplete(), false) ;
    s_fw.setSequenceFw(s) ; ASSERT_EQ(s_fw.isFwComplete(), false) ;
    s_fw.setPWDFw(k)      ; ASSERT_EQ(s_fw.isFwComplete(), false) ;
    s_fw.setIPDFw(k)      ; ASSERT_EQ(s_fw.isFwComplete(), true) ;
}


/*!
 * \brief Checks that the isRvComplete method works properly.
 */
TEST(KineticSignalTest, isRvComplete)
{   
    std::string s("AAAAA") ;
    std::vector<double> k({1., 2., 3., 4., 5.}) ;
    
    ngsai::KineticSignal s_rv ;

    // setting IPD then PWD then seq
    s_rv = ngsai::KineticSignal() ;
    ASSERT_EQ(s_rv.isRvComplete(), false) ;
    s_rv.setIPDRv(k)      ; ASSERT_EQ(s_rv.isRvComplete(), false) ;
    s_rv.setPWDRv(k)      ; ASSERT_EQ(s_rv.isRvComplete(), false) ;
    s_rv.setSequenceRv(s) ; ASSERT_EQ(s_rv.isRvComplete(), true) ;

    // setting IPD then seq then PWD
    s_rv = ngsai::KineticSignal() ;
    ASSERT_EQ(s_rv.isRvComplete(), false) ;
    s_rv.setIPDRv(k)      ; ASSERT_EQ(s_rv.isRvComplete(), false) ;
    s_rv.setSequenceRv(s) ; ASSERT_EQ(s_rv.isRvComplete(), false) ;
    s_rv.setPWDRv(k)      ; ASSERT_EQ(s_rv.isRvComplete(), true) ;

    // setting PWD then IPD then seq
    s_rv = ngsai::KineticSignal() ;
    ASSERT_EQ(s_rv.isRvComplete(), false) ;
    s_rv.setPWDRv(k)      ; ASSERT_EQ(s_rv.isRvComplete(), false) ;
    s_rv.setIPDRv(k)      ; ASSERT_EQ(s_rv.isRvComplete(), false) ;
    s_rv.setSequenceRv(s) ; ASSERT_EQ(s_rv.isRvComplete(), true) ;

    // setting PWD then seq then IPD
    s_rv = ngsai::KineticSignal() ;
    ASSERT_EQ(s_rv.isRvComplete(), false) ;
    s_rv.setPWDRv(k)      ; ASSERT_EQ(s_rv.isRvComplete(), false) ;
    s_rv.setSequenceRv(s) ; ASSERT_EQ(s_rv.isRvComplete(), false) ;
    s_rv.setIPDRv(k)      ; ASSERT_EQ(s_rv.isRvComplete(), true) ;

    // setting seq then IPD then PWD
    s_rv = ngsai::KineticSignal() ;
    ASSERT_EQ(s_rv.isRvComplete(), false) ;
    s_rv.setSequenceRv(s) ; ASSERT_EQ(s_rv.isRvComplete(), false) ;
    s_rv.setIPDRv(k)      ; ASSERT_EQ(s_rv.isRvComplete(), false) ;
    s_rv.setPWDRv(k)      ; ASSERT_EQ(s_rv.isRvComplete(), true) ;

    // setting seq then PWD then IPD
    s_rv = ngsai::KineticSignal() ;
    ASSERT_EQ(s_rv.isRvComplete(), false) ;
    s_rv.setSequenceRv(s) ; ASSERT_EQ(s_rv.isRvComplete(), false) ;
    s_rv.setPWDRv(k)      ; ASSERT_EQ(s_rv.isRvComplete(), false) ;
    s_rv.setIPDRv(k)      ; ASSERT_EQ(s_rv.isRvComplete(), true) ;
}

/*!
 * \brief Checks that the isComplete method works properly.
 */
TEST(KineticSignalTest, isComplete)
{   
    std::string seq("AAAAA") ;
    std::vector<double> k({1., 2., 3., 4., 5.}) ;
    
    ngsai::KineticSignal s ;

    // setting IPDfw then PWDfw then seqfw then
    //         IPDrv then PWDrv then seqrv
    s = ngsai::KineticSignal() ;
    ASSERT_EQ(s.isComplete(), false) ;
    s.setIPDFw(k)        ; ASSERT_EQ(s.isComplete(), false) ;
    s.setPWDFw(k)        ; ASSERT_EQ(s.isComplete(), false) ;
    s.setSequenceFw(seq) ; ASSERT_EQ(s.isComplete(), false) ;
    s.setIPDRv(k)        ; ASSERT_EQ(s.isComplete(), false) ;
    s.setPWDRv(k)        ; ASSERT_EQ(s.isComplete(), false) ;
    s.setSequenceRv(seq) ; ASSERT_EQ(s.isComplete(), true) ;

    // setting IPDfw then PWDfw then seqfw then
    //         IPDrv then seqrv then PWDrv
    s = ngsai::KineticSignal() ;
    ASSERT_EQ(s.isComplete(), false) ;
    s.setIPDFw(k)        ; ASSERT_EQ(s.isComplete(), false) ;
    s.setPWDFw(k)        ; ASSERT_EQ(s.isComplete(), false) ;
    s.setSequenceFw(seq) ; ASSERT_EQ(s.isComplete(), false) ;
    s.setIPDRv(k)        ; ASSERT_EQ(s.isComplete(), false) ;
    s.setSequenceRv(seq) ; ASSERT_EQ(s.isComplete(), false) ;
    s.setPWDRv(k)        ; ASSERT_EQ(s.isComplete(), true) ;

    // setting IPDfw then PWDfw then seqfw then
    //         PWDrv then IPDrv then seqrv
    s = ngsai::KineticSignal() ;
    ASSERT_EQ(s.isComplete(), false) ;
    s.setIPDFw(k)        ; ASSERT_EQ(s.isComplete(), false) ;
    s.setPWDFw(k)        ; ASSERT_EQ(s.isComplete(), false) ;
    s.setSequenceFw(seq) ; ASSERT_EQ(s.isComplete(), false) ;
    s.setPWDRv(k)        ; ASSERT_EQ(s.isComplete(), false) ;
    s.setIPDRv(k)        ; ASSERT_EQ(s.isComplete(), false) ;
    s.setSequenceRv(seq) ; ASSERT_EQ(s.isComplete(), true) ;

    // setting IPDfw then PWDfw then seqfw then
    //         PWDrv then seqrv then IPDrv
    s = ngsai::KineticSignal() ;
    ASSERT_EQ(s.isComplete(), false) ;
    s.setIPDFw(k)        ; ASSERT_EQ(s.isComplete(), false) ;
    s.setPWDFw(k)        ; ASSERT_EQ(s.isComplete(), false) ;
    s.setSequenceFw(seq) ; ASSERT_EQ(s.isComplete(), false) ;
    s.setPWDRv(k)        ; ASSERT_EQ(s.isComplete(), false) ;
    s.setSequenceRv(seq) ; ASSERT_EQ(s.isComplete(), false) ;
    s.setIPDRv(k)        ; ASSERT_EQ(s.isComplete(), true) ;

    // setting IPDfw then PWDfw then seqfw then
    //         seqrv then IPDrv then PWDrv
    s = ngsai::KineticSignal() ;
    ASSERT_EQ(s.isComplete(), false) ;
    s.setIPDFw(k)        ; ASSERT_EQ(s.isComplete(), false) ;
    s.setPWDFw(k)        ; ASSERT_EQ(s.isComplete(), false) ;
    s.setSequenceFw(seq) ; ASSERT_EQ(s.isComplete(), false) ;
    s.setSequenceRv(seq) ; ASSERT_EQ(s.isComplete(), false) ;
    s.setIPDRv(k)        ; ASSERT_EQ(s.isComplete(), false) ;
    s.setPWDRv(k)        ; ASSERT_EQ(s.isComplete(), true) ;

    // setting IPDfw then PWDfw then seqfw then
    //         seqrv then PWDrv then IPDrv
    s = ngsai::KineticSignal() ;
    ASSERT_EQ(s.isComplete(), false) ;
    s.setIPDFw(k)        ; ASSERT_EQ(s.isComplete(), false) ;
    s.setPWDFw(k)        ; ASSERT_EQ(s.isComplete(), false) ;
    s.setSequenceFw(seq) ; ASSERT_EQ(s.isComplete(), false) ;
    s.setSequenceRv(seq) ; ASSERT_EQ(s.isComplete(), false) ;
    s.setIPDRv(k)        ; ASSERT_EQ(s.isComplete(), false) ;
    s.setPWDRv(k)        ; ASSERT_EQ(s.isComplete(), true) ;


    // ---------------------------------------------------------------

    // setting IPDfw then seqfw then PWDfw then
    //         IPDrv then PWDrv then seqrv
    s = ngsai::KineticSignal() ;
    ASSERT_EQ(s.isComplete(), false) ;
    s.setIPDFw(k)        ; ASSERT_EQ(s.isComplete(), false) ;
    s.setSequenceFw(seq) ; ASSERT_EQ(s.isComplete(), false) ;
    s.setPWDFw(k)        ; ASSERT_EQ(s.isComplete(), false) ;
    s.setIPDRv(k)        ; ASSERT_EQ(s.isComplete(), false) ;
    s.setPWDRv(k)        ; ASSERT_EQ(s.isComplete(), false) ;
    s.setSequenceRv(seq) ; ASSERT_EQ(s.isComplete(), true) ;


    // setting IPDfw then seqfw then PWDfw then
    //         IPDrv then seqrv then PWDrv
    s = ngsai::KineticSignal() ;
    ASSERT_EQ(s.isComplete(), false) ;
    s.setIPDFw(k)        ; ASSERT_EQ(s.isComplete(), false) ;
    s.setSequenceFw(seq) ; ASSERT_EQ(s.isComplete(), false) ;
    s.setPWDFw(k)        ; ASSERT_EQ(s.isComplete(), false) ;
    s.setIPDRv(k)        ; ASSERT_EQ(s.isComplete(), false) ;
    s.setSequenceRv(seq) ; ASSERT_EQ(s.isComplete(), false) ;
    s.setPWDRv(k)        ; ASSERT_EQ(s.isComplete(), true) ;

    // setting IPDfw then seqfw then PWDfw then
    //         PWDrv then IPDrv then seqrv
    s = ngsai::KineticSignal() ;
    ASSERT_EQ(s.isComplete(), false) ;
    s.setIPDFw(k)        ; ASSERT_EQ(s.isComplete(), false) ;
    s.setSequenceFw(seq) ; ASSERT_EQ(s.isComplete(), false) ;
    s.setPWDFw(k)        ; ASSERT_EQ(s.isComplete(), false) ;
    s.setPWDRv(k)        ; ASSERT_EQ(s.isComplete(), false) ;
    s.setIPDRv(k)        ; ASSERT_EQ(s.isComplete(), false) ;
    s.setSequenceRv(seq) ; ASSERT_EQ(s.isComplete(), true) ;

    // setting IPDfw then seqfw then PWDfw then
    //         PWDrv then seqrv then IPDrv
    s = ngsai::KineticSignal() ;
    ASSERT_EQ(s.isComplete(), false) ;
    s.setIPDFw(k)        ; ASSERT_EQ(s.isComplete(), false) ;
    s.setSequenceFw(seq) ; ASSERT_EQ(s.isComplete(), false) ;
    s.setPWDFw(k)        ; ASSERT_EQ(s.isComplete(), false) ;
    s.setPWDRv(k)        ; ASSERT_EQ(s.isComplete(), false) ;
    s.setSequenceRv(seq) ; ASSERT_EQ(s.isComplete(), false) ;
    s.setIPDRv(k)        ; ASSERT_EQ(s.isComplete(), true) ;

    // setting IPDfw then seqfw then PWDfw then
    //         seqrv then IPDrv then PWDrv
    s = ngsai::KineticSignal() ;
    ASSERT_EQ(s.isComplete(), false) ;
    s.setIPDFw(k)        ; ASSERT_EQ(s.isComplete(), false) ;
    s.setSequenceFw(seq) ; ASSERT_EQ(s.isComplete(), false) ;
    s.setPWDFw(k)        ; ASSERT_EQ(s.isComplete(), false) ;
    s.setSequenceRv(seq) ; ASSERT_EQ(s.isComplete(), false) ;
    s.setIPDRv(k)        ; ASSERT_EQ(s.isComplete(), false) ;
    s.setPWDRv(k)        ; ASSERT_EQ(s.isComplete(), true) ;

    // setting IPDfw then seqfw then PWDfw then
    //         seqrv then PWDrv then IPDrv
    s = ngsai::KineticSignal() ;
    ASSERT_EQ(s.isComplete(), false) ;
    s.setIPDFw(k)        ; ASSERT_EQ(s.isComplete(), false) ;
    s.setSequenceFw(seq) ; ASSERT_EQ(s.isComplete(), false) ;
    s.setPWDFw(k)        ; ASSERT_EQ(s.isComplete(), false) ;
    s.setSequenceRv(seq) ; ASSERT_EQ(s.isComplete(), false) ;
    s.setPWDRv(k)        ; ASSERT_EQ(s.isComplete(), false) ;
    s.setIPDRv(k)        ; ASSERT_EQ(s.isComplete(), true) ;

    // ---------------------------------------------------------------

    // setting PWDfw then IPDfw then seqfw then
    //         IPDrv then PWDrv then seqrv
    s = ngsai::KineticSignal() ;
    ASSERT_EQ(s.isComplete(), false) ;
    s.setPWDFw(k)        ; ASSERT_EQ(s.isComplete(), false) ;
    s.setIPDFw(k)        ; ASSERT_EQ(s.isComplete(), false) ;
    s.setSequenceFw(seq) ; ASSERT_EQ(s.isComplete(), false) ;
    s.setIPDRv(k)        ; ASSERT_EQ(s.isComplete(), false) ;
    s.setPWDRv(k)        ; ASSERT_EQ(s.isComplete(), false) ;
    s.setSequenceRv(seq) ; ASSERT_EQ(s.isComplete(), true) ;

    // setting PWDfw then IPDfw then seqfw then
    //         IPDrv then seqrv then PWDrv
    s = ngsai::KineticSignal() ;
    ASSERT_EQ(s.isComplete(), false) ;
    s.setPWDFw(k)        ; ASSERT_EQ(s.isComplete(), false) ;
    s.setIPDFw(k)        ; ASSERT_EQ(s.isComplete(), false) ;
    s.setSequenceFw(seq) ; ASSERT_EQ(s.isComplete(), false) ;
    s.setIPDRv(k)        ; ASSERT_EQ(s.isComplete(), false) ;
    s.setSequenceRv(seq) ; ASSERT_EQ(s.isComplete(), false) ;
    s.setPWDRv(k)        ; ASSERT_EQ(s.isComplete(), true) ;

    // setting PWDfw then IPDfw then seqfw then
    //         PWDrv then IPDrv then seqrv
    s = ngsai::KineticSignal() ;
    ASSERT_EQ(s.isComplete(), false) ;
    s.setPWDFw(k)        ; ASSERT_EQ(s.isComplete(), false) ;
    s.setIPDFw(k)        ; ASSERT_EQ(s.isComplete(), false) ;
    s.setSequenceFw(seq) ; ASSERT_EQ(s.isComplete(), false) ;
    s.setPWDRv(k)        ; ASSERT_EQ(s.isComplete(), false) ;
    s.setIPDRv(k)        ; ASSERT_EQ(s.isComplete(), false) ;
    s.setSequenceRv(seq) ; ASSERT_EQ(s.isComplete(), true) ;

    // setting PWDfw then IPDfw then seqfw then
    //         PWDrv then seqrv then PWDrv
    s = ngsai::KineticSignal() ;
    ASSERT_EQ(s.isComplete(), false) ;
    s.setPWDFw(k)        ; ASSERT_EQ(s.isComplete(), false) ;
    s.setIPDFw(k)        ; ASSERT_EQ(s.isComplete(), false) ;
    s.setSequenceFw(seq) ; ASSERT_EQ(s.isComplete(), false) ;
    s.setPWDRv(k)        ; ASSERT_EQ(s.isComplete(), false) ;
    s.setSequenceRv(seq) ; ASSERT_EQ(s.isComplete(), false) ;
    s.setIPDRv(k)        ; ASSERT_EQ(s.isComplete(), true) ;

    // setting PWDfw then IPDfw then seqfw then
    //         seqrv then IPDrv then PWDrv
    s = ngsai::KineticSignal() ;
    ASSERT_EQ(s.isComplete(), false) ;
    s.setPWDFw(k)        ; ASSERT_EQ(s.isComplete(), false) ;
    s.setIPDFw(k)        ; ASSERT_EQ(s.isComplete(), false) ;
    s.setSequenceFw(seq) ; ASSERT_EQ(s.isComplete(), false) ;
    s.setSequenceRv(seq) ; ASSERT_EQ(s.isComplete(), false) ;
    s.setIPDRv(k)        ; ASSERT_EQ(s.isComplete(), false) ;
    s.setPWDRv(k)        ; ASSERT_EQ(s.isComplete(), true) ;

    // setting PWDfw then IPDfw then seqfw then
    //         seqrv then PWDrv then seqrv
    s = ngsai::KineticSignal() ;
    ASSERT_EQ(s.isComplete(), false) ;
    s.setPWDFw(k)        ; ASSERT_EQ(s.isComplete(), false) ;
    s.setIPDFw(k)        ; ASSERT_EQ(s.isComplete(), false) ;
    s.setSequenceFw(seq) ; ASSERT_EQ(s.isComplete(), false) ;
    s.setSequenceRv(seq) ; ASSERT_EQ(s.isComplete(), false) ;
    s.setPWDRv(k)        ; ASSERT_EQ(s.isComplete(), false) ;
    s.setIPDRv(k)        ; ASSERT_EQ(s.isComplete(), true) ;


    // ---------------------------------------------------------------

    // setting PWDfw then seqfw then IPDfw then
    //         IPDrv then PWDrv then seqrv
    s = ngsai::KineticSignal() ;
    ASSERT_EQ(s.isComplete(), false) ;
    s.setPWDFw(k)        ; ASSERT_EQ(s.isComplete(), false) ;
    s.setSequenceFw(seq) ; ASSERT_EQ(s.isComplete(), false) ;
    s.setIPDFw(k)        ; ASSERT_EQ(s.isComplete(), false) ;
    s.setIPDRv(k)        ; ASSERT_EQ(s.isComplete(), false) ;
    s.setPWDRv(k)        ; ASSERT_EQ(s.isComplete(), false) ;
    s.setSequenceRv(seq) ; ASSERT_EQ(s.isComplete(), true) ;

    // setting PWDfw then seqfw then IPDfw then
    //         IPDrv then seqrv then PWDrv
    s = ngsai::KineticSignal() ;
    ASSERT_EQ(s.isComplete(), false) ;
    s.setPWDFw(k)        ; ASSERT_EQ(s.isComplete(), false) ;
    s.setSequenceFw(seq) ; ASSERT_EQ(s.isComplete(), false) ;
    s.setIPDFw(k)        ; ASSERT_EQ(s.isComplete(), false) ;
    s.setIPDRv(k)        ; ASSERT_EQ(s.isComplete(), false) ;
    s.setSequenceRv(seq) ; ASSERT_EQ(s.isComplete(), false) ;
    s.setPWDRv(k)        ; ASSERT_EQ(s.isComplete(), true) ;

    // setting PWDfw then seqfw then IPDfw then
    //         PWDrv then IPDrv then seqrv
    s = ngsai::KineticSignal() ;
    ASSERT_EQ(s.isComplete(), false) ;
    s.setPWDFw(k)        ; ASSERT_EQ(s.isComplete(), false) ;
    s.setSequenceFw(seq) ; ASSERT_EQ(s.isComplete(), false) ;
    s.setIPDFw(k)        ; ASSERT_EQ(s.isComplete(), false) ;
    s.setPWDRv(k)        ; ASSERT_EQ(s.isComplete(), false) ;
    s.setIPDRv(k)        ; ASSERT_EQ(s.isComplete(), false) ;
    s.setSequenceRv(seq) ; ASSERT_EQ(s.isComplete(), true) ;

    // setting PWDfw then seqfw then IPDfw then
    //         PWDrv then seqrv then IPDrv
    s = ngsai::KineticSignal() ;
    ASSERT_EQ(s.isComplete(), false) ;
    s.setPWDFw(k)        ; ASSERT_EQ(s.isComplete(), false) ;
    s.setSequenceFw(seq) ; ASSERT_EQ(s.isComplete(), false) ;
    s.setIPDFw(k)        ; ASSERT_EQ(s.isComplete(), false) ;
    s.setPWDRv(k)        ; ASSERT_EQ(s.isComplete(), false) ;
    s.setSequenceRv(seq) ; ASSERT_EQ(s.isComplete(), false) ;
    s.setIPDRv(k)        ; ASSERT_EQ(s.isComplete(), true) ;

    // setting PWDfw then seqfw then IPDfw then
    //         seqrv then IPDrv then PWDrv
    s = ngsai::KineticSignal() ;
    ASSERT_EQ(s.isComplete(), false) ;
    s.setPWDFw(k)        ; ASSERT_EQ(s.isComplete(), false) ;
    s.setSequenceFw(seq) ; ASSERT_EQ(s.isComplete(), false) ;
    s.setIPDFw(k)        ; ASSERT_EQ(s.isComplete(), false) ;
    s.setSequenceRv(seq) ; ASSERT_EQ(s.isComplete(), false) ;
    s.setIPDRv(k)        ; ASSERT_EQ(s.isComplete(), false) ;
    s.setPWDRv(k)        ; ASSERT_EQ(s.isComplete(), true) ;

    // setting PWDfw then seqfw then IPDfw then
    //         seqrv then PWDrv then IPDrv
    s = ngsai::KineticSignal() ;
    ASSERT_EQ(s.isComplete(), false) ;
    s.setPWDFw(k)        ; ASSERT_EQ(s.isComplete(), false) ;
    s.setSequenceFw(seq) ; ASSERT_EQ(s.isComplete(), false) ;
    s.setIPDFw(k)        ; ASSERT_EQ(s.isComplete(), false) ;
    s.setSequenceRv(seq) ; ASSERT_EQ(s.isComplete(), false) ;
    s.setPWDRv(k)        ; ASSERT_EQ(s.isComplete(), false) ;
    s.setIPDRv(k)        ; ASSERT_EQ(s.isComplete(), true) ;


    // ---------------------------------------------------------------

    // setting seqfw then IPDfw then PWDfw then
    //         IPDrv then PWDrv then seqrv
    s = ngsai::KineticSignal() ;
    ASSERT_EQ(s.isComplete(), false) ;
    s.setSequenceFw(seq) ; ASSERT_EQ(s.isComplete(), false) ;
    s.setIPDFw(k)        ; ASSERT_EQ(s.isComplete(), false) ;
    s.setPWDFw(k)        ; ASSERT_EQ(s.isComplete(), false) ;
    s.setIPDRv(k)        ; ASSERT_EQ(s.isComplete(), false) ;
    s.setPWDRv(k)        ; ASSERT_EQ(s.isComplete(), false) ;
    s.setSequenceRv(seq) ; ASSERT_EQ(s.isComplete(), true) ;

    // setting seqfw then IPDfw then PWDfw then
    //         IPDrv then seqrv then PWDrv
    s = ngsai::KineticSignal() ;
    ASSERT_EQ(s.isComplete(), false) ;
    s.setSequenceFw(seq) ; ASSERT_EQ(s.isComplete(), false) ;
    s.setIPDFw(k)        ; ASSERT_EQ(s.isComplete(), false) ;
    s.setPWDFw(k)        ; ASSERT_EQ(s.isComplete(), false) ;
    s.setIPDRv(k)        ; ASSERT_EQ(s.isComplete(), false) ;
    s.setSequenceRv(seq) ; ASSERT_EQ(s.isComplete(), false) ;
    s.setPWDRv(k)        ; ASSERT_EQ(s.isComplete(), true) ;

    // setting seqfw then IPDfw then PWDfw then
    //         PWDrv then IPDrv then seqrv
    s = ngsai::KineticSignal() ;
    ASSERT_EQ(s.isComplete(), false) ;
    s.setSequenceFw(seq) ; ASSERT_EQ(s.isComplete(), false) ;
    s.setIPDFw(k)        ; ASSERT_EQ(s.isComplete(), false) ;
    s.setPWDFw(k)        ; ASSERT_EQ(s.isComplete(), false) ;
    s.setPWDRv(k)        ; ASSERT_EQ(s.isComplete(), false) ;
    s.setIPDRv(k)        ; ASSERT_EQ(s.isComplete(), false) ;
    s.setSequenceRv(seq) ; ASSERT_EQ(s.isComplete(), true) ;

    // setting seqfw then IPDfw then PWDfw then
    //         PWDrv then seqrv then IPDrv
    s = ngsai::KineticSignal() ;
    ASSERT_EQ(s.isComplete(), false) ;
    s.setSequenceFw(seq) ; ASSERT_EQ(s.isComplete(), false) ;
    s.setIPDFw(k)        ; ASSERT_EQ(s.isComplete(), false) ;
    s.setPWDFw(k)        ; ASSERT_EQ(s.isComplete(), false) ;
    s.setPWDRv(k)        ; ASSERT_EQ(s.isComplete(), false) ;
    s.setSequenceRv(seq) ; ASSERT_EQ(s.isComplete(), false) ;
    s.setIPDRv(k)        ; ASSERT_EQ(s.isComplete(), true) ;

    // setting seqfw then IPDfw then PWDfw then
    //         seqrv then IPDrv then PWDrv
    s = ngsai::KineticSignal() ;
    ASSERT_EQ(s.isComplete(), false) ;
    s.setSequenceFw(seq) ; ASSERT_EQ(s.isComplete(), false) ;
    s.setIPDFw(k)        ; ASSERT_EQ(s.isComplete(), false) ;
    s.setPWDFw(k)        ; ASSERT_EQ(s.isComplete(), false) ;
    s.setSequenceRv(seq) ; ASSERT_EQ(s.isComplete(), false) ;
    s.setIPDRv(k)        ; ASSERT_EQ(s.isComplete(), false) ;
    s.setPWDRv(k)        ; ASSERT_EQ(s.isComplete(), true) ;

    // setting seqfw then IPDfw then PWDfw then
    //         seqrv then PWDrv then IPDrv
    s = ngsai::KineticSignal() ;
    ASSERT_EQ(s.isComplete(), false) ;
    s.setSequenceFw(seq) ; ASSERT_EQ(s.isComplete(), false) ;
    s.setIPDFw(k)        ; ASSERT_EQ(s.isComplete(), false) ;
    s.setPWDFw(k)        ; ASSERT_EQ(s.isComplete(), false) ;
    s.setSequenceRv(seq) ; ASSERT_EQ(s.isComplete(), false) ;
    s.setPWDRv(k)        ; ASSERT_EQ(s.isComplete(), false) ;
    s.setIPDRv(k)        ; ASSERT_EQ(s.isComplete(), true) ;


    // ---------------------------------------------------------------

    // setting seqfw then PWDfw then IPDfw then
    //         IPDrv then PWDrv then seqrv
    s = ngsai::KineticSignal() ;
    ASSERT_EQ(s.isComplete(), false) ;
    s.setSequenceFw(seq) ; ASSERT_EQ(s.isComplete(), false) ;
    s.setPWDFw(k)        ; ASSERT_EQ(s.isComplete(), false) ;
    s.setIPDFw(k)        ; ASSERT_EQ(s.isComplete(), false) ;
    s.setIPDRv(k)        ; ASSERT_EQ(s.isComplete(), false) ;
    s.setPWDRv(k)        ; ASSERT_EQ(s.isComplete(), false) ;
    s.setSequenceRv(seq) ; ASSERT_EQ(s.isComplete(), true) ;

    // setting seqfw then PWDfw then IPDfw then
    //         IPDrv then seqrv then PWDrv
    s = ngsai::KineticSignal() ;
    ASSERT_EQ(s.isComplete(), false) ;
    s.setSequenceFw(seq) ; ASSERT_EQ(s.isComplete(), false) ;
    s.setPWDFw(k)        ; ASSERT_EQ(s.isComplete(), false) ;
    s.setIPDFw(k)        ; ASSERT_EQ(s.isComplete(), false) ;
    s.setIPDRv(k)        ; ASSERT_EQ(s.isComplete(), false) ;
    s.setSequenceRv(seq) ; ASSERT_EQ(s.isComplete(), false) ;
    s.setPWDRv(k)        ; ASSERT_EQ(s.isComplete(), true) ;

    // setting seqfw then PWDfw then IPDfw then
    //         PWDrv then IPDrv then seqrv
    s = ngsai::KineticSignal() ;
    ASSERT_EQ(s.isComplete(), false) ;
    s.setSequenceFw(seq) ; ASSERT_EQ(s.isComplete(), false) ;
    s.setPWDFw(k)        ; ASSERT_EQ(s.isComplete(), false) ;
    s.setIPDFw(k)        ; ASSERT_EQ(s.isComplete(), false) ;
    s.setPWDRv(k)        ; ASSERT_EQ(s.isComplete(), false) ;
    s.setIPDRv(k)        ; ASSERT_EQ(s.isComplete(), false) ;
    s.setSequenceRv(seq) ; ASSERT_EQ(s.isComplete(), true) ;

    // setting seqfw then PWDfw then IPDfw then
    //         PWDrv then seqrv then IPDrv
    s = ngsai::KineticSignal() ;
    ASSERT_EQ(s.isComplete(), false) ;
    s.setSequenceFw(seq) ; ASSERT_EQ(s.isComplete(), false) ;
    s.setPWDFw(k)        ; ASSERT_EQ(s.isComplete(), false) ;
    s.setIPDFw(k)        ; ASSERT_EQ(s.isComplete(), false) ;
    s.setPWDRv(k)        ; ASSERT_EQ(s.isComplete(), false) ;
    s.setSequenceRv(seq) ; ASSERT_EQ(s.isComplete(), false) ;
    s.setIPDRv(k)        ; ASSERT_EQ(s.isComplete(), true) ;

    // setting seqfw then PWDfw then IPDfw then
    //         seqrv then IPDrv then PWDrv
    s = ngsai::KineticSignal() ;
    ASSERT_EQ(s.isComplete(), false) ;
    s.setSequenceFw(seq) ; ASSERT_EQ(s.isComplete(), false) ;
    s.setPWDFw(k)        ; ASSERT_EQ(s.isComplete(), false) ;
    s.setIPDFw(k)        ; ASSERT_EQ(s.isComplete(), false) ;
    s.setSequenceRv(seq) ; ASSERT_EQ(s.isComplete(), false) ;
    s.setIPDRv(k)        ; ASSERT_EQ(s.isComplete(), false) ;
    s.setPWDRv(k)        ; ASSERT_EQ(s.isComplete(), true) ;

    // setting seqfw then PWDfw then IPDfw then
    //         seqrv then PWDrv then IPDrv
    s = ngsai::KineticSignal() ;
    ASSERT_EQ(s.isComplete(), false) ;
    s.setSequenceFw(seq) ; ASSERT_EQ(s.isComplete(), false) ;
    s.setPWDFw(k)        ; ASSERT_EQ(s.isComplete(), false) ;
    s.setIPDFw(k)        ; ASSERT_EQ(s.isComplete(), false) ;
    s.setSequenceRv(seq) ; ASSERT_EQ(s.isComplete(), false) ;
    s.setPWDRv(k)        ; ASSERT_EQ(s.isComplete(), false) ;
    s.setIPDRv(k)        ; ASSERT_EQ(s.isComplete(), true) ;

}

