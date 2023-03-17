#include <gtest/gtest.h>
#include <gmock/gmock.h>
#include <string>
#include <stdexcept>
#include <utility>        // std::move()

#include <ngsaipp/genome/CpGRegion.hpp>
#include <ngsaipp/genome/GenomeRegion.hpp>
#include <ngsaipp/genome/constants.hpp>


/*
 * This file contains tests for the ngsai::genome::CpGRegion class from 
 * src/genome/CpGRegion.cpp
 */


/*!
 * \brief Constructs a CpGRegion from a GenomeRegion.
 * \param region the GenomeRegion of interest.
 * \returns the CpGRegion.
 */
ngsai::genome::CpGRegion
construct(const ngsai::genome::GenomeRegion& region)
{   return ngsai::genome::CpGRegion(region) ; }


/*!
 * \brief Allocates a CpGRegion as a GenomeRegion.
 * \param chrom the chromosome.
 * \param start the start coordinate.
 * \param end the end coordinate.
 * \returns the pointer to the allocated CpGRegion.
 */
ngsai::genome::GenomeRegion* 
allocate(const std::string& chrom,
         size_t start,
         size_t end)
{   ngsai::genome::GenomeRegion* region =
                        new ngsai::genome::CpGRegion(chrom, start, end) ;
    return region ;
}


/*!
 * \brief Checks that the default constructor works properly.
 */
TEST(CpGRegionTest, constructor_default)
{
    ngsai::genome::CpGRegion region ;
    ASSERT_EQ(region.chrom, std::string("")) ;
    ASSERT_EQ(region.start, 0) ;
    ASSERT_EQ(region.end,   0) ;
    ASSERT_EQ(region.strand, ngsai::genome::strand::UNORIENTED) ;
}


/*!
 * \brief Checks that the parameter constructor works properly.
 */
TEST(CpGRegionTest, constructor_param)
{
    ngsai::genome::strand un = ngsai::genome::strand::UNORIENTED ;

    // empty chrom
    ngsai::genome::CpGRegion cpg1("", 1000, 1002) ;
    ASSERT_EQ(cpg1.chrom,  std::string("")) ;
    ASSERT_EQ(cpg1.start,  1000) ;
    ASSERT_EQ(cpg1.end,    1002) ;
    ASSERT_EQ(cpg1.strand, un) ;

    // regular un
    ngsai::genome::CpGRegion cpg2("chr1", 1000, 1002) ;
    ASSERT_EQ(cpg2.chrom,  std::string("chr1")) ;
    ASSERT_EQ(cpg2.start,  1000) ;
    ASSERT_EQ(cpg2.end,    1002) ;
    ASSERT_EQ(cpg2.strand, un) ;

    // should not work because start > end
    EXPECT_THROW(ngsai::genome::CpGRegion("chr1", 1001, 1000),
                 std::invalid_argument) ;
     // should not work because start == end
    EXPECT_THROW(ngsai::genome::CpGRegion("chr1", 1000, 1000),
                 std::invalid_argument) ;
    // should not work because length < 2
    EXPECT_THROW(ngsai::genome::CpGRegion("chr1", 1000, 1001),
                 std::invalid_argument) ;
    // should not work because length > 2
    EXPECT_THROW(ngsai::genome::CpGRegion("chr1", 1000, 1003),
                 std::invalid_argument) ;
}


/*!
 * \brief Checks that the copy constructor works properly.
 */
TEST(CpGRegionTest, constructor_copy)
{   
    ngsai::genome::strand fw = ngsai::genome::strand::FORWARD ;
    ngsai::genome::strand rv = ngsai::genome::strand::REVERSE ;
    ngsai::genome::strand un = ngsai::genome::strand::UNORIENTED ;

    // copying from CpGRegion
    ngsai::genome::CpGRegion cpg1("chr1", 1000, 1002) ;
    ngsai::genome::CpGRegion cpg2(cpg1) ;

    ASSERT_EQ(cpg1.chrom,  cpg2.chrom) ;
    ASSERT_EQ(cpg1.start,  cpg2.start) ;
    ASSERT_EQ(cpg1.end,    cpg2.end) ;
    ASSERT_EQ(cpg1.strand, cpg2.strand) ;

    // copying from GenomeRegion
    ngsai::genome::GenomeRegion region1("chr1", 1000, 1002, un) ;
    ngsai::genome::CpGRegion cpg3(region1) ;
    ASSERT_EQ(cpg3.chrom,  region1.chrom) ;
    ASSERT_EQ(cpg3.start,  region1.start) ;
    ASSERT_EQ(cpg3.end,    region1.end) ;
    ASSERT_EQ(cpg3.strand, region1.strand) ;

    // must fail because strand is not UNORIENTED
    ngsai::genome::GenomeRegion region2("chr1", 1000, 1002, fw) ;
    ngsai::genome::GenomeRegion region3("chr1", 1000, 1002, rv) ;
    EXPECT_THROW(construct(region2), std::invalid_argument) ;
    EXPECT_THROW(construct(region3), std::invalid_argument) ;

    // must fail because length < 2
    ngsai::genome::GenomeRegion region4("chr1", 1000, 1001, un) ;
    EXPECT_THROW(construct(region4), std::invalid_argument) ;

    // must fail because length > 2
    ngsai::genome::GenomeRegion region5("chr1", 1000, 1003, un) ;
    EXPECT_THROW(construct(region5), std::invalid_argument) ;
}


/*!
 * \brief Checks that the move constructor works properly.
 */
TEST(CpGRegionTest, constructor_move)
{
    ngsai::genome::CpGRegion cpg1("chr1", 1000, 1002) ;
    ngsai::genome::CpGRegion cpg2(cpg1) ;
    ngsai::genome::CpGRegion cpg3(std::move(cpg2)) ;

    ASSERT_EQ(cpg1.chrom,  cpg3.chrom) ;
    ASSERT_EQ(cpg1.start,  cpg3.start) ;
    ASSERT_EQ(cpg1.end,    cpg3.end) ;
    ASSERT_EQ(cpg1.strand, cpg3.strand) ;
}


/*!
 * \brief Checks that the assignment operator works properly.
 */
TEST(CpGRegionTest, assignment_operator)
{
    ngsai::genome::CpGRegion cpg1("chr1",
                                  1000,
                                  1002) ;
    ngsai::genome::CpGRegion cpg2 = cpg1 ;

    ASSERT_EQ(cpg1.chrom,  cpg2.chrom) ;
    ASSERT_EQ(cpg1.start,  cpg2.start) ;
    ASSERT_EQ(cpg1.end,    cpg2.end) ;
    ASSERT_EQ(cpg1.strand, cpg2.strand) ;
}


/*!
 * \brief Checks that the forwardCoordinates() method works properly.
 */
TEST(CpGRegionTest, forwardCoordinates)
{   
    ngsai::genome::CpGRegion cpg1("chr1", 1000, 1002) ;
    ngsai::genome::GenomeRegion cpg_fw("chr1",
                                       1000,
                                       1002,
                                       ngsai::genome::strand::FORWARD) ;
    
    ASSERT_EQ(cpg1.forwardCoordinates(), cpg_fw) ;
}


/*!
 * \brief Checks that the reverseCoordinates() method works properly.
 */
TEST(CpGRegionTest, reverseCoordinates)
{   
    ngsai::genome::CpGRegion cpg1("chr1", 1000, 1002) ;
    ngsai::genome::GenomeRegion cpg_rv("chr1",
                                       1000,
                                       1002,
                                       ngsai::genome::strand::REVERSE) ;
    
    ASSERT_EQ(cpg1.reverseCoordinates(), cpg_rv) ;
}


/*!
 * \brief Checks that polymorphism - allocation of CpGRegion as GenomeRegion - 
 * works properly.
 */
TEST(CpGRegionTest, polymorphism)
{   
    ngsai::genome::strand un = ngsai::genome::strand::UNORIENTED ;

    ngsai::genome::GenomeRegion* cpg = allocate("chr1", 1000, 1002) ;
    ASSERT_EQ(cpg->chrom,  std::string("chr1")) ;
    ASSERT_EQ(cpg->start,  1000) ;
    ASSERT_EQ(cpg->end,    1002) ;
    ASSERT_EQ(cpg->strand, un) ;
    delete cpg ;
    cpg = nullptr ;

    // this must fail because length < 2
    EXPECT_THROW(allocate("chr1", 1000, 1001), std::invalid_argument) ;

    // this must fail because length > 2
    EXPECT_THROW(allocate("chr1", 1000, 1003), std::invalid_argument) ;
}