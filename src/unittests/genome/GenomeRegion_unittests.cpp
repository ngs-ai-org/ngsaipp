#include <gtest/gtest.h>
#include <gmock/gmock.h>
#include <string>
#include <stdexcept>
#include <utility>        // std::move()

#include <ngsaipp/genome/GenomeRegion.hpp>
#include <ngsaipp/genome/constants.hpp>


/*
 * This file contains tests for the ngsai::genome::GenomeRegion class from 
 * src/genome/GenomeRegion.cpp
 */


/*!
 * \brief Checks that the default constructor works properly.
 */
TEST(GenomeRegionTest, constructor_default)
{
    ngsai::genome::GenomeRegion region ;
    ASSERT_EQ(region.chrom, std::string("")) ;
    ASSERT_EQ(region.start, 0) ;
    ASSERT_EQ(region.end,   0) ;
    ASSERT_EQ(region.strand, ngsai::genome::strand::UNORIENTED) ;
}


/*!
 * \brief Checks that the parameter constructor works properly.
 */
TEST(GenomeRegionTest, constructor_param)
{
    ngsai::genome::strand fw = ngsai::genome::strand::FORWARD ;
    ngsai::genome::strand rv = ngsai::genome::strand::REVERSE ;
    ngsai::genome::strand un = ngsai::genome::strand::UNORIENTED ;

    // empty chrom
    ngsai::genome::GenomeRegion region1("", 1000, 1500, fw) ;
    ASSERT_EQ(region1.chrom,  std::string("")) ;
    ASSERT_EQ(region1.start,  1000) ;
    ASSERT_EQ(region1.end,    1500) ;
    ASSERT_EQ(region1.strand, fw) ;

    // regular fw
    ngsai::genome::GenomeRegion region2("chr1", 1000, 1500, fw) ;
    ASSERT_EQ(region2.chrom,  std::string("chr1")) ;
    ASSERT_EQ(region2.start,  1000) ;
    ASSERT_EQ(region2.end,    1500) ;
    ASSERT_EQ(region2.strand, fw) ;

    // regular rv
    ngsai::genome::GenomeRegion region3("chr1", 1000, 1500, rv) ;
    ASSERT_EQ(region3.chrom,  std::string("chr1")) ;
    ASSERT_EQ(region3.start,  1000) ;
    ASSERT_EQ(region3.end,    1500) ;
    ASSERT_EQ(region3.strand, rv) ;

    // regular un
    ngsai::genome::GenomeRegion region4("chr1", 1000, 1500, un) ;
    ASSERT_EQ(region4.chrom,  std::string("chr1")) ;
    ASSERT_EQ(region4.start,  1000) ;
    ASSERT_EQ(region4.end,    1500) ;
    ASSERT_EQ(region4.strand, un) ;

    // single position fw
    ngsai::genome::GenomeRegion region5("chr1", 1000, 1001, fw) ;
    ASSERT_EQ(region5.chrom,  std::string("chr1")) ;
    ASSERT_EQ(region5.start,  1000) ;
    ASSERT_EQ(region5.end,    1001) ;
    ASSERT_EQ(region5.strand, fw) ;

    // single position rv
    ngsai::genome::GenomeRegion region6("chr1", 1000, 1001, rv) ;
    ASSERT_EQ(region6.chrom,  std::string("chr1")) ;
    ASSERT_EQ(region6.start,  1000) ;
    ASSERT_EQ(region6.end,    1001) ;
    ASSERT_EQ(region6.strand, rv) ;

    // single position un
    ngsai::genome::GenomeRegion region7("chr1", 1000, 1001, un) ;
    ASSERT_EQ(region7.chrom,  std::string("chr1")) ;
    ASSERT_EQ(region7.start,  1000) ;
    ASSERT_EQ(region7.end,    1001) ;
    ASSERT_EQ(region7.strand, un) ;

    // should fail because start == end
    EXPECT_THROW(ngsai::genome::GenomeRegion("chr1", 1000, 1000, fw),
                 std::invalid_argument) ;
    EXPECT_THROW(ngsai::genome::GenomeRegion("chr1", 1000, 1000, rv),
                 std::invalid_argument) ;
    EXPECT_THROW(ngsai::genome::GenomeRegion("chr1", 1000, 1000, un),
                 std::invalid_argument) ;

    // should fail because start > end
    EXPECT_THROW(ngsai::genome::GenomeRegion("chr1", 1001, 1000, fw),
                 std::invalid_argument) ;
    EXPECT_THROW(ngsai::genome::GenomeRegion("chr1", 1001, 1000, rv),
                 std::invalid_argument) ;
    EXPECT_THROW(ngsai::genome::GenomeRegion("chr1", 1001, 1000, un),
                 std::invalid_argument) ;
}


/*!
 * \brief Checks that the copy constructor works properly.
 */
TEST(GenomeRegionTest, constructor_copy)
{
    ngsai::genome::GenomeRegion region1("chr1",
                                        1000,
                                        1500,
                                        ngsai::genome::strand::FORWARD) ;
    ngsai::genome::GenomeRegion region2(region1) ;

    ASSERT_EQ(region1.chrom,  region2.chrom) ;
    ASSERT_EQ(region1.start,  region2.start) ;
    ASSERT_EQ(region1.end,    region2.end) ;
    ASSERT_EQ(region1.strand, region2.strand) ;
}


/*!
 * \brief Checks that the move constructor works properly.
 */
TEST(GenomeRegionTest, constructor_move)
{
    ngsai::genome::GenomeRegion region1("chr1",
                                        1000,
                                        1500,
                                        ngsai::genome::strand::FORWARD) ;
    ngsai::genome::GenomeRegion region2(region1) ;
    ngsai::genome::GenomeRegion region3(std::move(region2)) ;

    ASSERT_EQ(region1.chrom,  region3.chrom) ;
    ASSERT_EQ(region1.start,  region3.start) ;
    ASSERT_EQ(region1.end,    region3.end) ;
    ASSERT_EQ(region1.strand, region3.strand) ;
}


/*!
 * \brief Checks that the assignment operator works properly.
 */
TEST(GenomeRegionTest, assignment_operator)
{
    ngsai::genome::GenomeRegion region1("chr1",
                                        1000,
                                        1500,
                                        ngsai::genome::strand::FORWARD) ;
    ngsai::genome::GenomeRegion region2 = region1 ;

    ASSERT_EQ(region1.chrom,  region2.chrom) ;
    ASSERT_EQ(region1.start,  region2.start) ;
    ASSERT_EQ(region1.end,    region2.end) ;
    ASSERT_EQ(region1.strand, region2.strand) ;
}


/*!
 * \brief Checks that the move assignment operator works properly.
 */
TEST(GenomeRegionTest, assignment_move_operator)
{
    ngsai::genome::GenomeRegion region1("chr1",
                                        1000,
                                        1500,
                                        ngsai::genome::strand::FORWARD) ;
    ngsai::genome::GenomeRegion region2(region1) ;
    ngsai::genome::GenomeRegion region3 = std::move(region2) ;

    ASSERT_EQ(region1.chrom,  region3.chrom) ;
    ASSERT_EQ(region1.start,  region3.start) ;
    ASSERT_EQ(region1.end,    region3.end) ;
    ASSERT_EQ(region1.strand, region3.strand) ;
}


/*!
 * \brief Checks that the equality comparison operator works properly.
 */
TEST(GenomeRegionTest, equality_operator)
{   ngsai::genome::GenomeRegion region1("chr1",
                                        1000,
                                        1500,
                                        ngsai::genome::strand::FORWARD) ;
    ngsai::genome::GenomeRegion region2("chr1",
                                        1000,
                                        1500,
                                        ngsai::genome::strand::FORWARD) ;
    ngsai::genome::GenomeRegion region3("chr2",
                                        1000,
                                        1500,
                                        ngsai::genome::strand::FORWARD) ;
    ngsai::genome::GenomeRegion region4("chr1",
                                        1001,
                                        1500,
                                        ngsai::genome::strand::FORWARD) ;
    ngsai::genome::GenomeRegion region5("chr1",
                                        1000,
                                        1501,
                                        ngsai::genome::strand::FORWARD) ;
    ngsai::genome::GenomeRegion region6("chr1",
                                        1000,
                                        1500,
                                        ngsai::genome::strand::REVERSE) ;

    ASSERT_EQ(region1 == region2, true) ;
    ASSERT_EQ(region1 == region3, false) ;
    ASSERT_EQ(region1 == region4, false) ;
    ASSERT_EQ(region1 == region5, false) ;
    ASSERT_EQ(region1 == region6, false) ;
}


/*!
 * \brief Checks that the inequality comparison operator works properly.
 */
TEST(GenomeRegionTest, inequality_operator)
{   ngsai::genome::GenomeRegion region1("chr1",
                                        1000,
                                        1500,
                                        ngsai::genome::strand::FORWARD) ;
    ngsai::genome::GenomeRegion region2("chr1",
                                        1000,
                                        1500,
                                        ngsai::genome::strand::FORWARD) ;
    ngsai::genome::GenomeRegion region3("chr2",
                                        1000,
                                        1500,
                                        ngsai::genome::strand::FORWARD) ;
    ngsai::genome::GenomeRegion region4("chr1",
                                        1001,
                                        1500,
                                        ngsai::genome::strand::FORWARD) ;
    ngsai::genome::GenomeRegion region5("chr1",
                                        1000,
                                        1501,
                                        ngsai::genome::strand::FORWARD) ;
    ngsai::genome::GenomeRegion region6("chr1",
                                        1000,
                                        1500,
                                        ngsai::genome::strand::REVERSE) ;

    ASSERT_EQ(region1 != region2, false) ;
    ASSERT_EQ(region1 != region3, true) ;
    ASSERT_EQ(region1 != region4, true) ;
    ASSERT_EQ(region1 != region5, true) ;
    ASSERT_EQ(region1 != region6, true) ;
}


/*!
 * \brief Checks that the overlap operator works properly.
 */
TEST(GenomeRegionTest, overlap_operator)
{
    std::vector<ngsai::genome::strand> strands =
                    {ngsai::genome::strand::FORWARD,
                     ngsai::genome::strand::REVERSE,
                     ngsai::genome::strand::UNORIENTED} ;

    for(auto std1 : strands)
    {   for(auto std2 : strands)
        {
            // 2nd is upstream
            ASSERT_EQ(ngsai::genome::GenomeRegion("chr1", 100, 200, std1) | 
                      ngsai::genome::GenomeRegion("chr1",  50,  60, std2),
                      false) ;
            ASSERT_EQ(ngsai::genome::GenomeRegion("chr1", 100, 200, std1) | 
                      ngsai::genome::GenomeRegion("chr1",  90, 100, std2),
                      false) ;
            ASSERT_EQ(ngsai::genome::GenomeRegion("chr1", 100, 200, std1) | 
                      ngsai::genome::GenomeRegion("chr1",  99, 100, std2),
                      false) ;

            // 2nd overlaps start
            ASSERT_EQ(ngsai::genome::GenomeRegion("chr1", 100, 200, std1) | 
                      ngsai::genome::GenomeRegion("chr1",  90, 101, std2),
                      true) ;
            ASSERT_EQ(ngsai::genome::GenomeRegion("chr1", 100, 200, std1) | 
                      ngsai::genome::GenomeRegion("chr1",  90, 110, std2),
                      true) ;
            ASSERT_EQ(ngsai::genome::GenomeRegion("chr1", 100, 200, std1) | 
                      ngsai::genome::GenomeRegion("chr1", 100, 101, std2),
                      true) ;
            // 2nd overlaps interior
            ASSERT_EQ(ngsai::genome::GenomeRegion("chr1", 100, 200, std1) | 
                      ngsai::genome::GenomeRegion("chr1", 150, 160, std2),
                      true) ;
            ASSERT_EQ(ngsai::genome::GenomeRegion("chr1", 100, 200, std1) | 
                      ngsai::genome::GenomeRegion("chr1", 100, 200, std2),
                      true) ;    
            // 2nd overlap end
            ASSERT_EQ(ngsai::genome::GenomeRegion("chr1", 100, 200, std1) | 
                      ngsai::genome::GenomeRegion("chr1", 199, 200, std2),
                      true) ;
            ASSERT_EQ(ngsai::genome::GenomeRegion("chr1", 100, 200, std1) | 
                      ngsai::genome::GenomeRegion("chr1", 190, 210, std2),
                      true) ;
            ASSERT_EQ(ngsai::genome::GenomeRegion("chr1", 100, 200, std1) | 
                      ngsai::genome::GenomeRegion("chr1", 199, 210, std2),
                      true) ;
            
            // 2nd is downstream
            ASSERT_EQ(ngsai::genome::GenomeRegion("chr1", 100, 200, std1) | 
                      ngsai::genome::GenomeRegion("chr1", 200, 300, std2),
                      false) ;
            ASSERT_EQ(ngsai::genome::GenomeRegion("chr1", 100, 200, std1) | 
                      ngsai::genome::GenomeRegion("chr1", 200, 201, std2),
                      false) ;
            ASSERT_EQ(ngsai::genome::GenomeRegion("chr1", 100, 200, std1) | 
                      ngsai::genome::GenomeRegion("chr1", 250, 260, std2),
                      false) ;

            // not on same chromosome
            ASSERT_EQ(ngsai::genome::GenomeRegion("chr1", 100, 200, std1) | 
                      ngsai::genome::GenomeRegion("chr2", 150, 160, std2),
                      false) ;
        }
    }
}


/*!
 * \brief Checks that the downstream operator works properly.
 */
TEST(GenomeRegionTest, downstream_operator)
{
    std::vector<ngsai::genome::strand> strands =
                    {ngsai::genome::strand::FORWARD,
                     ngsai::genome::strand::REVERSE,
                     ngsai::genome::strand::UNORIENTED} ;

    for(auto std1 : strands)
    {   for(auto std2 : strands)
        {
            // 2nd is upstream
            ASSERT_EQ(ngsai::genome::GenomeRegion("chr1", 100, 200, std1) < 
                      ngsai::genome::GenomeRegion("chr1",  50,  60, std2),
                      ngsai::bool_ext::False) ;
            ASSERT_EQ(ngsai::genome::GenomeRegion("chr1", 100, 200, std1) < 
                      ngsai::genome::GenomeRegion("chr1",  90, 100, std2),
                      ngsai::bool_ext::False) ;
            ASSERT_EQ(ngsai::genome::GenomeRegion("chr1", 100, 200, std1) < 
                      ngsai::genome::GenomeRegion("chr1",  99, 100, std2),
                      ngsai::bool_ext::False) ;

            // 2nd overlaps start
            ASSERT_EQ(ngsai::genome::GenomeRegion("chr1", 100, 200, std1) < 
                      ngsai::genome::GenomeRegion("chr1",  90, 101, std2),
                      ngsai::bool_ext::False) ;
            ASSERT_EQ(ngsai::genome::GenomeRegion("chr1", 100, 200, std1) < 
                      ngsai::genome::GenomeRegion("chr1",  90, 110, std2),
                      ngsai::bool_ext::False) ;
            ASSERT_EQ(ngsai::genome::GenomeRegion("chr1", 100, 200, std1) < 
                      ngsai::genome::GenomeRegion("chr1", 100, 101, std2),
                      ngsai::bool_ext::False) ;
            // 2nd overlaps interior
            ASSERT_EQ(ngsai::genome::GenomeRegion("chr1", 100, 200, std1) < 
                      ngsai::genome::GenomeRegion("chr1", 150, 160, std2),
                      ngsai::bool_ext::False) ;
            ASSERT_EQ(ngsai::genome::GenomeRegion("chr1", 100, 200, std1) < 
                      ngsai::genome::GenomeRegion("chr1", 100, 200, std2),
                      ngsai::bool_ext::False) ;    
            // 2nd overlap end
            ASSERT_EQ(ngsai::genome::GenomeRegion("chr1", 100, 200, std1) < 
                      ngsai::genome::GenomeRegion("chr1", 199, 200, std2),
                      ngsai::bool_ext::False) ;
            ASSERT_EQ(ngsai::genome::GenomeRegion("chr1", 100, 200, std1) < 
                      ngsai::genome::GenomeRegion("chr1", 190, 210, std2),
                      ngsai::bool_ext::False) ;
            ASSERT_EQ(ngsai::genome::GenomeRegion("chr1", 100, 200, std1) < 
                      ngsai::genome::GenomeRegion("chr1", 199, 210, std2),
                      ngsai::bool_ext::False) ;
            
            // 2nd is downstream
            ASSERT_EQ(ngsai::genome::GenomeRegion("chr1", 100, 200, std1) < 
                      ngsai::genome::GenomeRegion("chr1", 200, 300, std2),
                      ngsai::bool_ext::True) ;
            ASSERT_EQ(ngsai::genome::GenomeRegion("chr1", 100, 200, std1) < 
                      ngsai::genome::GenomeRegion("chr1", 200, 201, std2),
                      ngsai::bool_ext::True) ;
            ASSERT_EQ(ngsai::genome::GenomeRegion("chr1", 100, 200, std1) < 
                      ngsai::genome::GenomeRegion("chr1", 250, 260, std2),
                      ngsai::bool_ext::True) ;

            // not on same chromosome
            ASSERT_EQ(ngsai::genome::GenomeRegion("chr1", 100, 200, std1) < 
                      ngsai::genome::GenomeRegion("chr2", 150, 160, std2),
                      ngsai::bool_ext::Undefined) ;
        }
    }
}


/*!
 * \brief Checks that the upstream operator works properly.
 */
TEST(GenomeRegionTest, upstream_operator)
{
    std::vector<ngsai::genome::strand> strands =
                    {ngsai::genome::strand::FORWARD,
                     ngsai::genome::strand::REVERSE,
                     ngsai::genome::strand::UNORIENTED} ;

    for(auto std1 : strands)
    {   for(auto std2 : strands)
        {
            // 2nd is upstream
            ASSERT_EQ(ngsai::genome::GenomeRegion("chr1", 100, 200, std1) >
                      ngsai::genome::GenomeRegion("chr1",  50,  60, std2),
                      ngsai::bool_ext::True) ;
            ASSERT_EQ(ngsai::genome::GenomeRegion("chr1", 100, 200, std1) > 
                      ngsai::genome::GenomeRegion("chr1",  90, 100, std2),
                      ngsai::bool_ext::True) ;
            ASSERT_EQ(ngsai::genome::GenomeRegion("chr1", 100, 200, std1) > 
                      ngsai::genome::GenomeRegion("chr1",  99, 100, std2),
                      ngsai::bool_ext::True) ;

            // 2nd overlaps start
            ASSERT_EQ(ngsai::genome::GenomeRegion("chr1", 100, 200, std1) > 
                      ngsai::genome::GenomeRegion("chr1",  90, 101, std2),
                      ngsai::bool_ext::False) ;
            ASSERT_EQ(ngsai::genome::GenomeRegion("chr1", 100, 200, std1) > 
                      ngsai::genome::GenomeRegion("chr1",  90, 110, std2),
                      ngsai::bool_ext::False) ;
            ASSERT_EQ(ngsai::genome::GenomeRegion("chr1", 100, 200, std1) > 
                      ngsai::genome::GenomeRegion("chr1", 100, 101, std2),
                      ngsai::bool_ext::False) ;
            // 2nd overlaps interior
            ASSERT_EQ(ngsai::genome::GenomeRegion("chr1", 100, 200, std1) > 
                      ngsai::genome::GenomeRegion("chr1", 150, 160, std2),
                      ngsai::bool_ext::False) ;
            ASSERT_EQ(ngsai::genome::GenomeRegion("chr1", 100, 200, std1) > 
                      ngsai::genome::GenomeRegion("chr1", 100, 200, std2),
                      ngsai::bool_ext::False) ;    
            // 2nd overlap end
            ASSERT_EQ(ngsai::genome::GenomeRegion("chr1", 100, 200, std1) > 
                      ngsai::genome::GenomeRegion("chr1", 199, 200, std2),
                      ngsai::bool_ext::False) ;
            ASSERT_EQ(ngsai::genome::GenomeRegion("chr1", 100, 200, std1) > 
                      ngsai::genome::GenomeRegion("chr1", 190, 210, std2),
                      ngsai::bool_ext::False) ;
            ASSERT_EQ(ngsai::genome::GenomeRegion("chr1", 100, 200, std1) > 
                      ngsai::genome::GenomeRegion("chr1", 199, 210, std2),
                      ngsai::bool_ext::False) ;
            
            // 2nd is downstream
            ASSERT_EQ(ngsai::genome::GenomeRegion("chr1", 100, 200, std1) > 
                      ngsai::genome::GenomeRegion("chr1", 200, 300, std2),
                      ngsai::bool_ext::False) ;
            ASSERT_EQ(ngsai::genome::GenomeRegion("chr1", 100, 200, std1) > 
                      ngsai::genome::GenomeRegion("chr1", 200, 201, std2),
                      ngsai::bool_ext::False) ;
            ASSERT_EQ(ngsai::genome::GenomeRegion("chr1", 100, 200, std1) > 
                      ngsai::genome::GenomeRegion("chr1", 250, 260, std2),
                      ngsai::bool_ext::False) ;

            // not on same chromosome
            ASSERT_EQ(ngsai::genome::GenomeRegion("chr1", 100, 200, std1) > 
                      ngsai::genome::GenomeRegion("chr2", 150, 160, std2),
                      ngsai::bool_ext::Undefined) ;
        }
    }
}


/*!
 * \brief Checks that the size() method works properly.
 */
TEST(GenomeRegionTest, size)
{
    ngsai::genome::strand fw = ngsai::genome::strand::FORWARD ;
    ngsai::genome::strand rv = ngsai::genome::strand::REVERSE ;
    ngsai::genome::strand un = ngsai::genome::strand::UNORIENTED ;

    ASSERT_EQ(ngsai::genome::GenomeRegion().size(), 0) ;

    ASSERT_EQ(ngsai::genome::GenomeRegion("chr1", 1000, 1001, fw).size(), 1) ;
    ASSERT_EQ(ngsai::genome::GenomeRegion("chr1", 1000, 1001, rv).size(), 1) ;
    ASSERT_EQ(ngsai::genome::GenomeRegion("chr1", 1000, 1001, un).size(), 1) ;

    ASSERT_EQ(ngsai::genome::GenomeRegion("chr1", 1000, 1100, fw).size(), 100) ;
    ASSERT_EQ(ngsai::genome::GenomeRegion("chr1", 1000, 1100, rv).size(), 100) ;
    ASSERT_EQ(ngsai::genome::GenomeRegion("chr1", 1000, 1100, un).size(), 100) ;
}


/*!
 * \brief Checks that the upstream operator works properly.
 */
TEST(GenomeRegionTest, overlap_len)
{
    std::vector<ngsai::genome::strand> strands =
                    {ngsai::genome::strand::FORWARD,
                     ngsai::genome::strand::REVERSE,
                     ngsai::genome::strand::UNORIENTED} ;

    for(auto std1 : strands)
    {   for(auto std2 : strands)
        {
            // 2nd is upstream
            ASSERT_EQ(ngsai::genome::GenomeRegion("chr1", 100, 200, std1).overlap_len(
                      ngsai::genome::GenomeRegion("chr1",  50,  60, std2)),
                      0) ;
            ASSERT_EQ(ngsai::genome::GenomeRegion("chr1", 100, 200, std1).overlap_len( 
                      ngsai::genome::GenomeRegion("chr1",  90, 100, std2)),
                      0) ;
            ASSERT_EQ(ngsai::genome::GenomeRegion("chr1", 100, 200, std1).overlap_len(
                      ngsai::genome::GenomeRegion("chr1",  99, 100, std2)),
                      0) ;

            // 2nd overlaps start
            ASSERT_EQ(ngsai::genome::GenomeRegion("chr1", 100, 200, std1).overlap_len(
                      ngsai::genome::GenomeRegion("chr1",  90, 101, std2)),
                      1) ;
            ASSERT_EQ(ngsai::genome::GenomeRegion("chr1", 100, 200, std1).overlap_len(
                      ngsai::genome::GenomeRegion("chr1",  90, 110, std2)),
                      10) ;
            ASSERT_EQ(ngsai::genome::GenomeRegion("chr1", 100, 200, std1).overlap_len(
                      ngsai::genome::GenomeRegion("chr1", 100, 101, std2)),
                      1) ;
            // 2nd overlaps interior
            ASSERT_EQ(ngsai::genome::GenomeRegion("chr1", 100, 200, std1).overlap_len(
                      ngsai::genome::GenomeRegion("chr1", 150, 160, std2)),
                      10) ;
            ASSERT_EQ(ngsai::genome::GenomeRegion("chr1", 100, 200, std1).overlap_len(
                      ngsai::genome::GenomeRegion("chr1", 100, 200, std2)),
                      100) ;    
            // 2nd overlap end
            ASSERT_EQ(ngsai::genome::GenomeRegion("chr1", 100, 200, std1).overlap_len(
                      ngsai::genome::GenomeRegion("chr1", 199, 200, std2)),
                      1) ;
            ASSERT_EQ(ngsai::genome::GenomeRegion("chr1", 100, 200, std1).overlap_len(
                      ngsai::genome::GenomeRegion("chr1", 190, 210, std2)),
                      10) ;
            ASSERT_EQ(ngsai::genome::GenomeRegion("chr1", 100, 200, std1).overlap_len(
                      ngsai::genome::GenomeRegion("chr1", 199, 210, std2)),
                      1) ;
            
            // 2nd is downstream
            ASSERT_EQ(ngsai::genome::GenomeRegion("chr1", 100, 200, std1).overlap_len(
                      ngsai::genome::GenomeRegion("chr1", 200, 300, std2)),
                      0) ;
            ASSERT_EQ(ngsai::genome::GenomeRegion("chr1", 100, 200, std1).overlap_len(
                      ngsai::genome::GenomeRegion("chr1", 200, 201, std2)),
                      0) ;
            ASSERT_EQ(ngsai::genome::GenomeRegion("chr1", 100, 200, std1).overlap_len(
                      ngsai::genome::GenomeRegion("chr1", 250, 260, std2)),
                      0) ;

            // not on same chromosome
            ASSERT_EQ(ngsai::genome::GenomeRegion("chr1", 100, 200, std1).overlap_len(
                      ngsai::genome::GenomeRegion("chr2", 150, 160, std2)),
                      0) ;
        }
    }
}


/*!
 * \brief Checks that the toString() method works properly.
 */
TEST(GenomeRegionTest, toString)
{
    ngsai::genome::GenomeRegion region1("chr1",
                                        1000,
                                        1500,
                                        ngsai::genome::strand::FORWARD) ;
    ngsai::genome::GenomeRegion region2("chr1",
                                        1000,
                                        1500,
                                        ngsai::genome::strand::REVERSE) ;
    ngsai::genome::GenomeRegion region3("chr1",
                                        1000,
                                        1500,
                                        ngsai::genome::strand::UNORIENTED) ;

    ASSERT_EQ(region1.toString(), std::string("chr1\t1000\t1500\t+")) ;
    ASSERT_EQ(region2.toString(), std::string("chr1\t1000\t1500\t-")) ;
    ASSERT_EQ(region3.toString(), std::string("chr1\t1000\t1500\t.")) ;
}