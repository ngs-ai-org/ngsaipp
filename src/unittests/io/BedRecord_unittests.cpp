#include <gtest/gtest.h>
#include <gmock/gmock.h>
#include <string>
#include <stdexcept>
#include <utility>        // std::move()

#include <ngsaipp/io/BedRecord.hpp>


/*
 * This file contains tests for the ngsai::BedRecord class from 
 * src/io/BedRecord.cpp
 */


/*!
 * \brief Checks that the default constructor works properly.
 */
TEST(BedRecordTest, constructor_default)
{
    ngsai::BedRecord region ;
    ASSERT_EQ(region.chrom,  std::string("")) ;
    ASSERT_EQ(region.start,  0) ;
    ASSERT_EQ(region.end,    0) ;
    ASSERT_EQ(region.strand, ngsai::genome::strand::UNORIENTED) ;
    ASSERT_EQ(region.name,   std::string("")) ;
    ASSERT_EQ(region.score,  0) ;
}


/*!
 * \brief Checks that the parameter constructor works properly.
 */
TEST(BedRecordTest, constructor_param)
{
    ngsai::genome::strand fw = ngsai::genome::strand::FORWARD ;
    ngsai::genome::strand rv = ngsai::genome::strand::REVERSE ;
    ngsai::genome::strand un = ngsai::genome::strand::UNORIENTED ;

    // empty chrom
    ngsai::BedRecord  region1("", 1000, 1500, "CpG", 99, fw) ;
    ASSERT_EQ(region1.chrom,  std::string("")) ;
    ASSERT_EQ(region1.start,  1000) ;
    ASSERT_EQ(region1.end,    1500) ;
    ASSERT_EQ(region1.strand, fw) ;
    ASSERT_EQ(region1.name,   std::string("CpG")) ;
    ASSERT_EQ(region1.score,  99) ;

    // regular fw
    ngsai::BedRecord region2("chr1", 1000, 1500, "CpG", 99, fw) ;
    ASSERT_EQ(region2.chrom,  std::string("chr1")) ;
    ASSERT_EQ(region2.start,  1000) ;
    ASSERT_EQ(region2.end,    1500) ;
    ASSERT_EQ(region2.strand, fw) ;
    ASSERT_EQ(region2.name,   std::string("CpG")) ;
    ASSERT_EQ(region2.score,  99) ;

    // regular rv
    ngsai::BedRecord region3("chr1", 1000, 1500, "CpG", 99, rv) ;
    ASSERT_EQ(region3.chrom,  std::string("chr1")) ;
    ASSERT_EQ(region3.start,  1000) ;
    ASSERT_EQ(region3.end,    1500) ;
    ASSERT_EQ(region3.strand, rv) ;
    ASSERT_EQ(region3.name,   std::string("CpG")) ;
    ASSERT_EQ(region3.score,  99) ;

    // regular un
    ngsai::BedRecord region4("chr1", 1000, 1500, "CpG", 99, un) ;
    ASSERT_EQ(region4.chrom,  std::string("chr1")) ;
    ASSERT_EQ(region4.start,  1000) ;
    ASSERT_EQ(region4.end,    1500) ;
    ASSERT_EQ(region4.strand, un) ;
    ASSERT_EQ(region4.name,   std::string("CpG")) ;
    ASSERT_EQ(region4.score,  99) ;


    // should fail because start == end
    EXPECT_THROW(ngsai::BedRecord("chr1", 1000, 1000, "CpG", 99, fw),
                 std::invalid_argument) ;
    EXPECT_THROW(ngsai::BedRecord("chr1", 1000, 1000, "CpG", 99, rv),
                 std::invalid_argument) ;
    EXPECT_THROW(ngsai::BedRecord("chr1", 1000, 1000, "CpG", 99, un),
                 std::invalid_argument) ;

    // should fail because start > end
    EXPECT_THROW(ngsai::BedRecord("chr1", 1001, 1000, "CpG", 99, fw),
                 std::invalid_argument) ;
    EXPECT_THROW(ngsai::BedRecord("chr1", 1001, 1000, "CpG", 99, rv),
                 std::invalid_argument) ;
    EXPECT_THROW(ngsai::BedRecord("chr1", 1001, 1000, "CpG", 99, un),
                 std::invalid_argument) ;
}


/*!
 * \brief Checks that the copy constructor works properly.
 */
TEST(BedRecordTest, constructor_copy)
{
    ngsai::BedRecord region1("chr1",
                             1000,
                             1500,
                             "CpG",
                             99,
                             ngsai::genome::strand::FORWARD) ;
    ngsai::BedRecord region2(region1) ;

    ASSERT_EQ(region1.chrom,  region2.chrom) ;
    ASSERT_EQ(region1.start,  region2.start) ;
    ASSERT_EQ(region1.end,    region2.end) ;
    ASSERT_EQ(region1.strand, region2.strand) ;
    ASSERT_EQ(region1.name,   region2.name) ;
    ASSERT_EQ(region1.score,  region2.score) ;
}


/*!
 * \brief Checks that the move constructor works properly.
 */
TEST(BedRecordTest, constructor_move)
{
    ngsai::BedRecord region1("chr1",
                             1000,
                             1500,
                             "CpG",
                             99,
                             ngsai::genome::strand::FORWARD) ;
    ngsai::BedRecord region2(region1) ;
    ngsai::BedRecord region3(std::move(region2)) ;

    ASSERT_EQ(region1.chrom,  region3.chrom) ;
    ASSERT_EQ(region1.start,  region3.start) ;
    ASSERT_EQ(region1.end,    region3.end) ;
    ASSERT_EQ(region1.strand, region3.strand) ;
    ASSERT_EQ(region1.name,   region3.name) ;
    ASSERT_EQ(region1.score,  region3.score) ;
}


/*!
 * \brief Checks that the assignment operator works properly.
 */
TEST(BedRecordTest, assignment_operator)
{
    ngsai::BedRecord region1("chr1",
                             1000,
                             1500,
                             "CpG",
                             99,
                             ngsai::genome::strand::FORWARD) ;
    ngsai::BedRecord region2 = region1 ;

    ASSERT_EQ(region1.chrom,  region2.chrom) ;
    ASSERT_EQ(region1.start,  region2.start) ;
    ASSERT_EQ(region1.end,    region2.end) ;
    ASSERT_EQ(region1.strand, region2.strand) ;
    ASSERT_EQ(region1.name,   region2.name) ;
    ASSERT_EQ(region1.score,  region2.score) ;
}


/*!
 * \brief Checks that the move assignment operator works properly.
 */
TEST(BedRecordTest, assignment_move_operator)
{
    ngsai::BedRecord region1("chr1",
                             1000,
                             1500,
                             "CpG",
                             99,
                             ngsai::genome::strand::FORWARD) ;
    ngsai::BedRecord region2 = region1 ;
    ngsai::BedRecord region3(std::move(region1)) ;

    ASSERT_EQ(region3.chrom,  region2.chrom) ;
    ASSERT_EQ(region3.start,  region2.start) ;
    ASSERT_EQ(region3.end,    region2.end) ;
    ASSERT_EQ(region3.strand, region2.strand) ;
    ASSERT_EQ(region3.name,   region2.name) ;
    ASSERT_EQ(region3.score,  region2.score) ;
}


/*!
 * \brief Checks that the equality comparison operator works properly.
 */
TEST(BedRecordTest, equality_operator)
{   ngsai::BedRecord region1("chr1",
                             1000,
                             1500,
                             "CpG",
                             99,
                             ngsai::genome::strand::FORWARD) ;
    ngsai::BedRecord region2("chr1",
                              1000,
                              1500,
                              "CpG",
                              99,
                              ngsai::genome::strand::FORWARD) ;
    ngsai::BedRecord region3("chr2",
                             1000,
                             1500,
                             "CpG",
                             99,
                             ngsai::genome::strand::FORWARD) ;
    ngsai::BedRecord region4("chr1",
                             1001,
                             1500,
                             "CpG",
                             99,
                             ngsai::genome::strand::FORWARD) ;
    ngsai::BedRecord region5("chr1",
                             1000,
                             1501,
                             "CpG",
                             99,
                             ngsai::genome::strand::FORWARD) ;
    ngsai::BedRecord region6("chr1",
                             1000,
                             1500,
                             "CpG",
                             99,
                             ngsai::genome::strand::REVERSE) ;
    ngsai::BedRecord region7("chr1",
                             1000,
                             1500,
                             "gene",
                             99,
                             ngsai::genome::strand::FORWARD) ;
    ngsai::BedRecord region8("chr1",
                              1000,
                              1500,
                              "CpG",
                              11,
                              ngsai::genome::strand::FORWARD) ;

    ASSERT_EQ(region1 == region2, true) ;
    ASSERT_EQ(region1 == region3, false) ;
    ASSERT_EQ(region1 == region4, false) ;
    ASSERT_EQ(region1 == region5, false) ;
    ASSERT_EQ(region1 == region6, false) ;
    ASSERT_EQ(region1 == region7, false) ;
    ASSERT_EQ(region1 == region8, false) ;
}


/*!
 * \brief Checks that the inequality comparison operator works properly.
 */
TEST(BedRecordTest, inequality_operator)
{   ngsai::BedRecord region1("chr1",
                             1000,
                             1500,
                             "CpG",
                             99,
                             ngsai::genome::strand::FORWARD) ;
    ngsai::BedRecord region2("chr1",
                              1000,
                              1500,
                              "CpG",
                              99,
                              ngsai::genome::strand::FORWARD) ;
    ngsai::BedRecord region3("chr2",
                             1000,
                             1500,
                             "CpG",
                             99,
                             ngsai::genome::strand::FORWARD) ;
    ngsai::BedRecord region4("chr1",
                             1001,
                             1500,
                             "CpG",
                             99,
                             ngsai::genome::strand::FORWARD) ;
    ngsai::BedRecord region5("chr1",
                             1000,
                             1501,
                             "CpG",
                             99,
                             ngsai::genome::strand::FORWARD) ;
    ngsai::BedRecord region6("chr1",
                             1000,
                             1500,
                             "CpG",
                             99,
                             ngsai::genome::strand::REVERSE) ;
    ngsai::BedRecord region7("chr1",
                             1000,
                             1500,
                             "gene",
                             99,
                             ngsai::genome::strand::FORWARD) ;
    ngsai::BedRecord region8("chr1",
                              1000,
                              1500,
                              "CpG",
                              11,
                              ngsai::genome::strand::FORWARD) ;

    ASSERT_EQ(region1 != region2, false) ;
    ASSERT_EQ(region1 != region3, true) ;
    ASSERT_EQ(region1 != region4, true) ;
    ASSERT_EQ(region1 != region5, true) ;
    ASSERT_EQ(region1 != region6, true) ;
    ASSERT_EQ(region1 != region7, true) ;
    ASSERT_EQ(region1 != region8, true) ;
}


/*!
 * \brief Checks that the toString() method works properly.
 */
TEST(BedRecordTest, toString)
{
    ngsai::BedRecord region1("chr1",
                             1000,
                             1500,
                             "CpG",
                             99,
                             ngsai::genome::strand::FORWARD) ;
    ngsai::BedRecord region2("chr1",
                             1000,
                             1500,
                             "CpG",
                             99,
                             ngsai::genome::strand::REVERSE) ;
    ngsai::BedRecord region3("chr1",
                             1000,
                             1500,
                             "CpG",
                             99,
                             ngsai::genome::strand::UNORIENTED) ;

    ASSERT_EQ(region1.toString(), std::string("chr1\t1000\t1500\tCpG\t99\t+")) ;
    ASSERT_EQ(region2.toString(), std::string("chr1\t1000\t1500\tCpG\t99\t-")) ;
    ASSERT_EQ(region3.toString(), std::string("chr1\t1000\t1500\tCpG\t99\t.")) ;
}