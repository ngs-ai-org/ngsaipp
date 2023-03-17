#include <gtest/gtest.h>
#include <gmock/gmock.h>
#include <string>

#include <ngsaipp/genome/constants.hpp>


/*
 * This file contains tests for the functions defined in 
 * src/genome/constants.hpp
 */


/*!
 * \brief Checks that the function char_to_strand works properly.
 */
TEST(constantsTest, char_to_strand)
{
    ngsai::genome::strand fw = ngsai::genome::strand::FORWARD ;
    ngsai::genome::strand rv = ngsai::genome::strand::REVERSE ;
    ngsai::genome::strand un = ngsai::genome::strand::UNORIENTED ;
    
    ASSERT_EQ(ngsai::genome::char_to_strand('+'), fw) ;
    ASSERT_EQ(ngsai::genome::char_to_strand('-'), rv) ;
    ASSERT_EQ(ngsai::genome::char_to_strand('.'), un) ;
    
    // must fail because only '+', '-', '.' char are accepted
    EXPECT_THROW(ngsai::genome::char_to_strand('\n'),
                 std::invalid_argument) ;
    EXPECT_THROW(ngsai::genome::char_to_strand('z'),
                 std::invalid_argument) ;
}


/*!
 * \brief Checks that the function strand_to_char works properly.
 */
TEST(constantsTest, strand_to_char)
{
    ngsai::genome::strand fw = ngsai::genome::strand::FORWARD ;
    ngsai::genome::strand rv = ngsai::genome::strand::REVERSE ;
    ngsai::genome::strand un = ngsai::genome::strand::UNORIENTED ;
    
    ASSERT_EQ(ngsai::genome::strand_to_char(fw), '+') ;
    ASSERT_EQ(ngsai::genome::strand_to_char(rv), '-') ;
    ASSERT_EQ(ngsai::genome::strand_to_char(un), '.') ;
}
