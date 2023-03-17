#include <gtest/gtest.h>
#include <gmock/gmock.h>
#include <string>
#include <stdexcept>
#include <utility>        // std::move()

#include <ngsaipp/io/BedReader.hpp>
#include <ngsaipp/io/BedRecord.hpp>
#include <ngsaipp/genome/constants.hpp>


/*
 * This file contains tests for the ngsai::BedReader class from 
 * src/io/BedReader.cpp
 */


std::string path_good_bed6("../data/io/good_bed6.bed") ;

std::string path_empty("../data/io/empty.bed") ;
std::string path_empty_eol("../data/io/empty_eol.bed") ;

std::string path_empty_1st_line("../data/io/bad_empty_1st_line.bed") ;
std::string path_empty_mid_line("../data/io/bad_empty_mid_line.bed") ;
std::string path_empty_last_line("../data/io/bad_empty_last_line.bed") ;

std::string path_bad_format_line1_1("../data/io/bad_format_line1_1.bed") ;
std::string path_bad_format_line1_2("../data/io/bad_format_line1_2.bed") ;
std::string path_bad_format_line1_3("../data/io/bad_format_line1_3.bed") ;

std::string path_bad_format_line2_1("../data/io/bad_format_line2_1.bed") ;
std::string path_bad_format_line2_2("../data/io/bad_format_line2_2.bed") ;
std::string path_bad_format_line2_3("../data/io/bad_format_line2_3.bed") ;

std::string path_bad_format_line3_1("../data/io/bad_format_line3_1.bed") ;
std::string path_bad_format_line3_2("../data/io/bad_format_line3_2.bed") ;
std::string path_bad_format_line3_3("../data/io/bad_format_line3_3.bed") ;


ngsai::genome::strand fw = ngsai::genome::strand::FORWARD ;
ngsai::genome::strand rv = ngsai::genome::strand::REVERSE ;
ngsai::genome::strand un = ngsai::genome::strand::UNORIENTED ;

ngsai::BedRecord line_1("chr1", 1000, 5000, "cloneA", 960, fw) ;
ngsai::BedRecord line_2("chr1", 2000, 6000, "cloneA", 960, rv) ;
ngsai::BedRecord line_3("chr1", 3000, 7000, "cloneA", 960, un) ;


/*!
 * \brief The BedReaderTest class is a class made for testing purpose. It 
 * provides an API to get references to the internal private fields of 
 * BedReader.
 */
class BedReaderTest : public ngsai::BedReader
{   public:
        BedReaderTest(const std::string& path)
            : BedReader(path)
        { ; }

        ~BedReaderTest()
        { ; }

        std::string& getMPath()
        {   return m_path ; }

        std::ifstream& getMFBed()
        {   return m_f_bed ; }

        bool& getMIsOpen()
        {   return m_is_open ; }

} ;

/*!
 * \brief Checks that the default constructor works properly.
 */
TEST(BedReaderTest, constructor_argument)
{
    BedReaderTest reader(path_good_bed6) ;
    ASSERT_EQ(reader.getMPath(),          path_good_bed6) ;
    ASSERT_EQ(reader.getMFBed().good(),   true) ;
    ASSERT_EQ(reader.getMFBed().tellg(),   0) ;
    ASSERT_EQ(reader.getMIsOpen(),        true) ;
}



/*!
 * \brief Checks it can read a well formatted bed6 file.
 */
TEST(BedReaderTest, read_good_bed6)
{   
    BedReaderTest reader(path_good_bed6) ;
    
    ngsai::BedRecord line ;
    std::vector<ngsai::BedRecord> lines ;
    while(reader.getNext(line))
    {   lines.push_back(line) ; }

    ASSERT_EQ(lines.size(), 3) ;
    ASSERT_EQ(lines[0], line_1) ;
    ASSERT_EQ(lines[1], line_2) ;
    ASSERT_EQ(lines[2], line_3) ;
}


/*!
 * \brief Checks that nothing happens when trying to read an empty file.
 */
TEST(BedReaderTest, read_empty)
{   
    BedReaderTest reader(path_empty) ;
    ngsai::BedRecord line ;

    ASSERT_EQ(reader.getNext(line), false) ;
}


/*!
 * \brief Checks that it throws an error when reading an empty file containing 
 * only an EOL char.
 */
TEST(BedReaderTest, read_empty_eol)
{   
    BedReaderTest reader(path_empty_eol) ;
    ngsai::BedRecord line ;

    EXPECT_THROW(reader.getNext(line), std::runtime_error) ;
}


/*!
 * \brief Checks that it throws an error when reading a file containing 
 * a 1st empty line (otherwise well formatted).
 */
TEST(BedReaderTest, read_empty_1st_line)
{   
    BedReaderTest reader(path_empty_1st_line) ;
    ngsai::BedRecord line ;

    EXPECT_THROW(reader.getNext(line), std::runtime_error) ;
}


/*!
 * \brief Checks that it throws an error when reading a file containing 
 * a mid empty line (otherwise well formatted).
 */
TEST(BedReaderTest, read_empty_mid_line)
{   
    ngsai::BedRecord line ;
    BedReaderTest reader(path_empty_mid_line) ;
    
    reader.getNext(line) ;
    ASSERT_EQ(line, line_1) ;

    EXPECT_THROW(reader.getNext(line), std::runtime_error) ;
}


/*!
 * \brief Checks that it throws an error when reading a file containing 
 * a mid empty line (otherwise well formatted).
 */
TEST(BedReaderTest, read_empty_last_line)
{   
    ngsai::BedRecord line ;
    BedReaderTest reader(path_empty_last_line) ;
    
    reader.getNext(line) ;
    ASSERT_EQ(line, line_1) ;

    reader.getNext(line) ;
    ASSERT_EQ(line, line_2) ;

    reader.getNext(line) ;
    ASSERT_EQ(line, line_3) ;
    
    EXPECT_THROW(reader.getNext(line), std::runtime_error) ;
}


/*!
 * \brief Checks that it throws an error when reading a file containing 
 * a 1st line with a format error (otherwise well formatted).
 */
TEST(BedReaderTest, read_bad_format_line1_1)
{   
    BedReaderTest reader(path_bad_format_line1_1) ;
    ngsai::BedRecord line ;

    EXPECT_THROW(reader.getNext(line), std::runtime_error) ;
}


/*!
 * \brief Checks that it throws an error when reading a file containing 
 * a 1st line with a format error (otherwise well formatted).
 */
TEST(BedReaderTest, read_bad_format_line1_2)
{   
    BedReaderTest reader(path_bad_format_line1_1) ;
    ngsai::BedRecord line ;

    EXPECT_THROW(reader.getNext(line), std::runtime_error) ;
}


/*!
 * \brief Checks that it throws an error when reading a file containing 
 * a 1st line with a format error (otherwise well formatted).
 */
TEST(BedReaderTest, read_bad_format_line1_3)
{   
    BedReaderTest reader(path_bad_format_line1_3) ;
    ngsai::BedRecord line ;

    EXPECT_THROW(reader.getNext(line), std::runtime_error) ;
}


/*!
 * \brief Checks that it throws an error when reading a file containing 
 * a 2nd line with a format error (otherwise well formatted).
 */
TEST(BedReaderTest, read_bad_format_line2_1)
{   
    BedReaderTest reader(path_bad_format_line2_1) ;
    
    ngsai::BedRecord line ;

    reader.getNext(line) ;
    ASSERT_EQ(line, line_1) ;

    EXPECT_THROW(reader.getNext(line), std::runtime_error) ;
}


/*!
 * \brief Checks that it throws an error when reading a file containing 
 * a 2nd line with a format error (otherwise well formatted).
 */
TEST(BedReaderTest, read_bad_format_line2_2)
{   
    BedReaderTest reader(path_bad_format_line2_2) ;
    
    ngsai::BedRecord line ;

    reader.getNext(line) ;
    ASSERT_EQ(line, line_1) ;

    EXPECT_THROW(reader.getNext(line), std::runtime_error) ;
}


/*!
 * \brief Checks that it throws an error when reading a file containing 
 * a 2nd line with a format error (otherwise well formatted).
 */
TEST(BedReaderTest, read_bad_format_line2_3)
{   
    BedReaderTest reader(path_bad_format_line2_3) ;
    
    ngsai::BedRecord line ;

    reader.getNext(line) ;
    ASSERT_EQ(line, line_1) ;

    EXPECT_THROW(reader.getNext(line), std::runtime_error) ;
}


/*!
 * \brief Checks that it throws an error when reading a file containing 
 * a 3rd line with a format error (otherwise well formatted).
 */
TEST(BedReaderTest, read_bad_format_line3_1)
{   
    BedReaderTest reader(path_bad_format_line3_1) ;

    ngsai::BedRecord line ;

    reader.getNext(line) ;
    ASSERT_EQ(line, line_1) ;

    reader.getNext(line) ;
    ASSERT_EQ(line, line_2) ;

    EXPECT_THROW(reader.getNext(line), std::runtime_error) ;
}


/*!
 * \brief Checks that it throws an error when reading a file containing 
 * a 3rd line with a format error (otherwise well formatted).
 */
TEST(BedReaderTest, read_bad_format_line3_2)
{   
    BedReaderTest reader(path_bad_format_line3_2) ;

    ngsai::BedRecord line ;

    reader.getNext(line) ;
    ASSERT_EQ(line, line_1) ;

    reader.getNext(line) ;
    ASSERT_EQ(line, line_2) ;

    EXPECT_THROW(reader.getNext(line), std::runtime_error) ;
}


/*!
 * \brief Checks that it throws an error when reading a file containing 
 * a 3rd line with a format error (otherwise well formatted).
 */
TEST(BedReaderTest, read_bad_format_line3_3)
{   
    BedReaderTest reader(path_bad_format_line3_3) ;

    ngsai::BedRecord line ;

    reader.getNext(line) ;
    ASSERT_EQ(line, line_1) ;

    reader.getNext(line) ;
    ASSERT_EQ(line, line_2) ;

    EXPECT_THROW(reader.getNext(line), std::runtime_error) ;
}