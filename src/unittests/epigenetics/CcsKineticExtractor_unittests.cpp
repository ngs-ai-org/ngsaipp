#include <gtest/gtest.h>
#include <gmock/gmock.h>
#include <string>
#include <vector>
#include <pbbam/GenomicInterval.h>
#include <pbbam/CompositeBamReader.h>  // GenomicIntervalCompositeBamReader
#include <pbbam/BamRecord.h>

#include <ngsaipp/io/BedReader.hpp>
#include <ngsaipp/epigenetics/CcsKineticExtractor.hpp>
#include <ngsaipp/genome/constants.hpp>      // ngsai::genome::strand



using testing::ElementsAre ;
using testing::ContainerEq ;



/*
 * This file ontains tests for the ngsai::CcsKineticExtract class from 
 * src/epigenetics/CcsKineticExtractor.cpp
 */


/*!
 * \brief Tests that constructor works properly and that no data have been 
 * extracted yet.
 */
TEST(CcsKineticExtractorTest, constructor)
{
    ngsai::CcsKineticExtractor extractor ;
    
    ASSERT_EQ(extractor.getSequence(), std::string()) ;
    ASSERT_EQ(extractor.getIPD(), std::vector<uint16_t>()) ;
    ASSERT_EQ(extractor.getPWD(), std::vector<uint16_t>()) ;
}


/*!
 * \brief Checks whether the extract() method works well. A 3bp window is 
 * used. This generates alignment mismatches at the really 1st position of 
 * the window.
 */
TEST(CcsKineticExtractorTest, extract_3bp)
{  // bam readers
    std::vector<std::string> paths_bam({"../data/epigenetics/reads.bam"}) ;
    PacBio::BAM::GenomicIntervalCompositeBamReader bam_reader(paths_bam) ;
    
    // 1st CpG information with window size 3
    ngsai::BedRecord cpg1_p("chr1", 18, 21, ".", 1000, ngsai::genome::FORWARD) ;
    ngsai::BedRecord cpg1_m("chr1", 19, 22, ".", 1000, ngsai::genome::REVERSE) ;
    PacBio::BAM::GenomicInterval cpg1(cpg1_p.chrom, 
                                          cpg1_p.start,
                                          cpg1_p.end) ;
    /* get reads overlapping 1st CpG
                       CCS m1_2_3/1/ccs maps fw  
                       CCCAAATCGAAATCCC
                          |||X||||||
                       CCS m1_2_3/1m/ccs maps fw
                       CCCAAAACGAAATCCC
                          ||||||||||
    REF 5' AAAAAAAAAAAAAAAAAAACGAAATAAAAAAAAAAAAAACGAAATAAAAAAAAAAAAAAA 3' fw
        3' TTTTTTTTTTTTTTTTTTTGCTTTATTTTTTTTTTTTTTGCTTTATTTTTTTTTTTTTTT 5' rv
                          ||||||||||
                       CCCTTTTGCTTTACCC
                       CCS 2 maps rv
                          ||||||X|||
                       CCCTTTTGCATTACCC
                       CCS 2m maps rv
    */
    PacBio::BAM::BamRecord read ;
    std::vector<PacBio::BAM::BamRecord> reads ;
    bam_reader.Interval(cpg1) ;
    while(bam_reader.GetNext(read))
    {   reads.push_back(read) ; }
    ASSERT_EQ(reads.size(), 4) ;
    ASSERT_EQ(reads[0].FullName(), "m1_2_3/1/ccs") ;
    ASSERT_EQ(reads[1].FullName(), "m1_2_3/1m/ccs") ;
    ASSERT_EQ(reads[2].FullName(), "m1_2_3/2/ccs") ;
    ASSERT_EQ(reads[3].FullName(), "m1_2_3/2m/ccs") ;
    
    // extractor
    ngsai::CcsKineticExtractor extractor ;
    
    // // extract kinetics of m1_2_3/1/ccs
    // // on + strand    
    std::string seq("ACG") ;
    std::vector<uint16_t> ipds = {1,9,4} ;
    std::vector<uint16_t> pwds = {1,4,3} ;
    ASSERT_EQ(extractor.extract(reads[0], cpg1_p), true) ;
    EXPECT_THAT(extractor.getSequence(), seq) ;
    EXPECT_THAT(extractor.getIPD(), ContainerEq(ipds)) ;
    EXPECT_THAT(extractor.getPWD(), ContainerEq(pwds)) ;
    // on - strand    
    seq = ("TCG") ;
    ipds = {1,9,4} ;
    pwds = {1,4,3} ;
    ASSERT_EQ(extractor.extract(reads[0], cpg1_m), true) ;
    EXPECT_THAT(extractor.getSequence(), seq) ;
    EXPECT_THAT(extractor.getIPD(), ContainerEq(ipds)) ;
    EXPECT_THAT(extractor.getPWD(), ContainerEq(pwds)) ;
 
    // extract kinetics m1_2_3/1m/ccs
    // CCS on + strand
    // this must fail, there is a missmatch at 1st pos in window alignment  
    // 5'- T C G -3' read
    //     X | |
    // 5'- A C G -3' ref (fw)
    seq = "" ;
    ipds = std::vector<uint16_t>() ;
    pwds = std::vector<uint16_t>() ;
    ASSERT_EQ(extractor.extract(reads[1], cpg1_p), false) ;
    EXPECT_THAT(extractor.getSequence(), seq) ;
    EXPECT_THAT(extractor.getIPD(), ContainerEq(ipds)) ;
    EXPECT_THAT(extractor.getPWD(), ContainerEq(pwds)) ;
    // on - strand
    seq = ("TCG") ;
    ipds = {1,9,4} ;
    pwds = {1,4,3} ;
    ASSERT_EQ(extractor.extract(reads[0], cpg1_m), true) ;
    EXPECT_THAT(extractor.getSequence(), seq) ;
    EXPECT_THAT(extractor.getIPD(), ContainerEq(ipds)) ;
    EXPECT_THAT(extractor.getPWD(), ContainerEq(pwds)) ;
 
    // extract kinetics of m1_2_3/2/ccs
    // on + strand    
    seq = "ACG" ;
    ipds = {1,9,4} ;
    pwds = {1,4,3} ;
    ASSERT_EQ(extractor.extract(reads[2], cpg1_p), true) ;
    EXPECT_THAT(extractor.getSequence(), seq) ;
    EXPECT_THAT(extractor.getIPD(), ContainerEq(ipds)) ;
    EXPECT_THAT(extractor.getPWD(), ContainerEq(pwds)) ;
    // on - strand    
    seq = ("TCG") ;
    ipds = {1,9,4} ;
    pwds = {1,4,3} ;
    ASSERT_EQ(extractor.extract(reads[2], cpg1_m), true) ;
    EXPECT_THAT(extractor.getSequence(), seq) ;
    EXPECT_THAT(extractor.getIPD(), ContainerEq(ipds)) ;
    EXPECT_THAT(extractor.getPWD(), ContainerEq(pwds)) ;
 
    // extract kinetics of m1_2_3/2m/ccs
    // on + strand    
    seq = "ACG" ;
    ipds = {1,9,4} ;
    pwds = {1,4,3} ;
    ASSERT_EQ(extractor.extract(reads[3], cpg1_p), true) ;
    EXPECT_THAT(extractor.getSequence(), seq) ;
    EXPECT_THAT(extractor.getIPD(), ContainerEq(ipds)) ;
    EXPECT_THAT(extractor.getPWD(), ContainerEq(pwds)) ;
    // on - strand
    // this must fail, there is a missmatch at 1st pos in window alignment  
    // 5'- A C G -3' read
    //     X | |
    // 5'- T C G -3' ref (rv)
    seq = "" ;
    ipds = std::vector<uint16_t>() ;
    pwds = std::vector<uint16_t>() ;
    ASSERT_EQ(extractor.extract(reads[3], cpg1_m), false) ;
    EXPECT_THAT(extractor.getSequence(), seq) ;
    EXPECT_THAT(extractor.getIPD(), ContainerEq(ipds)) ;
    EXPECT_THAT(extractor.getPWD(), ContainerEq(pwds)) ;
 
 
    // 2nd CpG information with window size 3
    ngsai::BedRecord cpg2_p("chr1", 38, 41, ".", 1000, ngsai::genome::FORWARD) ;
    ngsai::BedRecord cpg2_m("chr1", 39, 42, ".", 1000, ngsai::genome::REVERSE) ;
    PacBio::BAM::GenomicInterval cpg2(cpg2_p.chrom, 
                                      cpg2_p.start,
                                      cpg2_p.end) ;
    /* get reads overlapping 2nd CpG
                                           CCS m1_2_3/3m/ccs maps fw
                                           CCCAAATCGAAATCCC
                                              |||X||||||
                                           CCS m1_2_3/3/ccs maps fw
                                           CCCAAAACGAAATCCC
                                              ||||||||||
    REF 5' AAAAAAAAAAAAAAAAAAACGAAATAAAAAAAAAAAAAACGAAATAAAAAAAAAAAAAAA 3' fw
        3' TTTTTTTTTTTTTTTTTTTGCTTTATTTTTTTTTTTTTTGCTTTATTTTTTTTTTTTTTT 5' rv
                                              ||||||||||
                                           CCCTTTTGCTTTACCC
                                           CCS m1_2_3/4/ccs maps rv
                                              ||||||X|||
                                           CCCTTTTGCATTACCC
                                           CCS m1_2_3/4m/ccs maps rv
    */
    reads.clear() ;
    bam_reader.Interval(cpg2) ;
    while(bam_reader.GetNext(read))
    {   reads.push_back(read) ; }
    ASSERT_EQ(reads.size(), 4) ;
    ASSERT_EQ(reads[0].FullName(), "m1_2_3/3/ccs") ;
    ASSERT_EQ(reads[1].FullName(), "m1_2_3/3m/ccs") ;
    ASSERT_EQ(reads[2].FullName(), "m1_2_3/4/ccs") ;
    ASSERT_EQ(reads[3].FullName(), "m1_2_3/4m/ccs") ;
 
    // extract kinetics of m1_2_3/3/ccs
    // on + strand    
    seq = "ACG" ;
    ipds = {1,9,4} ;
    pwds = {1,4,3} ;
    ASSERT_EQ(extractor.extract(reads[0], cpg2_p), true) ;
    EXPECT_THAT(extractor.getSequence(), seq) ;
    EXPECT_THAT(extractor.getIPD(), ContainerEq(ipds)) ;
    EXPECT_THAT(extractor.getPWD(), ContainerEq(pwds)) ;
    // on - strand    
    seq = ("TCG") ;
    ipds = {1,9,4} ;
    pwds = {1,4,3} ;
    ASSERT_EQ(extractor.extract(reads[0], cpg2_m), true) ;
    EXPECT_THAT(extractor.getSequence(), seq) ;
    EXPECT_THAT(extractor.getIPD(), ContainerEq(ipds)) ;
    EXPECT_THAT(extractor.getPWD(), ContainerEq(pwds)) ;
 
    // extract kinetics of m1_2_3/3m/ccs
    // on + strand    
    // this must fail, there is a missmatch at 1st pos in window alignment  
    // 5'- T C G -3' read
    //     X | |
    // 5'- A C G -3' ref (fw)
    seq = "" ;
    ipds = std::vector<uint16_t>() ;
    pwds = std::vector<uint16_t>() ;
    ASSERT_EQ(extractor.extract(reads[1], cpg2_p), false) ;
    EXPECT_THAT(extractor.getSequence(), seq) ;
    EXPECT_THAT(extractor.getIPD(), ContainerEq(ipds)) ;
    EXPECT_THAT(extractor.getPWD(), ContainerEq(pwds)) ;
    // on - strand    
    seq = ("TCG") ;
    ipds = {1,9,4} ;
    pwds = {1,4,3} ;
    ASSERT_EQ(extractor.extract(reads[1], cpg2_m), true) ;
    EXPECT_THAT(extractor.getSequence(), seq) ;
    EXPECT_THAT(extractor.getIPD(), ContainerEq(ipds)) ;
    EXPECT_THAT(extractor.getPWD(), ContainerEq(pwds)) ;
 
    // extract kinetics of m1_2_3/4/ccs
    // on + strand    
    seq = "ACG" ;
    ipds = {1,9,4} ;
    pwds = {1,4,3} ;
    ASSERT_EQ(extractor.extract(reads[2], cpg2_p), true) ;
    EXPECT_THAT(extractor.getSequence(), seq) ;
    EXPECT_THAT(extractor.getIPD(), ContainerEq(ipds)) ;
    EXPECT_THAT(extractor.getPWD(), ContainerEq(pwds)) ;
    // on - strand    
    seq = ("TCG") ;
    ipds = {1,9,4} ;
    pwds = {1,4,3} ;
    ASSERT_EQ(extractor.extract(reads[2], cpg2_m), true) ;
    EXPECT_THAT(extractor.getSequence(), seq) ;
    EXPECT_THAT(extractor.getIPD(), ContainerEq(ipds)) ;
    EXPECT_THAT(extractor.getPWD(), ContainerEq(pwds)) ;
 
    // extract kinetics of m1_2_3/4m/ccs
    // on + strand    
    seq = "ACG" ;
    ipds = {1,9,4} ;
    pwds = {1,4,3} ;
    ASSERT_EQ(extractor.extract(reads[3], cpg2_p), true) ;
    EXPECT_THAT(extractor.getSequence(), seq) ;
    EXPECT_THAT(extractor.getIPD(), ContainerEq(ipds)) ;
    EXPECT_THAT(extractor.getPWD(), ContainerEq(pwds)) ;
    // on - strand
    // this must fail, there is a missmatch at 1st pos in window alignment  
    // 5'- a C G -3' read
    //     X | |
    // 5'- T C G -3' ref (rv)
    seq = "" ;
    ipds = std::vector<uint16_t>() ;
    pwds = std::vector<uint16_t>() ;
    ASSERT_EQ(extractor.extract(reads[3], cpg2_m), false) ;
    EXPECT_THAT(extractor.getSequence(), seq) ;
    EXPECT_THAT(extractor.getIPD(), ContainerEq(ipds)) ;
    EXPECT_THAT(extractor.getPWD(), ContainerEq(pwds)) ;
 }


/*!
 * \brief Checks whether the extract() method works well. A 9bp window is 
 * used. This generates alignment mismatches inside the window.
 */
TEST(CcsKineticExtractorTest, extract_9bp)
{  
    // bam readers
    std::vector<std::string> paths_bam({"../data/epigenetics/reads.bam"}) ;
    PacBio::BAM::GenomicIntervalCompositeBamReader bam_reader(paths_bam) ;
    
    // 1st CpG information with window size 9
    ngsai::BedRecord cpg1_p("chr1", 15, 24, ".", 1000, ngsai::genome::FORWARD) ;
    ngsai::BedRecord cpg1_m("chr1", 16, 25, ".", 1000, ngsai::genome::REVERSE) ;
    PacBio::BAM::GenomicInterval cpg1(cpg1_p.chrom, 
                                          cpg1_p.start,
                                          cpg1_p.end) ;
    /* get reads overlapping 1st CpG
                       CCS m1_2_3/1/ccs maps fw  
                       CCCAAATCGAAATCCC
                          |||X||||||
                       CCS m1_2_3/1m/ccs maps fw
                       CCCAAAACGAAATCCC
                          ||||||||||
    REF 5' AAAAAAAAAAAAAAAAAAACGAAATAAAAAAAAAAAAAACGAAATAAAAAAAAAAAAAAA 3' fw
        3' TTTTTTTTTTTTTTTTTTTGCTTTATTTTTTTTTTTTTTGCTTTATTTTTTTTTTTTTTT 5' rv
                          ||||||||||
                       CCCTTTTGCTTTACCC
                       CCS 2 maps rv
                          ||||||X|||
                       CCCTTTTGCATTACCC
                       CCS 2m maps rv
    */
    PacBio::BAM::BamRecord read ;
    std::vector<PacBio::BAM::BamRecord> reads ;
    bam_reader.Interval(cpg1) ;
    while(bam_reader.GetNext(read))
    {   reads.push_back(read) ; }
    ASSERT_EQ(reads.size(), 4) ;
    ASSERT_EQ(reads[0].FullName(), "m1_2_3/1/ccs") ;
    ASSERT_EQ(reads[1].FullName(), "m1_2_3/1m/ccs") ;
    ASSERT_EQ(reads[2].FullName(), "m1_2_3/2/ccs") ;
    ASSERT_EQ(reads[3].FullName(), "m1_2_3/2m/ccs") ;
    
    // extractor
    ngsai::CcsKineticExtractor extractor ;

    // extract kinetics of m1_2_3/1/ccs
    // on + strand    
    std::string seq("AAAACGAAA") ;
    std::vector<uint16_t> ipds = {1,1,1,1,9,4,2,2,2} ;
    std::vector<uint16_t> pwds = {1,1,1,1,4,3,2,2,2} ;
    ASSERT_EQ(extractor.extract(reads[0], cpg1_p), true) ;
    EXPECT_THAT(extractor.getSequence(), seq) ;
    EXPECT_THAT(extractor.getIPD(), ContainerEq(ipds)) ;
    EXPECT_THAT(extractor.getPWD(), ContainerEq(pwds)) ;
    // on - strand    
    seq = "ATTTCGTTT" ;
    ipds = {1,1,1,1,9,4,2,2,2} ;
    pwds = {1,1,1,1,4,3,2,2,2} ;
    ASSERT_EQ(extractor.extract(reads[0], cpg1_m), true) ;
    EXPECT_THAT(extractor.getSequence(), seq) ;
    EXPECT_THAT(extractor.getIPD(), ContainerEq(ipds)) ;
    EXPECT_THAT(extractor.getPWD(), ContainerEq(pwds)) ;

    // extract kinetics m1_2_3/1m/ccs
    // CCS on + strand
    // this must fail, there is a missmatch in window alignment  
    // 5'- A A A T C G A A A -3' read
    //     | | | X | | | | |
    // 5'- A A A A C G A A A-3' ref (fw)
    seq = "" ;
    ipds = std::vector<uint16_t>() ;
    pwds = std::vector<uint16_t>() ;
    ASSERT_EQ(extractor.extract(reads[1], cpg1_p), false) ;
    EXPECT_THAT(extractor.getSequence(), seq) ;
    EXPECT_THAT(extractor.getIPD(), ContainerEq(ipds)) ;
    EXPECT_THAT(extractor.getPWD(), ContainerEq(pwds)) ;
    // on - strand
    // this must fail, there is a missmatch in window alignment  
    // 5'- A T T T C G A T T -3' read (rv) 
    //     | | | X | | | | |
    // 5'- A T T T C G T T T-3' ref (rv)
    seq = "" ;
    ipds = std::vector<uint16_t>() ;
    pwds = std::vector<uint16_t>() ;
    ASSERT_EQ(extractor.extract(reads[1], cpg1_m), false) ;
    EXPECT_THAT(extractor.getSequence(), seq) ;
    EXPECT_THAT(extractor.getIPD(), ContainerEq(ipds)) ;
    EXPECT_THAT(extractor.getPWD(), ContainerEq(pwds)) ;

    
    // extract kinetics of m1_2_3/2/ccs
    // on + strand    
    seq = "AAAACGAAA" ;
    ipds = {1,1,1,1,9,4,2,2,2} ;
    pwds = {1,1,1,1,4,3,2,2,2} ;
    ASSERT_EQ(extractor.extract(reads[2], cpg1_p), true) ;
    EXPECT_THAT(extractor.getSequence(), seq) ;
    EXPECT_THAT(extractor.getIPD(), ContainerEq(ipds)) ;
    EXPECT_THAT(extractor.getPWD(), ContainerEq(pwds)) ;
    // on - strand    
    seq = "ATTTCGTTT" ;
    ipds = {1,1,1,1,9,4,2,2,2} ;
    pwds = {1,1,1,1,4,3,2,2,2} ;
    ASSERT_EQ(extractor.extract(reads[2], cpg1_m), true) ;
    EXPECT_THAT(extractor.getSequence(), seq) ;
    EXPECT_THAT(extractor.getIPD(), ContainerEq(ipds)) ;
    EXPECT_THAT(extractor.getPWD(), ContainerEq(pwds)) ;


    // extract kinetics of m1_2_3/2m/ccs
    // on + strand
    // CCS on + strand
    // this must fail, there is a missmatch in window alignment  
    // 5'- A A A T C G A A A -3' read
    //     | | | X | | | | |
    // 5'- A A A A C G A A A-3' ref (fw) 
    seq = "" ;
    ipds = std::vector<uint16_t>() ;
    pwds = std::vector<uint16_t>() ;
    ASSERT_EQ(extractor.extract(reads[3], cpg1_p), false) ;
    EXPECT_THAT(extractor.getSequence(), seq) ;
    EXPECT_THAT(extractor.getIPD(), ContainerEq(ipds)) ;
    EXPECT_THAT(extractor.getPWD(), ContainerEq(pwds)) ;
    // on - strand
   // this must fail, there is a missmatch in window alignment  
    // 5'- A T T T C G A T T -3' read (rv) 
    //     | | | X | | | | |
    // 5'- A T T T C G T T T-3' ref (rv)
    seq = "" ;
    ipds = std::vector<uint16_t>() ;
    pwds = std::vector<uint16_t>() ;
    ASSERT_EQ(extractor.extract(reads[3], cpg1_m), false) ;
    EXPECT_THAT(extractor.getSequence(), seq) ;
    EXPECT_THAT(extractor.getIPD(), ContainerEq(ipds)) ;
    EXPECT_THAT(extractor.getPWD(), ContainerEq(pwds)) ;
    

    // 2nd CpG information with window size 9
    ngsai::BedRecord cpg2_p("chr1", 35, 44, ".", 1000, ngsai::genome::FORWARD) ;
    ngsai::BedRecord cpg2_m("chr1", 36, 45, ".", 1000, ngsai::genome::REVERSE) ;
    PacBio::BAM::GenomicInterval cpg2(cpg2_p.chrom, 
                                      cpg2_p.start,
                                      cpg2_p.end) ;
    /* get reads overlapping 2nd CpG
                                           CCS m1_2_3/3m/ccs maps fw
                                           CCCAAATCGAAATCCC
                                              |||X||||||
                                           CCS m1_2_3/3/ccs maps fw
                                           CCCAAAACGAAATCCC
                                              ||||||||||
    REF 5' AAAAAAAAAAAAAAAAAAACGAAATAAAAAAAAAAAAAACGAAATAAAAAAAAAAAAAAA 3' fw
        3' TTTTTTTTTTTTTTTTTTTGCTTTATTTTTTTTTTTTTTGCTTTATTTTTTTTTTTTTTT 5' rv
                                              ||||||||||
                                           CCCTTTTGCTTTACCC
                                           CCS m1_2_3/4/ccs maps rv
                                              ||||||X|||
                                           CCCTTTTGCATTACCC
                                           CCS m1_2_3/4m/ccs maps rv

    */
    reads.clear() ;
    bam_reader.Interval(cpg2) ;
    while(bam_reader.GetNext(read))
    {   reads.push_back(read) ; }
    ASSERT_EQ(reads.size(), 4) ;
    ASSERT_EQ(reads[0].FullName(), "m1_2_3/3/ccs") ;
    ASSERT_EQ(reads[1].FullName(), "m1_2_3/3m/ccs") ;
    ASSERT_EQ(reads[2].FullName(), "m1_2_3/4/ccs") ;
    ASSERT_EQ(reads[3].FullName(), "m1_2_3/4m/ccs") ;

    // extract kinetics of m1_2_3/3/ccs
    // on + strand    
    seq = "AAAACGAAA" ;
    ipds = {1,1,1,1,9,4,2,2,2} ;
    pwds = {1,1,1,1,4,3,2,2,2} ;
    ASSERT_EQ(extractor.extract(reads[0], cpg2_p), true) ;
    EXPECT_THAT(extractor.getSequence(), seq) ;
    EXPECT_THAT(extractor.getIPD(), ContainerEq(ipds)) ;
    EXPECT_THAT(extractor.getPWD(), ContainerEq(pwds)) ;
    // on - strand    
    seq = "ATTTCGTTT" ;
    ipds = {1,1,1,1,9,4,2,2,2} ;
    pwds = {1,1,1,1,4,3,2,2,2} ;
    ASSERT_EQ(extractor.extract(reads[0], cpg2_m), true) ;
    EXPECT_THAT(extractor.getSequence(), seq) ;
    EXPECT_THAT(extractor.getIPD(), ContainerEq(ipds)) ;
    EXPECT_THAT(extractor.getPWD(), ContainerEq(pwds)) ;

    
    // extract kinetics m1_2_3/3m/ccs
    // CCS on + strand
    // this must fail, there is a missmatch in window alignment  
    // 5'- A A A T C G A A A -3' read
    //     | | | X | | | | |
    // 5'- A A A A C G A A A-3' ref (fw)
    seq = "" ;
    ipds = std::vector<uint16_t>() ;
    pwds = std::vector<uint16_t>() ;
    ASSERT_EQ(extractor.extract(reads[1], cpg2_p), false) ;
    EXPECT_THAT(extractor.getSequence(), seq) ;
    EXPECT_THAT(extractor.getIPD(), ContainerEq(ipds)) ;
    EXPECT_THAT(extractor.getPWD(), ContainerEq(pwds)) ;
    // on - strand
    // this must fail, there is a missmatch in window alignment  
    // 5'- A T T T C G A T T -3' read (rv) 
    //     | | | X | | | | |
    // 5'- A T T T C G T T T-3' ref (rv)
    seq = "" ;
    ipds = std::vector<uint16_t>() ;
    pwds = std::vector<uint16_t>() ;
    ASSERT_EQ(extractor.extract(reads[1], cpg2_m), false) ;
    EXPECT_THAT(extractor.getSequence(), seq) ;
    EXPECT_THAT(extractor.getIPD(), ContainerEq(ipds)) ;
    EXPECT_THAT(extractor.getPWD(), ContainerEq(pwds)) ;
    
    // extract kinetics of m1_2_3/4/ccs
    // on + strand    
    // on + strand    
    seq = "AAAACGAAA" ;
    ipds = {1,1,1,1,9,4,2,2,2} ;
    pwds = {1,1,1,1,4,3,2,2,2} ;
    ASSERT_EQ(extractor.extract(reads[2], cpg2_p), true) ;
    EXPECT_THAT(extractor.getSequence(), seq) ;
    EXPECT_THAT(extractor.getIPD(), ContainerEq(ipds)) ;
    EXPECT_THAT(extractor.getPWD(), ContainerEq(pwds)) ;
    // on - strand    
    seq = "ATTTCGTTT" ;
    ipds = {1,1,1,1,9,4,2,2,2} ;
    pwds = {1,1,1,1,4,3,2,2,2} ;
    ASSERT_EQ(extractor.extract(reads[2], cpg2_m), true) ;
    EXPECT_THAT(extractor.getSequence(), seq) ;
    EXPECT_THAT(extractor.getIPD(), ContainerEq(ipds)) ;
    EXPECT_THAT(extractor.getPWD(), ContainerEq(pwds)) ;

    
    // extract kinetics of m1_2_3/4m/ccs
    // on + strand
    // CCS on + strand
    // this must fail, there is a missmatch in window alignment  
    // 5'- A A A T C G A A A -3' read
    //     | | | X | | | | |
    // 5'- A A A A C G A A A-3' ref (fw) 
    seq = "" ;
    ipds = std::vector<uint16_t>() ;
    pwds = std::vector<uint16_t>() ;
    ASSERT_EQ(extractor.extract(reads[3], cpg2_p), false) ;
    EXPECT_THAT(extractor.getSequence(), seq) ;
    EXPECT_THAT(extractor.getIPD(), ContainerEq(ipds)) ;
    EXPECT_THAT(extractor.getPWD(), ContainerEq(pwds)) ;
    // on - strand
   // this must fail, there is a missmatch in window alignment  
    // 5'- A T T T C G A T T -3' read (rv) 
    //     | | | X | | | | |
    // 5'- A T T T C G T T T-3' ref (rv)
    seq = "" ;
    ipds = std::vector<uint16_t>() ;
    pwds = std::vector<uint16_t>() ;
    ASSERT_EQ(extractor.extract(reads[3], cpg2_m), false) ;
    EXPECT_THAT(extractor.getSequence(), seq) ;
    EXPECT_THAT(extractor.getIPD(), ContainerEq(ipds)) ;
    EXPECT_THAT(extractor.getPWD(), ContainerEq(pwds)) ;
}