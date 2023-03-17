#include <gtest/gtest.h>
#include <gmock/gmock.h>
#include <string>
#include <vector>
#include <stdexcept>
#include <limits>
#include <cmath>

#include "utility.hpp"  // kmerMapEqual()
#include <ngsaipp/epigenetics/NormalizedKineticModel.hpp>
#include <ngsaipp/epigenetics/KineticSignal.hpp>


using testing::ElementsAre ;
using testing::ContainerEq ;


/*
 * This file ontains tests for the ngsai::NormalizedKineticModel class from 
 * src/epigenetics/NormalizedKineticModel.cpp
 */


/*!
 * \brief A child class of NormlaizedKineticModel allowing to access the 
 * protected fields of NormalizedKineticModel for testing purposes.
 */
class NormalizedKineticModelToy : public ngsai::NormalizedKineticModel
{
    public:

        NormalizedKineticModelToy(const ngsai::KmerMap& bckg_model)
            : NormalizedKineticModel(bckg_model)
        { ; }

        NormalizedKineticModelToy(const ngsai::NormalizedKineticModel& other)
            : NormalizedKineticModel(other)
        { ; }

        NormalizedKineticModelToy(ngsai::NormalizedKineticModel&& other)
            : NormalizedKineticModel(other)
        { ; }

        ~NormalizedKineticModelToy()
        { ; }

        // ugly workaround to test NormalizedKineticModel.copy() using a
        // NormalizedKineticModelToy instance
        NormalizedKineticModelToy* copy() const
        {   NormalizedKineticModel* copy = NormalizedKineticModel::copy() ;
            NormalizedKineticModelToy* copy2 = 
                                        new NormalizedKineticModelToy(*copy) ;
            delete copy ;
            return copy2 ;
        }

        std::vector<hist_1d_double>& getHistogramsIPD()
        {   return m_histograms_ipd ; }

        std::vector<hist_1d_double>& getHistogramsPWD()
        {   return m_histograms_pwd ; }

        ngsai::KmerMap& getBckgModel()
        {   return m_bckg_model ; }

        bool& getIsInit()
        {   return m_is_init ; }

        bool& getIsDensity()
        {   return m_is_density ; }

        bool& getIsLog()
        {   return m_is_log ; }

} ;


/*!
 * \brief Contains the necessary KmerMap that will be used by each and every 
 * NormalizedKineticModelToy instance.
 * Rule is simple, the center position get an expected IPD/PDW of :
 * 1 if there is A in the center
 * 2 if there is C in the center
 * 3 if there is G in the center
 * 4 if there is T in the center
 */
class NormalizedKineticModelTest : public ::testing::Test
{   
    protected:

        virtual void 
        SetUp() override
        { ; }

        virtual void 
        TearDown() override
        { ; }

    protected:

        static ngsai::KmerMap* kmermap ;

        static void 
        SetUpTestSuite()
        {   
            if(kmermap == nullptr)
            {   std::vector<uint32_t> a({1,1,1}) ;  // values if A is in center
                std::vector<uint32_t> c({1,2,1}) ;  // values if C is in center
                std::vector<uint32_t> g({1,3,1}) ;  // values if G is in center
                std::vector<uint32_t> t({1,4,1}) ;  // values if T is in center

                kmermap = new ngsai::KmerMap(3) ;
                kmermap->insert(ngsai::KmerData(std::string("AAA"), a, a)) ;
                kmermap->insert(ngsai::KmerData(std::string("AAC"), a, a)) ;
                kmermap->insert(ngsai::KmerData(std::string("AAG"), a, a)) ;
                kmermap->insert(ngsai::KmerData(std::string("AAT"), a, a)) ;
                kmermap->insert(ngsai::KmerData(std::string("ACA"), c, c)) ;
                kmermap->insert(ngsai::KmerData(std::string("ACC"), c, c)) ;
                kmermap->insert(ngsai::KmerData(std::string("ACG"), c, c)) ;
                kmermap->insert(ngsai::KmerData(std::string("ACT"), c, c)) ;
                kmermap->insert(ngsai::KmerData(std::string("AGA"), g, g)) ;
                kmermap->insert(ngsai::KmerData(std::string("AGC"), g, g)) ;
                kmermap->insert(ngsai::KmerData(std::string("AGG"), g, g)) ;
                kmermap->insert(ngsai::KmerData(std::string("AGT"), g, g)) ;
                kmermap->insert(ngsai::KmerData(std::string("ATA"), t, t)) ;
                kmermap->insert(ngsai::KmerData(std::string("ATC"), t, t)) ;
                kmermap->insert(ngsai::KmerData(std::string("ATG"), t, t)) ;
                kmermap->insert(ngsai::KmerData(std::string("ATT"), t, t)) ;
                kmermap->insert(ngsai::KmerData(std::string("CAA"), a, a)) ;
                kmermap->insert(ngsai::KmerData(std::string("CAC"), a, a)) ;
                kmermap->insert(ngsai::KmerData(std::string("CAG"), a, a)) ;
                kmermap->insert(ngsai::KmerData(std::string("CAT"), a, a)) ;
                kmermap->insert(ngsai::KmerData(std::string("CCA"), c, c)) ;
                kmermap->insert(ngsai::KmerData(std::string("CCC"), c, c)) ;
                kmermap->insert(ngsai::KmerData(std::string("CCG"), c, c)) ;
                kmermap->insert(ngsai::KmerData(std::string("CCT"), c, c)) ;
                kmermap->insert(ngsai::KmerData(std::string("CGA"), g, g)) ;
                kmermap->insert(ngsai::KmerData(std::string("CGC"), g, g)) ;
                kmermap->insert(ngsai::KmerData(std::string("CGG"), g, g)) ;
                kmermap->insert(ngsai::KmerData(std::string("CGT"), g, g)) ;
                kmermap->insert(ngsai::KmerData(std::string("CTA"), t, t)) ;
                kmermap->insert(ngsai::KmerData(std::string("CTC"), t, t)) ;
                kmermap->insert(ngsai::KmerData(std::string("CTG"), t, t)) ;
                kmermap->insert(ngsai::KmerData(std::string("CTT"), t, t)) ;
                kmermap->insert(ngsai::KmerData(std::string("GAA"), a, a)) ;
                kmermap->insert(ngsai::KmerData(std::string("GAC"), a, a)) ;
                kmermap->insert(ngsai::KmerData(std::string("GAG"), a, a)) ;
                kmermap->insert(ngsai::KmerData(std::string("GAT"), a, a)) ;
                kmermap->insert(ngsai::KmerData(std::string("GCA"), c, c)) ;
                kmermap->insert(ngsai::KmerData(std::string("GCC"), c, c)) ;
                kmermap->insert(ngsai::KmerData(std::string("GCG"), c, c)) ;
                kmermap->insert(ngsai::KmerData(std::string("GCT"), c, c)) ;
                kmermap->insert(ngsai::KmerData(std::string("GGA"), g, g)) ;
                kmermap->insert(ngsai::KmerData(std::string("GGC"), g, g)) ;
                kmermap->insert(ngsai::KmerData(std::string("GGG"), g, g)) ;
                kmermap->insert(ngsai::KmerData(std::string("GGT"), g, g)) ;
                kmermap->insert(ngsai::KmerData(std::string("GTA"), t, t)) ;
                kmermap->insert(ngsai::KmerData(std::string("GTC"), t, t)) ;
                kmermap->insert(ngsai::KmerData(std::string("GTG"), t, t)) ;
                kmermap->insert(ngsai::KmerData(std::string("GTT"), t, t)) ;
                kmermap->insert(ngsai::KmerData(std::string("TAA"), a, a)) ;
                kmermap->insert(ngsai::KmerData(std::string("TAC"), a, a)) ;
                kmermap->insert(ngsai::KmerData(std::string("TAG"), a, a)) ;
                kmermap->insert(ngsai::KmerData(std::string("TAT"), a, a)) ;
                kmermap->insert(ngsai::KmerData(std::string("TCA"), c, c)) ;
                kmermap->insert(ngsai::KmerData(std::string("TCC"), c, c)) ;
                kmermap->insert(ngsai::KmerData(std::string("TCG"), c, c)) ;
                kmermap->insert(ngsai::KmerData(std::string("TCT"), c, c)) ;
                kmermap->insert(ngsai::KmerData(std::string("TGA"), g, g)) ;
                kmermap->insert(ngsai::KmerData(std::string("TGC"), g, g)) ;
                kmermap->insert(ngsai::KmerData(std::string("TGG"), g, g)) ;
                kmermap->insert(ngsai::KmerData(std::string("TGT"), g, g)) ;
                kmermap->insert(ngsai::KmerData(std::string("TTA"), t, t)) ;
                kmermap->insert(ngsai::KmerData(std::string("TTC"), t, t)) ;
                kmermap->insert(ngsai::KmerData(std::string("TTG"), t, t)) ;
                kmermap->insert(ngsai::KmerData(std::string("TTT"), t, t)) ;
            }
        }

        static void
        TearDownTestSuite()
        {   if(kmermap != nullptr)
            {   delete kmermap ; 
                kmermap = nullptr ;
            }
        }
}  ;


// initialise pointer
ngsai::KmerMap* NormalizedKineticModelTest::kmermap = nullptr ;


/*!
 * \brief Tests that constructor works properly.
 */
TEST_F(NormalizedKineticModelTest, constructor)
{
    NormalizedKineticModelToy model(*kmermap) ;

    ASSERT_EQ(model.getHistogramsIPD().size(), 0) ;
    ASSERT_EQ(model.getHistogramsPWD().size(), 0) ;
    ASSERT_EQ(kmerMapEqual(model.getBckgModel(), *kmermap), true) ;
    ASSERT_EQ(model.isInit(), false) ;
    ASSERT_EQ(model.isDensity(), false) ;
    ASSERT_EQ(model.isLog(), false) ;
}


/*!
 * \brief Tests that getKineticSignalRequiredSize() method works properly.
 */
TEST_F(NormalizedKineticModelTest, getKineticSignalRequiredSize)
{   
    size_t size(5) ;
    double min(-2.) ;
    double max(2.) ;
    size_t n_bin(4) ;
    double pseudocounts(10.) ;

    std::vector<size_t> kmer_sizes({3,5,7}) ;
    for(const auto& kmer_size : kmer_sizes)
    {   // size of kinetic signal must be longer than model, normalization 
        // trims its ends by kmer_size/2
        size_t size_exp = size + 2*(kmer_size / 2) ;
        
        
        ngsai::KmerMap map(kmer_size) ;
        NormalizedKineticModelToy model(map) ;
        
        // uninitialised
        ASSERT_EQ(model.getKineticSignalRequiredSize(), 0) ;

        // initialised
        model.setParameters(size,
                            min,
                            max,
                            n_bin,
                            pseudocounts) ;
        ASSERT_EQ(model.getKineticSignalRequiredSize(), size_exp) ;

    }
}


/*!
 * \brief Tests that the add() method with KineticSignal works properly.
 */
TEST_F(NormalizedKineticModelTest, add_kineticsignal)
{
    double pseudocounts(2.) ;

    double inf = std::numeric_limits<double>::infinity() ;

    // some kinetic stretches
    std::vector<double> values_1(7, -inf) ;
    std::vector<double> values_2(7, -2.) ;
    std::vector<double> values_3(7, -1.) ;
    std::vector<double> values_4(7, 0.) ;
    std::vector<double> values_5(7, 1.) ;
    std::vector<double> values_6(7, 2.) ;
    std::vector<double> values_7(7, inf) ;

    // some underlying sequences
    std::string seq_a("AAAAAAA") ;
    std::string seq_c("CCCCCCC") ;
    std::string seq_g("GGGGGGG") ;
    std::string seq_t("TTTTTTT") ;

    // ngsai::KineticSignal k1_a, k1_c, k1_g, k1_t ;
    // ngsai::KineticSignal k2_a, k2_c, k2_g, k2_t ;
    // ngsai::KineticSignal k3_a, k3_c, k3_g, k3_t ;
    // ngsai::KineticSignal k4_a, k4_c, k4_g, k4_t ;
    // ngsai::KineticSignal k5_a, k5_c, k5_g, k5_t ;
    // ngsai::KineticSignal k6_a, k6_c, k6_g, k6_t ;
    // ngsai::KineticSignal k7_a, k7_c, k7_g, k7_t ;

    ngsai::KineticSignal k1_a(seq_a, seq_a, values_1, values_1, values_1, values_1) ;
    ngsai::KineticSignal k1_c(seq_c, seq_c, values_1, values_1, values_1, values_1) ;
    ngsai::KineticSignal k1_g(seq_g, seq_g, values_1, values_1, values_1, values_1) ;
    ngsai::KineticSignal k1_t(seq_t, seq_t, values_1, values_1, values_1, values_1) ;
    ngsai::KineticSignal k2_a(seq_a, seq_a, values_2, values_2, values_2, values_2) ;
    ngsai::KineticSignal k2_c(seq_c, seq_c, values_2, values_2, values_2, values_2) ;
    ngsai::KineticSignal k2_g(seq_g, seq_g, values_2, values_2, values_2, values_2) ;
    ngsai::KineticSignal k2_t(seq_t, seq_t, values_2, values_2, values_2, values_2) ;
    ngsai::KineticSignal k3_a(seq_a, seq_a, values_3, values_3, values_3, values_3) ;
    ngsai::KineticSignal k3_c(seq_c, seq_c, values_3, values_3, values_3, values_3) ;
    ngsai::KineticSignal k3_g(seq_g, seq_g, values_3, values_3, values_3, values_3) ;
    ngsai::KineticSignal k3_t(seq_t, seq_t, values_3, values_3, values_3, values_3) ;
    ngsai::KineticSignal k4_a(seq_a, seq_a, values_4, values_4, values_4, values_4) ;
    ngsai::KineticSignal k4_c(seq_c, seq_c, values_4, values_4, values_4, values_4) ;
    ngsai::KineticSignal k4_g(seq_g, seq_g, values_4, values_4, values_4, values_4) ;
    ngsai::KineticSignal k4_t(seq_t, seq_t, values_4, values_4, values_4, values_4) ;
    ngsai::KineticSignal k5_a(seq_a, seq_a, values_5, values_5, values_5, values_5) ;
    ngsai::KineticSignal k5_c(seq_c, seq_c, values_5, values_5, values_5, values_5) ;
    ngsai::KineticSignal k5_g(seq_g, seq_g, values_5, values_5, values_5, values_5) ;
    ngsai::KineticSignal k5_t(seq_t, seq_t, values_5, values_5, values_5, values_5) ;
    ngsai::KineticSignal k6_a(seq_a, seq_a, values_6, values_6, values_6, values_6) ;
    ngsai::KineticSignal k6_c(seq_c, seq_c, values_6, values_6, values_6, values_6) ;
    ngsai::KineticSignal k6_g(seq_g, seq_g, values_6, values_6, values_6, values_6) ;
    ngsai::KineticSignal k6_t(seq_t, seq_t, values_6, values_6, values_6, values_6) ;
    ngsai::KineticSignal k7_a(seq_a, seq_a, values_7, values_7, values_7, values_7) ;
    ngsai::KineticSignal k7_c(seq_c, seq_c, values_7, values_7, values_7, values_7) ;
    ngsai::KineticSignal k7_g(seq_g, seq_g, values_7, values_7, values_7, values_7) ;
    ngsai::KineticSignal k7_t(seq_t, seq_t, values_7, values_7, values_7, values_7) ;


    NormalizedKineticModelToy model(*kmermap) ;
    
    // model is not initialised yet, must fail
    EXPECT_THROW(model.add(k2_a), std::runtime_error) ;
    // 4+2 bins (-inf,-2), [-2,-1), [-1,0), [0,1), [1,2), [2,inf)
    model.setParameters(5, -2, 2, 4, pseudocounts) ;
    // -inf / 1. -> -inf, add +1 count in (-inf, 2) -> curr. total 3
    model.add(k1_a) ;
    // -inf / 2. -> -inf, add +1 count in (-inf, 2) -> curr. total 4
    model.add(k1_c) ;
    // -inf / 3. -> -inf, add +1 count in (-inf, 2) -> curr. total 5
    model.add(k1_g) ;
    // -inf / 4. -> -inf, add +1 count in (-inf, 2) -> curr. total 6
    model.add(k1_t) ;
    // -2 / 1. -> -2.0, add +1 count in [-2,-1) -> curr. total 3
    model.add(k2_a) ;
    // -2 / 2. -> -1.0, add +1 count in [-1,0)  -> curr. total 3
    model.add(k2_c) ;
    // -2 / 3. -> -0.66, add +1 count in [-1,0) -> curr. total 4
    model.add(k2_g) ;
    // -2 / 4. -> -0.25, add +1 count in [-1,0) -> curr. total 5
    model.add(k2_t) ;
    // -1 / 1. -> -1.0, add +1 count in [-1,0) -> curr. total 6
    model.add(k3_a) ;
    // -1 / 2. -> -0.5, add +1 count in [-1,0) -> curr. total 7
    model.add(k3_c) ;
    // -1 / 3. -> -0.33, add +1 count in [-1,0) -> curr. total 8 
    model.add(k3_g) ;
    // -1 / 4. -> -0.25, add +1 count in [-1,0) -> curr. total 9 
    model.add(k3_t) ;
    // 0 / 1. -> 0, add +1 count in [0,1) -> curr. total 3
    model.add(k4_a) ;
    // 0 / 2. -> 0, add +1 count in [0,1) -> curr. total 4
    model.add(k4_c) ;
    // 0 / 3. -> 0, add +1 count in [0,1) -> curr. total 5 
    model.add(k4_g) ;
    // 0 / 4. -> 0, add +1 count in [0,1) -> curr. total 6 
    model.add(k4_t) ;
    // 1 / 1. -> 1.0, add +1 count in [1,2) -> curr. total 3
    model.add(k5_a) ;
    // 1 / 2. -> 0.5, add +1 count in [0,1) -> curr. total 7
    model.add(k5_c) ;
    // 1 / 3. -> 0.33, add +1 count in [0,1) -> curr. total 8 
    model.add(k5_g) ;
    // 1 / 4. -> 0.25, add +1 count in [0,1) -> curr. total 9 
    model.add(k5_t) ;
    // 2 / 1. -> 2.0, add +1 count in (2,inf) -> curr. total 3
    model.add(k6_a) ;
    // 2 / 2. -> 1.0, add +1 count in [1,2) -> curr. total 4
    model.add(k6_c) ;
    // 2 / 3. -> 0.66, add +1 count in [0,1) -> curr. total 10
    model.add(k6_g) ;
    // 2 / 4. -> 0.5, add +1 count in [0,1) -> curr. total 11
    model.add(k6_t) ;
    // inf / 1. -> inf, add +1 count in (2,inf) -> curr. total 4
    model.add(k7_a) ;
    // inf / 2. -> inf, add +1 count in (2,inf) -> curr. total 5
    model.add(k7_c) ;
    // inf / 3. -> inf, add +1 count in (2,inf) -> curr. total 6
    model.add(k7_g) ;
    // nf / 4. -> inf, add +1 count in (2,inf) -> curr. total 7
    model.add(k7_t) ;
    // if you don believe it, uncomment :-)
    // std::cerr << model << std::endl ;

    // expected counts in each IPD bin
    std::vector<double> ipd_count_expected({6,3,9,11,4,7}) ;
    auto& histograms_ipd = model.getHistogramsIPD() ;
    for(auto& histogram : histograms_ipd)
    {   size_t i=0 ;
        for(const auto& x : bh::indexed(histogram, bh::coverage::all))
        {   ASSERT_EQ(*x, ipd_count_expected[i]) ;
            i++ ;
        }
    }

    // expected counts in each PWD bin (same as IPD)
    std::vector<double> pwd_count_expected({6,3,9,11,4,7}) ;
    auto& histograms_pwd = model.getHistogramsPWD() ;
    for(auto& histogram : histograms_pwd)
    {   size_t i=0 ;
        for(const auto& x : bh::indexed(histogram, bh::coverage::all))
        {   ASSERT_EQ(*x, pwd_count_expected[i]) ;
            i++ ;
        }
    }

    std::string s_err1("A") ;
    std::string s_err2("AA") ;
    std::string s_err3("AAA") ;
    std::string s_err4("AAAA") ;
    std::string s_err6("AAAAAA") ;
    std::vector<double> v_err1({1.}) ;
    std::vector<double> v_err2({1.,2.}) ;
    std::vector<double> v_err3({1.,2.,3.}) ;
    std::vector<double> v_err4({1.,2.,3.,4.}) ;
    std::vector<double> v_err6({1.,2.,3.,4.,5.,6.}) ;
    ngsai::KineticSignal k_err0 ;
    ngsai::KineticSignal k_err1(s_err1, s_err1, v_err1, v_err1, v_err1, v_err1) ;
    ngsai::KineticSignal k_err2(s_err2, s_err2, v_err2, v_err2, v_err2, v_err2) ;
    ngsai::KineticSignal k_err3(s_err3, s_err3, v_err3, v_err3, v_err3, v_err3) ;
    ngsai::KineticSignal k_err4(s_err4, s_err4, v_err4, v_err4, v_err4, v_err4) ;
    ngsai::KineticSignal k_err6(s_err6, s_err6, v_err6, v_err6, v_err6, v_err6) ;

    // must fail because not same size as model
    EXPECT_THROW(model.add(k_err0), std::invalid_argument) ; // too short
    EXPECT_THROW(model.add(k_err1), std::invalid_argument) ; // too short
    EXPECT_THROW(model.add(k_err2), std::invalid_argument) ; // too short
    EXPECT_THROW(model.add(k_err3), std::invalid_argument) ; // too short
    EXPECT_THROW(model.add(k_err4), std::invalid_argument) ; // too short
    EXPECT_THROW(model.add(k_err6), std::invalid_argument) ; // too long
}


/*!
 * \brief Tests that the copy constructor works properly.
 */
TEST_F(NormalizedKineticModelTest, constructor_copy)
{
    // copy empty model
    NormalizedKineticModelToy model1(*kmermap) ;
    NormalizedKineticModelToy model2(*kmermap) ;
    ASSERT_EQ(model1.getHistogramsIPD().size(), 
              model2.getHistogramsIPD().size()) ;
    ASSERT_EQ(model1.getHistogramsPWD().size(), 
              model2.getHistogramsPWD().size()) ;
    ASSERT_EQ(kmerMapEqual(model1.getBckgModel(), model2.getBckgModel()),
              true) ;
    ASSERT_EQ(model1.isInit(),    model2.isInit()) ;
    ASSERT_EQ(model1.isDensity(), model2.isDensity()) ;
    ASSERT_EQ(model1.isLog(),     model2.isLog()) ;

    // copy model with data
    // [-inf,-1), [-1,0),  [0,1), [1,inf)
    ngsai::KineticSignal k(std::string("CCC"),
                           std::string("CCC"),
                           std::vector<double>({1., 1., 1.}),
                           std::vector<double>({1., 1., 1.}),
                           std::vector<double>({1., 1., 1.}),
                           std::vector<double>({1., 1., 1.})) ;
    model1.setParameters(1, -1, 1, 2, 1) ;
    model1.add(k) ; // 1. / 2. -> 0.5, add +1 to bin [0,1)

    NormalizedKineticModelToy model3(model1) ;

    ASSERT_EQ(model1.getHistogramsIPD().size(), 
              model3.getHistogramsIPD().size()) ;
    ASSERT_EQ(model1.getHistogramsPWD().size(), 
              model3.getHistogramsPWD().size()) ;
    ASSERT_EQ(kmerMapEqual(model1.getBckgModel(), model3.getBckgModel()),
              true) ;
    ASSERT_EQ(model1.isInit(),    model3.isInit()) ;
    ASSERT_EQ(model1.isDensity(), model3.isDensity()) ;
    ASSERT_EQ(model1.isLog(),     model3.isLog()) ;
    auto& histograms1_ipd = model1.getHistogramsIPD() ;
    auto& histograms1_pwd = model1.getHistogramsPWD() ;
    auto& histograms3_ipd = model3.getHistogramsIPD() ;
    auto& histograms3_pwd = model3.getHistogramsPWD() ;
    EXPECT_THAT(histograms1_ipd, ContainerEq(histograms3_ipd)) ;
    EXPECT_THAT(histograms1_pwd, ContainerEq(histograms3_pwd)) ;
}


/*!
 * \brief Tests that the move constructor works properly.
 */
TEST_F(NormalizedKineticModelTest, constructor_move)
{
    // assign empty model
    NormalizedKineticModelToy model1(*kmermap) ;
    NormalizedKineticModelToy model2(model1) ;
    NormalizedKineticModelToy model3(std::move(model2)) ;
    ASSERT_EQ(model1.getHistogramsIPD().size(), 
              model3.getHistogramsPWD().size()) ;
    ASSERT_EQ(kmerMapEqual(model1.getBckgModel(), model3.getBckgModel()),
              true) ;
    ASSERT_EQ(model1.isInit(),    model3.isInit()) ;
    ASSERT_EQ(model1.isDensity(), model3.isDensity()) ;
    ASSERT_EQ(model1.isLog(),     model3.isLog()) ;

    // copy model with data
    // [-inf,-1), [-1,0),  [0,1), [1,inf)
    ngsai::KineticSignal k(std::string("CCC"),
                           std::string("CCC"),
                           std::vector<double>({1., 1., 1.}),
                           std::vector<double>({1., 1., 1.}),
                           std::vector<double>({1., 1., 1.}),
                           std::vector<double>({1., 1., 1.})) ;
    model1.setParameters(1, -1, 1, 2, 1) ;
    model1.add(k) ; // 1. / 2. -> 0.5, add +1 to bin [0,1)

    NormalizedKineticModelToy model4(model1) ;
    NormalizedKineticModelToy model5(std::move(model4)) ;

    ASSERT_EQ(model1.getHistogramsIPD().size(),
              model5.getHistogramsIPD().size()) ;
    ASSERT_EQ(model1.getHistogramsPWD().size(),
              model5.getHistogramsPWD().size()) ;
    ASSERT_EQ(kmerMapEqual(model1.getBckgModel(), model5.getBckgModel()),
              true) ;
    ASSERT_EQ(model1.isInit(),    model5.isInit()) ;
    ASSERT_EQ(model1.isDensity(), model5.isDensity()) ;
    ASSERT_EQ(model1.isLog(),     model5.isLog()) ;
    auto& histograms1_ipd = model1.getHistogramsIPD() ;
    auto& histograms1_pwd = model1.getHistogramsPWD() ;
    auto& histograms5_ipd = model5.getHistogramsIPD() ;
    auto& histograms5_pwd = model5.getHistogramsPWD() ;
    EXPECT_THAT(histograms1_ipd, ContainerEq(histograms5_ipd)) ;
    EXPECT_THAT(histograms1_pwd, ContainerEq(histograms5_pwd)) ;
}


/*!
 * \brief Tests that the assignment operator works properly.
 */
TEST_F(NormalizedKineticModelTest, assignment_operator)
{
    // assign empty model
    NormalizedKineticModelToy model1(*kmermap) ;
    NormalizedKineticModelToy model2 = model1 ;
    ASSERT_EQ(model1.getHistogramsIPD().size(),
              model2.getHistogramsIPD().size()) ;
    ASSERT_EQ(model1.getHistogramsPWD().size(),
              model2.getHistogramsPWD().size()) ;
    ASSERT_EQ(kmerMapEqual(model1.getBckgModel(), model2.getBckgModel()),
              true) ;
    ASSERT_EQ(model1.isInit(),    model2.isInit()) ;
    ASSERT_EQ(model1.isDensity(), model2.isDensity()) ;
    ASSERT_EQ(model1.isLog(),     model2.isLog()) ;

    // assign model with data
    // [-inf,-1), [-1,0),  [0,1), [1,inf)
    ngsai::KineticSignal k(std::string("CCC"),
                           std::string("CCC"),
                           std::vector<double>({1., 1., 1.}),
                           std::vector<double>({1., 1., 1.}),
                           std::vector<double>({1., 1., 1.}),
                           std::vector<double>({1., 1., 1.})) ;
    model1.setParameters(1, -1, 1, 2, 1) ;
    model1.add(k) ; // 1. / 2. -> 0.5, add +1 to bin [0,1)

    model2 = model1 ;

    ASSERT_EQ(model1.getHistogramsIPD().size(),
              model2.getHistogramsIPD().size()) ;
    ASSERT_EQ(model1.getHistogramsPWD().size(),
              model2.getHistogramsPWD().size()) ;
    ASSERT_EQ(kmerMapEqual(model1.getBckgModel(), model2.getBckgModel()),
              true) ;
    ASSERT_EQ(model1.isInit(),    model2.isInit()) ;
    ASSERT_EQ(model1.isDensity(), model2.isDensity()) ;
    ASSERT_EQ(model1.isLog(),     model2.isLog()) ;
    auto& histograms1_ipd = model1.getHistogramsIPD() ;
    auto& histograms2_ipd = model2.getHistogramsIPD() ;
    auto& histograms1_pwd = model1.getHistogramsPWD() ;
    auto& histograms2_pwd = model2.getHistogramsPWD() ;
    EXPECT_THAT(histograms1_ipd, ContainerEq(histograms2_ipd)) ;
    EXPECT_THAT(histograms1_pwd, ContainerEq(histograms2_pwd)) ;
}


/*!
 * \brief Tests that the move assignment operator works properly.
 */
TEST_F(NormalizedKineticModelTest, assignment_move_operator)
{
    // assign empty model
    NormalizedKineticModelToy model1(*kmermap) ;
    NormalizedKineticModelToy model2 = model1 ;
    NormalizedKineticModelToy model3(std::move(model2)) ;
    ASSERT_EQ(model1.getHistogramsIPD().size(), 
              model3.getHistogramsIPD().size()) ;
    ASSERT_EQ(model1.getHistogramsPWD().size(), 
              model3.getHistogramsPWD().size()) ;
    ASSERT_EQ(kmerMapEqual(model1.getBckgModel(), model3.getBckgModel()),
              true) ;
    ASSERT_EQ(model1.isInit(),    model3.isInit()) ;
    ASSERT_EQ(model1.isDensity(), model3.isDensity()) ;
    ASSERT_EQ(model1.isLog(),     model3.isLog()) ;

    // assign model with data
    // [-inf,-1), [-1,0),  [0,1), [1,inf)
    ngsai::KineticSignal k(std::string("CCC"),
                           std::string("CCC"),
                           std::vector<double>({1., 1., 1.}),
                           std::vector<double>({1., 1., 1.}),
                           std::vector<double>({1., 1., 1.}),
                           std::vector<double>({1., 1., 1.})) ;
    model1.setParameters(1, -1, 1, 2, 1) ;
    model1.add(k) ; // 1. / 2. -> 0.5, add +1 to bin [0,1)

    model3 = model1 ;

    NormalizedKineticModelToy model4(std::move(model3)) ;

    ASSERT_EQ(model1.getHistogramsIPD().size(), 
              model4.getHistogramsIPD().size()) ;
    ASSERT_EQ(model1.getHistogramsPWD().size(), 
              model4.getHistogramsPWD().size()) ;
    ASSERT_EQ(kmerMapEqual(model1.getBckgModel(), model4.getBckgModel()),
              true) ;
    ASSERT_EQ(model1.isInit(),    model4.isInit()) ;
    ASSERT_EQ(model1.isDensity(), model4.isDensity()) ;
    ASSERT_EQ(model1.isLog(),     model4.isLog()) ;
    auto& histograms1_ipd = model1.getHistogramsIPD() ;
    auto& histograms1_pwd = model1.getHistogramsPWD() ;
    auto& histograms4_ipd = model4.getHistogramsIPD() ;
    auto& histograms4_pwd = model4.getHistogramsPWD() ;
    EXPECT_THAT(histograms1_ipd, ContainerEq(histograms4_ipd)) ;
    EXPECT_THAT(histograms1_pwd, ContainerEq(histograms4_pwd)) ;
}


/*!
 * \brief Tests that the copy() method works properly.
 */
TEST_F(NormalizedKineticModelTest, copy)
{
    // copy empty model
    NormalizedKineticModelToy model1(*kmermap) ;
    NormalizedKineticModelToy* copy(nullptr) ;
    copy = model1.copy() ;
    ASSERT_EQ(model1.getHistogramsIPD().size(), 
              copy->getHistogramsIPD().size()) ;
    ASSERT_EQ(model1.getHistogramsPWD().size(), 
              copy->getHistogramsPWD().size()) ;
    ASSERT_EQ(kmerMapEqual(model1.getBckgModel(), copy->getBckgModel()),
              true) ;
    ASSERT_EQ(model1.isInit(),    copy->isInit()) ;
    ASSERT_EQ(model1.isDensity(), copy->isDensity()) ;
    ASSERT_EQ(model1.isLog(),     copy->isLog()) ;
    delete copy ;
    copy = nullptr ;

    // assign model with data
    // [-inf,-1), [-1,0),  [0,1), [1,inf)
    ngsai::KineticSignal k(std::string("CCC"),
                           std::string("CCC"),
                           std::vector<double>({1., 1., 1.}),
                           std::vector<double>({1., 1., 1.}),
                           std::vector<double>({1., 1., 1.}),
                           std::vector<double>({1., 1., 1.})) ;
    model1.setParameters(1, -1, 1, 2, 1) ;
    model1.add(k) ; // 1. / 2. -> 0.5, add +1 to bin [0,1)

    copy = model1.copy() ;

    ASSERT_EQ(model1.getHistogramsIPD().size(), 
              copy->getHistogramsIPD().size()) ;
    ASSERT_EQ(model1.getHistogramsPWD().size(), 
              copy->getHistogramsPWD().size()) ;
    ASSERT_EQ(kmerMapEqual(model1.getBckgModel(), copy->getBckgModel()),
              true) ;
    ASSERT_EQ(model1.isInit(),    copy->isInit()) ;
    ASSERT_EQ(model1.isDensity(), copy->isDensity()) ;
    ASSERT_EQ(model1.isLog(),     copy->isLog()) ;
    auto& histograms1_ipd = model1.getHistogramsIPD() ;
    auto& histograms1_pwd = model1.getHistogramsPWD() ;
    auto& histogramsc_ipd = copy->getHistogramsIPD() ;
    auto& histogramsc_pwd = copy->getHistogramsPWD() ;
    EXPECT_THAT(histograms1_ipd, ContainerEq(histogramsc_ipd)) ;
    EXPECT_THAT(histograms1_pwd, ContainerEq(histogramsc_pwd)) ;

    delete copy ;
    copy = nullptr ;
}


/*!
 * \brief Tests that the logLikelihood() method works properly.
 */
TEST_F(NormalizedKineticModelTest, logLikelihood)
{   
    // 5 values -> 3 values normalized with 3bp normalization model
    std::vector<double> k_1({0., 1., 2., 3., 0.}) ;  // train model

    // some underlying sequences
    std::string seq_a("AAAAA") ;
    std::string seq_c("CCCCC") ;
    std::string seq_g("GGGGG") ;
    std::string seq_t("TTTTT") ;

    ngsai::KineticSignal s1_a(seq_a, seq_a, k_1, k_1, k_1, k_1) ;
    ngsai::KineticSignal s1_c(seq_c, seq_c, k_1, k_1, k_1, k_1) ;
    ngsai::KineticSignal s1_g(seq_g, seq_g, k_1, k_1, k_1, k_1) ;
    ngsai::KineticSignal s1_t(seq_t, seq_t, k_1, k_1, k_1, k_1) ;

    NormalizedKineticModelToy model(*kmermap) ;

    // not initialized
    EXPECT_THROW(model.logLikelihood(s1_a), std::runtime_error) ;

    // 3+2 bins (-inf,0), [0,1), [1,2), [2,3), [3,inf)
    model.setParameters(3, 0., 3., 3, 1.) ;

    // model contains
    // pos 0      1       4      2      1      1
    // pos 1      1       3      2      2      1
    // pos 2      1       2      3      1      2
    // bins   (-inf,0), [0,1), [1,2), [2,3), [3,inf)
    model.add(s1_a) ;
    model.add(s1_c) ;
    model.add(s1_g) ;
    model.add(s1_t) ;
    // std::cout << model << std::endl << std::endl ;;

    model.density() ;
    // model contains
    // pos 0     1/9     4/9    2/9    1/9    1/9
    // pos 1     1/9     3/9    2/9    2/9    1/9
    // pos 2     1/9     2/9    3/9    1/9    2/9
    // bins   (-inf,0), [0,1), [1,2), [2,3), [3,inf)
    // std::cout << model << std::endl << std::endl ;

    // not log density
    EXPECT_THROW(model.logLikelihood(s1_a), std::runtime_error) ;

    model.log() ;

    
    // computing the log likelihood of all these vectors will ensure to 
    // query each bin of all histograms 1x (if using a AAAAA sequence)
    std::vector<double> k_2({0.,  -1.,  -1.,  -1,  0.}) ;
    std::vector<double> k_3({0.,   0.,   0.,   0., 0.}) ;
    std::vector<double> k_4({0.,   1.,   1.,   1., 0.}) ;
    std::vector<double> k_5({0.,   2.,   2.,   2., 0.}) ;
    std::vector<double> k_6({0.,   3.,   3.,   3., 0.}) ;

    ngsai::KineticSignal s2_a(seq_a, seq_a, k_2, k_2, k_2, k_2) ;
    ngsai::KineticSignal s3_a(seq_a, seq_a, k_3, k_3, k_3, k_3) ;
    ngsai::KineticSignal s4_a(seq_a, seq_a, k_4, k_4, k_4, k_4) ;
    ngsai::KineticSignal s5_a(seq_a, seq_a, k_5, k_5, k_5, k_5) ;
    ngsai::KineticSignal s6_a(seq_a, seq_a, k_6, k_6, k_6, k_6) ;

    ASSERT_DOUBLE_EQ(model.logLikelihood(s2_a),   // IPD norm into -1, -1, -1
                                                  // PWD norm into -1, -1, -1
                     log(1./9.) + log(1./9.) + log(1./9.) +    // IPDR
                     log(1./9.) + log(1./9.) + log(1./9.)) ;   // PWDR

    ASSERT_DOUBLE_EQ(model.logLikelihood(s3_a),   // IPD norm into 0, 0, 0
                                                  // PWD norm into 0, 0, 0    
                     log(4./9.) + log(3./9.) + log(2./9.) +    // IPDR
                     log(4./9.) + log(3./9.) + log(2./9.)) ;   // PWDR
    
    ASSERT_DOUBLE_EQ(model.logLikelihood(s4_a),   // IPD norm into 1, 1, 1
                                                  // PWD norm into 1, 1, 1
                     log(2./9.) + log(2./9.) + log(3./9.) +    // IPDR
                     log(2./9.) + log(2./9.) + log(3./9.)) ;   // PWDR
    
    ASSERT_DOUBLE_EQ(model.logLikelihood(s5_a),   // IPD norm into 2, 2, 2
                                                  // PWD norm into 2, 2, 2
                     log(1./9.) + log(2./9.) + log(1./9.) +    // IPDR
                     log(1./9.) + log(2./9.) + log(1./9.)) ;   // PWDR
    
    ASSERT_DOUBLE_EQ(model.logLikelihood(s6_a),   // IPD norm into 3, 3, 3
                                                  // PWD norm into 3, 3, 3
                     log(1./9.) + log(1./9.) + log(2./9.) +    // IPDR
                     log(1./9.) + log(1./9.) + log(2./9.)) ;   // PWDR


    // must return nan because data have same size as model but missing data
    ngsai::KineticSignal s1_inc, s2_inc, s3_inc, s4_inc, s5_inc, s6_inc ;
    // misses seq fw
                                   s1_inc.setSequenceRv(seq_a) ;
    s1_inc.setIPDFw(k_2)         ; s1_inc.setIPDRv(k_2) ;
    s1_inc.setPWDFw(k_2)         ; s1_inc.setPWDRv(k_2) ;
    // misses seq rv
    s2_inc.setSequenceFw(seq_a) ;
    s2_inc.setIPDFw(k_2)        ; s2_inc.setIPDRv(k_2) ;
    s2_inc.setPWDFw(k_2)        ; s2_inc.setPWDRv(k_2) ;
    // misses IPD fw
    s3_inc.setSequenceFw(seq_a) ; s3_inc.setSequenceRv(seq_a) ;
                                  s3_inc.setIPDRv(k_2)    ;
    s3_inc.setPWDFw(k_2)        ; s3_inc.setPWDFw(k_2) ;
    // misses IPD rv
    s4_inc.setSequenceFw(seq_a) ; s4_inc.setSequenceRv(seq_a) ;
    s4_inc.setIPDFw(k_2)        ;
    s4_inc.setPWDFw(k_2)        ; s4_inc.setPWDRv(k_2) ;
    // misses PWD fw
    s5_inc.setSequenceFw(seq_a) ; s5_inc.setSequenceRv(seq_a) ;
    s5_inc.setIPDFw(k_2)        ; s5_inc.setIPDRv(k_2)    ;
                                  s5_inc.setPWDRv(k_2) ;
    // misses PWD rv
    s6_inc.setSequenceFw(seq_a) ; s6_inc.setSequenceRv(seq_a) ;
    s6_inc.setIPDFw(k_2)         ; s6_inc.setIPDRv(k_2)    ;
    s6_inc.setPWDFw(k_2)         ;
    // must return nan because incomplete
    ASSERT_EQ(std::isnan(model.logLikelihood(s1_inc)), true) ;
    ASSERT_EQ(std::isnan(model.logLikelihood(s2_inc)), true) ;
    ASSERT_EQ(std::isnan(model.logLikelihood(s3_inc)), true) ;
    ASSERT_EQ(std::isnan(model.logLikelihood(s4_inc)), true) ;
    ASSERT_EQ(std::isnan(model.logLikelihood(s5_inc)), true) ;
    ASSERT_EQ(std::isnan(model.logLikelihood(s6_inc)), true) ;


    // must fail, size of kinetic signal does not match size of model
    std::string seq_e1("A") ;
    std::string seq_e2("AA") ;
    std::string seq_e4("AAAA") ;
    std::string seq_e6("aAAAA") ;
    std::vector<double> val_e1({1}) ;
    std::vector<double> val_e2({1,1}) ;
    std::vector<double> val_e4({1,1,1,1}) ;
    std::vector<double> val_e6({1,1,1,1,1}) ;
    ngsai::KineticSignal s1_err(seq_e1,seq_e1,val_e1,val_e1,val_e1,val_e1) ;
    ngsai::KineticSignal s2_err(seq_e2,seq_e2,val_e2,val_e2,val_e2,val_e2) ;
    ngsai::KineticSignal s4_err(seq_e4,seq_e4,val_e4,val_e4,val_e4,val_e4) ;
    ngsai::KineticSignal s6_err(seq_e6,seq_e6,val_e6,val_e6,val_e6,val_e6) ;
    EXPECT_THROW(model.logLikelihood(s1_err), std::invalid_argument) ; // tpo short
    EXPECT_THROW(model.logLikelihood(s2_err), std::invalid_argument) ; // tpo short
    EXPECT_THROW(model.logLikelihood(s4_err), std::invalid_argument) ; // too short
    EXPECT_THROW(model.logLikelihood(s6_err), std::invalid_argument) ; // too long
}


/*!
 * \brief Tests that the add() with KineticModel methods 
 * works properly.
 */
TEST_F(NormalizedKineticModelTest, add_kineticmodel)
{
    double pseudocounts(2.) ;

    NormalizedKineticModelToy model1(*kmermap) ;
    NormalizedKineticModelToy model2(*kmermap) ;
    NormalizedKineticModelToy model3(*kmermap) ;
    NormalizedKineticModelToy model4(*kmermap) ;

    double inf = std::numeric_limits<double>::infinity() ;

    // some kinetic stretches
    std::vector<double> values_1(7, -inf) ;
    std::vector<double> values_2(7, -2.) ;
    std::vector<double> values_3(7, -1.) ;
    std::vector<double> values_4(7, 0.) ;
    std::vector<double> values_5(7, 1.) ;
    std::vector<double> values_6(7, 2.) ;
    std::vector<double> values_7(7, inf) ;

    // some underlying sequences
    std::string seq_a("AAAAAAA") ;
    std::string seq_c("CCCCCCC") ;
    std::string seq_g("GGGGGGG") ;
    std::string seq_t("TTTTTTT") ;

    ngsai::KineticSignal k1_a(seq_a, seq_a, values_1, values_1, values_1, values_1) ;
    ngsai::KineticSignal k1_c(seq_c, seq_c, values_1, values_1, values_1, values_1) ;
    ngsai::KineticSignal k1_g(seq_g, seq_g, values_1, values_1, values_1, values_1) ;
    ngsai::KineticSignal k1_t(seq_t, seq_t, values_1, values_1, values_1, values_1) ;
    ngsai::KineticSignal k2_a(seq_a, seq_a, values_2, values_2, values_2, values_2) ;
    ngsai::KineticSignal k2_c(seq_c, seq_c, values_2, values_2, values_2, values_2) ;
    ngsai::KineticSignal k2_g(seq_g, seq_g, values_2, values_2, values_2, values_2) ;
    ngsai::KineticSignal k2_t(seq_t, seq_t, values_2, values_2, values_2, values_2) ;
    ngsai::KineticSignal k3_a(seq_a, seq_a, values_3, values_3, values_3, values_3) ;
    ngsai::KineticSignal k3_c(seq_c, seq_c, values_3, values_3, values_3, values_3) ;
    ngsai::KineticSignal k3_g(seq_g, seq_g, values_3, values_3, values_3, values_3) ;
    ngsai::KineticSignal k3_t(seq_t, seq_t, values_3, values_3, values_3, values_3) ;
    ngsai::KineticSignal k4_a(seq_a, seq_a, values_4, values_4, values_4, values_4) ;
    ngsai::KineticSignal k4_c(seq_c, seq_c, values_4, values_4, values_4, values_4) ;
    ngsai::KineticSignal k4_g(seq_g, seq_g, values_4, values_4, values_4, values_4) ;
    ngsai::KineticSignal k4_t(seq_t, seq_t, values_4, values_4, values_4, values_4) ;
    ngsai::KineticSignal k5_a(seq_a, seq_a, values_5, values_5, values_5, values_5) ;
    ngsai::KineticSignal k5_c(seq_c, seq_c, values_5, values_5, values_5, values_5) ;
    ngsai::KineticSignal k5_g(seq_g, seq_g, values_5, values_5, values_5, values_5) ;
    ngsai::KineticSignal k5_t(seq_t, seq_t, values_5, values_5, values_5, values_5) ;
    ngsai::KineticSignal k6_a(seq_a, seq_a, values_6, values_6, values_6, values_6) ;
    ngsai::KineticSignal k6_c(seq_c, seq_c, values_6, values_6, values_6, values_6) ;
    ngsai::KineticSignal k6_g(seq_g, seq_g, values_6, values_6, values_6, values_6) ;
    ngsai::KineticSignal k6_t(seq_t, seq_t, values_6, values_6, values_6, values_6) ;
    ngsai::KineticSignal k7_a(seq_a, seq_a, values_7, values_7, values_7, values_7) ;
    ngsai::KineticSignal k7_c(seq_c, seq_c, values_7, values_7, values_7, values_7) ;
    ngsai::KineticSignal k7_g(seq_g, seq_g, values_7, values_7, values_7, values_7) ;
    ngsai::KineticSignal k7_t(seq_t, seq_t, values_7, values_7, values_7, values_7) ;

    // must fail, model1 not initialised
    EXPECT_THROW(model1.add(model4), std::invalid_argument) ;

    // 4+2 bins (-inf,-2), [-2,-1), [-1,0), [0,1), [1,2), [2,inf)
    model1.setParameters(5, -2, 2, 4, pseudocounts) ;
    model2.setParameters(5, -2, 2, 4, 0.) ;
    model3.setParameters(5, -2, 2, 4, 0.) ;

    // must fail, model4 not initialised
    EXPECT_THROW(model1.add(model4), std::invalid_argument) ;

    // 4+2 bins (-inf,-2), [-2,-1), [-1,0), [0,1), [1,2), [2,inf)
    model4.setParameters(6, -2, 2, 4, pseudocounts) ;

    // must fail, not same size
    EXPECT_THROW(model1.add(model4), std::invalid_argument) ;

    // 2+2 bins (-inf,-1), [-1,0), [0,1), [1,inf)
    model4.setParameters(5, -1, 1, 1, pseudocounts) ;

    // must fail, not same nb of bins
    EXPECT_THROW(model1.add(model4), std::runtime_error) ;

    
    // -inf / 1. -> -inf, add +1 count in (-inf, 2) -> curr. total 3
    model1.add(k1_a) ;
    // -inf / 2. -> -inf, add +1 count in (-inf, 2) -> curr. total 4
    model1.add(k1_c) ;
    // -inf / 3. -> -inf, add +1 count in (-inf, 2) -> curr. total 5
    model1.add(k1_g) ;
    // -inf / 4. -> -inf, add +1 count in (-inf, 2) -> curr. total 6
    model1.add(k1_t) ;
    // -2 / 1. -> -2.0, add +1 count in [-2,-1) -> curr. total 3
    model1.add(k2_a) ;
    // -2 / 2. -> -1.0, add +1 count in [-1,0)  -> curr. total 3
    model1.add(k2_c) ;
    // -2 / 3. -> -0.66, add +1 count in [-1,0) -> curr. total 4
    model1.add(k2_g) ;
    // -2 / 4. -> -0.25, add +1 count in [-1,0) -> curr. total 5
    model1.add(k2_t) ;
    // -1 / 1. -> -1.0, add +1 count in [-1,0) -> curr. total 6
    model1.add(k3_a) ;
    // -1 / 2. -> -0.5, add +1 count in [-1,0) -> curr. total 7
    model1.add(k3_c) ;
    // -1 / 3. -> -0.33, add +1 count in [-1,0) -> curr. total 8 
    model1.add(k3_g) ;
    // -1 / 4. -> -0.25, add +1 count in [-1,0) -> curr. total 9 
    model1.add(k3_t) ;
    // 0 / 1. -> 0, add +1 count in [0,1) -> curr. total 3
    model1.add(k4_a) ;
    // 0 / 2. -> 0, add +1 count in [0,1) -> curr. total 4
    model1.add(k4_c) ;
    // 0 / 3. -> 0, add +1 count in [0,1) -> curr. total 5 
    model1.add(k4_g) ;
    // 0 / 4. -> 0, add +1 count in [0,1) -> curr. total 6 
    model1.add(k4_t) ;
    // 1 / 1. -> 1.0, add +1 count in [1,2) -> curr. total 3
    model1.add(k5_a) ;
    // 1 / 2. -> 0.5, add +1 count in [0,1) -> curr. total 7
    model1.add(k5_c) ;
    // 1 / 3. -> 0.33, add +1 count in [0,1) -> curr. total 8 
    model1.add(k5_g) ;
    // 1 / 4. -> 0.25, add +1 count in [0,1) -> curr. total 9 
    model1.add(k5_t) ;
    // 2 / 1. -> 2.0, add +1 count in (2,inf) -> curr. total 3
    model1.add(k6_a) ;
    // 2 / 2. -> 1.0, add +1 count in [1,2) -> curr. total 4
    model1.add(k6_c) ;
    // 2 / 3. -> 0.66, add +1 count in [0,1) -> curr. total 10
    model1.add(k6_g) ;
    // 2 / 4. -> 0.5, add +1 count in [0,1) -> curr. total 11
    model1.add(k6_t) ;
    // inf / 1. -> inf, add +1 count in (2,inf) -> curr. total 4
    model1.add(k7_a) ;
    // inf / 2. -> inf, add +1 count in (2,inf) -> curr. total 5
    model1.add(k7_c) ;
    // inf / 3. -> inf, add +1 count in (2,inf) -> curr. total 6
    model1.add(k7_g) ;
    // nf / 4. -> inf, add +1 count in (2,inf) -> curr. total 7
    model1.add(k7_t) ;

    // expected counts in each IPD bin
    std::vector<double> count_expected({6,3,9,11,4,7}) ;

    // model2 += (model1 + model1) ;
    model2.add(model1) ;
    model2.add(model1) ;
    // model3 += (model1 + model1 + model1) ;
    model3.add(model1) ;
    model3.add(model1) ;
    model3.add(model1) ;

    // model2 is 2x model1
    auto& histograms_ipd = model2.getHistogramsIPD() ;
    for(auto& histogram : histograms_ipd)
    {   size_t i=0 ;
        for(const auto& x : bh::indexed(histogram, bh::coverage::all))
        {   ASSERT_DOUBLE_EQ(*x, 2*count_expected[i]) ;
            i++ ;
        }
    }
    auto& histograms_pwd = model2.getHistogramsPWD() ;
    for(auto& histogram : histograms_pwd)
    {   size_t i=0 ;
        for(const auto& x : bh::indexed(histogram, bh::coverage::all))
        {   ASSERT_DOUBLE_EQ(*x, 2*count_expected[i]) ;
            i++ ;
        }
    }
    ASSERT_EQ(model2.size(), 5) ;
    ASSERT_EQ(model2.getIsInit(), true) ;
    ASSERT_EQ(model2.getIsDensity(), false) ;
    ASSERT_EQ(model2.getIsLog(), false) ;
    ASSERT_EQ(kmerMapEqual(model2.getBckgModel(), *kmermap), true) ;
    // model3 is 3x model1
    histograms_ipd = model3.getHistogramsIPD() ;
    for(auto& histogram : histograms_ipd)
    {   size_t i=0 ;
        for(const auto& x : bh::indexed(histogram, bh::coverage::all))
        {   ASSERT_DOUBLE_EQ(*x, 3*count_expected[i]) ;
            i++ ;
        }
    }
    histograms_pwd = model3.getHistogramsPWD() ;
    for(auto& histogram : histograms_pwd)
    {   size_t i=0 ;
        for(const auto& x : bh::indexed(histogram, bh::coverage::all))
        {   ASSERT_DOUBLE_EQ(*x, 3*count_expected[i]) ;
            i++ ;
        }
    }
    ASSERT_EQ(model3.size(), 5) ;
    ASSERT_EQ(model3.getIsInit(), true) ;
    ASSERT_EQ(model3.getIsDensity(), false) ;
    ASSERT_EQ(model3.getIsLog(), false) ;
    ASSERT_EQ(kmerMapEqual(model3.getBckgModel(), *kmermap), true) ;


    
    // sum density models
    model1.density() ;

    NormalizedKineticModelToy model5(*kmermap) ;
    NormalizedKineticModelToy model6(*kmermap) ;
    model5.setParameters(5, -2, 2, 4, 0.) ;
    model6.setParameters(5, -2, 2, 4, 0.) ;

    // model5 += (model1 + model1) ;
    model5.add(model1) ;
    model5.add(model1) ;
    // model6 += (model1 + model1 + model1) ;
    model6.add(model1) ;
    model6.add(model1) ;
    model6.add(model1) ;

    // model5 is 2x model1
    histograms_ipd = model5.getHistogramsIPD() ;
    for(auto& histogram : histograms_ipd)
    {   size_t i=0 ;
        for(const auto& x : bh::indexed(histogram, bh::coverage::all))
        {   ASSERT_DOUBLE_EQ(*x, 2.*(count_expected[i] / 40.)) ;
            i++ ;
        }
    }
    histograms_pwd = model5.getHistogramsPWD() ;
    for(auto& histogram : histograms_pwd)
    {   size_t i=0 ;
        for(const auto& x : bh::indexed(histogram, bh::coverage::all))
        {   ASSERT_DOUBLE_EQ(*x, 2.*(count_expected[i] / 40.)) ;
            i++ ;
        }
    }
    ASSERT_EQ(model5.size(), 5) ;
    ASSERT_EQ(model5.getIsInit(), true) ;
    ASSERT_EQ(model5.getIsDensity(), false) ;
    ASSERT_EQ(model5.getIsLog(), false) ;
    ASSERT_EQ(kmerMapEqual(model5.getBckgModel(), *kmermap), true) ;
    // model6 is 3x model1
    histograms_ipd = model6.getHistogramsIPD() ;
    for(auto& histogram : histograms_ipd)
    {   size_t i=0 ;
        for(const auto& x : bh::indexed(histogram, bh::coverage::all))
        {   ASSERT_DOUBLE_EQ(*x, 3.*(count_expected[i] / 40.)) ;
            i++ ;
        }
    }
    histograms_pwd = model6.getHistogramsPWD() ;
    for(auto& histogram : histograms_pwd)
    {   size_t i=0 ;
        for(const auto& x : bh::indexed(histogram, bh::coverage::all))
        {   ASSERT_DOUBLE_EQ(*x, 3.*(count_expected[i] / 40.)) ;
            i++ ;
        }
    }
    ASSERT_EQ(model6.size(), 5) ;
    ASSERT_EQ(model6.getIsInit(), true) ;
    ASSERT_EQ(model6.getIsDensity(), false) ;
    ASSERT_EQ(model6.getIsLog(), false) ;
    ASSERT_EQ(kmerMapEqual(model6.getBckgModel(), *kmermap), true) ;


    // sum log density models
    model1.log() ;

    NormalizedKineticModelToy model7(*kmermap) ;
    NormalizedKineticModelToy model8(*kmermap) ;
    model7.setParameters(5, -2, 2, 4, 0.) ;
    model8.setParameters(5, -2, 2, 4, 0.) ;

    // model7 += (model1 + model1) ;
    model7.add(model1) ;
    model7.add(model1) ;
    // model8 += (model1 + model1 + model1) ;
    model8.add(model1) ;
    model8.add(model1) ;
    model8.add(model1) ;

    // model7 is 2x model1
    histograms_ipd = model7.getHistogramsIPD() ;
    for(auto& histogram : histograms_ipd)
    {   size_t i=0 ;
        for(const auto& x : bh::indexed(histogram, bh::coverage::all))
        {   ASSERT_DOUBLE_EQ(*x, 2.*log((count_expected[i] / 40.))) ;
            i++ ;
        }
    }
    histograms_pwd = model7.getHistogramsPWD() ;
    for(auto& histogram : histograms_pwd)
    {   size_t i=0 ;
        for(const auto& x : bh::indexed(histogram, bh::coverage::all))
        {   ASSERT_DOUBLE_EQ(*x, 2.*log((count_expected[i] / 40.))) ;
            i++ ;
        }
    }
    ASSERT_EQ(model7.size(), 5) ;
    ASSERT_EQ(model7.getIsInit(), true) ;
    ASSERT_EQ(model7.getIsDensity(), false) ;
    ASSERT_EQ(model7.getIsLog(), false) ;
    ASSERT_EQ(kmerMapEqual(model7.getBckgModel(), *kmermap), true) ;
    // model8 is 3x model1
    histograms_ipd = model8.getHistogramsIPD() ;
    for(auto& histogram : histograms_ipd)
    {   size_t i=0 ;
        for(const auto& x : bh::indexed(histogram, bh::coverage::all))
        {   ASSERT_DOUBLE_EQ(*x, 3.*log((count_expected[i] / 40.))) ;
            i++ ;
        }
    }
    histograms_pwd = model8.getHistogramsPWD() ;
    for(auto& histogram : histograms_pwd)
    {   size_t i=0 ;
        for(const auto& x : bh::indexed(histogram, bh::coverage::all))
        {   ASSERT_DOUBLE_EQ(*x, 3.*log((count_expected[i] / 40.))) ;
            i++ ;
        }
    }
    ASSERT_EQ(model8.size(), 5) ;
    ASSERT_EQ(model8.getIsInit(), true) ;
    ASSERT_EQ(model8.getIsDensity(), false) ;
    ASSERT_EQ(model8.getIsLog(), false) ;
    ASSERT_EQ(kmerMapEqual(model8.getBckgModel(), *kmermap), true) ;
}