#include <gtest/gtest.h>
#include <gmock/gmock.h>
#include <string>
#include <vector>
#include <set>
#include <stdexcept>
#include <limits>
#include <cmath>

#include <ngsaipp/epigenetics/DiPositionNormalizedKineticModel.hpp>
#include "utility.hpp"


using testing::ElementsAre ;
using testing::ContainerEq ;


/*
 * This file ontains tests for the ngsai::DiPositionNormalizedKineticModel class from 
 * src/epigenetics/DiPositionNormalizedKineticModel.cpp
 */



/*!
 * \brief A child class of RawKineticModel allowing to access the protected 
 * fields of DiPositionNormalizedKineticModel for testing purposes.
 */
class DiPositionNormalizedKineticModelToy : public ngsai::DiPositionNormalizedKineticModel
{
    public:

        DiPositionNormalizedKineticModelToy()
            : DiPositionNormalizedKineticModel()
        { ; }

        DiPositionNormalizedKineticModelToy(const ngsai::DiPositionNormalizedKineticModel& other)
            : DiPositionNormalizedKineticModel(other)
        { ; }

        DiPositionNormalizedKineticModelToy(ngsai::DiPositionNormalizedKineticModel&& other)
            : DiPositionNormalizedKineticModel(other)
        { ; }

        ~DiPositionNormalizedKineticModelToy()
        { ; }

        // ugly workaround to test RawKineticModel.copy() using a
        // RawKineticModelTest instance
        DiPositionNormalizedKineticModelToy* copy() const
        {   DiPositionNormalizedKineticModel* copy = DiPositionNormalizedKineticModel::copy() ;
            DiPositionNormalizedKineticModelToy* copy2 = new DiPositionNormalizedKineticModelToy(*copy) ;
            delete copy ;
            return copy2 ;
        }

        size_t getPairNb()
        {   return m_pair_nb ; }

        std::vector<std::pair<size_t,size_t>>& getPairIndex() 
        {   return m_pair_index ; }

        std::vector<std::vector<hist_2d_double>>& getHistogramsIPD()
        {   return m_histograms_ipd ; }

        std::vector<std::vector<hist_2d_double>>& getHistogramsPWD()
        {   return m_histograms_pwd ; }

        ngsai::KmerMap& getBckgModel()
        {   return m_bckg_model ; }

        bool& hasKmerMap()
        {   return m_has_kmermap ; }

        bool& getIsInit()
        {   return m_is_init ; }

        bool& getIsDensity()
        {   return m_is_density ; }

        bool& getIsLog()
        {   return m_is_log ; }

} ;


/*!
 * \brief Contains the necessary KmerMap that will be used by each and every 
 * PairWiseNormalizedKineticModelToy instance.
 * Rule is simple, the center position get an expected IPD/PDW of :
 * 1 if there is A in the center
 * 2 if there is C in the center
 * 3 if there is G in the center
 * 4 if there is T in the center
 */
class DiPositionNormalizedKineticModelTest : public ::testing::Test
{   
    protected:

        virtual void 
        SetUp() override
        { ; }

        virtual void 
        TearDown() override
        { ; }

    protected:

        // create a model with 3 position, 2 bins, xmin=0., xmax=2., 
        // pseudocounts=0.(-> with 2+2 bins on both axes)
        // bins are (-inf,0), [0,1), [0,2), [2,inf)
        // adding all these kinetics will make the model contain, for IPD and 
        // PWD, for all POSn vs POSn+1 histogram
        //          (-inf,0), [0,1), [0,2), [2,inf)
        // [2,inf)      5       5      5       5
        // [0,2)        5       1      1       1
        // [0,1)        5       1      1       1
        // (-inf,0)     9       1      1       1
        static std::vector<ngsai::KineticSignal> kinetics ;
        
        // all kinetics will get this sequence as fw and rv sequence
        static std::string sequence ;

        // this model will normalize kinetics using the following schema :
        // if central position is A -> expected signal is 1
        // if central position is C -> expected signal is 2
        // if central position is G -> expected signal is 3
        // if central position is T -> expected signal is 4
        static ngsai::KmerMap* kmermap ;

        static void 
        SetUpTestSuite()
        {   
            // setup kinetics
            if(kinetics.size() == 0)
            {   
                std::vector<std::vector<double>> data ;
                data.push_back(std::vector<double>({0.,-1.,-1.,-1.,0})) ;
                data.push_back(std::vector<double>({0., 0.,-1.,-1.,0})) ;
                data.push_back(std::vector<double>({0., 1.,-1.,-1.,0})) ;
                data.push_back(std::vector<double>({0., 2.,-1.,-1.,0})) ;
                data.push_back(std::vector<double>({0.,-1., 0.,-1.,0})) ;
                data.push_back(std::vector<double>({0., 0., 0.,-1.,0})) ;
                data.push_back(std::vector<double>({0., 1., 0.,-1.,0})) ;
                data.push_back(std::vector<double>({0., 2., 0.,-1.,0})) ;
                data.push_back(std::vector<double>({0.,-1., 1.,-1.,0})) ;
                data.push_back(std::vector<double>({0., 0., 1.,-1.,0})) ;
                data.push_back(std::vector<double>({0., 1., 1.,-1.,0})) ;
                data.push_back(std::vector<double>({0., 2., 1.,-1.,0})) ;
                data.push_back(std::vector<double>({0.,-1., 2.,-1.,0})) ;
                data.push_back(std::vector<double>({0., 0., 2.,-1.,0})) ;
                data.push_back(std::vector<double>({0., 1., 2.,-1.,0})) ;
                data.push_back(std::vector<double>({0., 2., 2.,-1.,0})) ;
                
                data.push_back(std::vector<double>({0,-1, -1.,-1.,0})) ;
                data.push_back(std::vector<double>({0, 0, -1.,-1.,0})) ;
                data.push_back(std::vector<double>({0, 1, -1.,-1.,0})) ;
                data.push_back(std::vector<double>({0, 2, -1.,-1.,0})) ;
                data.push_back(std::vector<double>({0,-1, -1., 0.,0})) ;
                data.push_back(std::vector<double>({0, 0, -1., 0.,0})) ;
                data.push_back(std::vector<double>({0, 1, -1., 0.,0})) ;
                data.push_back(std::vector<double>({0, 2, -1., 0.,0})) ;
                data.push_back(std::vector<double>({0,-1, -1., 1.,0})) ;
                data.push_back(std::vector<double>({0, 0, -1., 1.,0})) ;
                data.push_back(std::vector<double>({0, 1, -1., 1.,0})) ;
                data.push_back(std::vector<double>({0, 2, -1., 1.,0})) ;
                data.push_back(std::vector<double>({0,-1, -1., 2.,0})) ;
                data.push_back(std::vector<double>({0, 0, -1., 2.,0})) ;
                data.push_back(std::vector<double>({0, 1, -1., 2.,0})) ;
                data.push_back(std::vector<double>({0, 2, -1., 2.,0})) ;
                
                data.push_back(std::vector<double>({0,-1, -1.,-1.,0})) ;
                data.push_back(std::vector<double>({0,-1,  0.,-1.,0})) ;
                data.push_back(std::vector<double>({0,-1,  1.,-1.,0})) ;
                data.push_back(std::vector<double>({0,-1,  2.,-1.,0})) ;
                data.push_back(std::vector<double>({0,-1, -1., 0.,0})) ;
                data.push_back(std::vector<double>({0,-1,  0., 0.,0})) ;
                data.push_back(std::vector<double>({0,-1,  1., 0.,0})) ;
                data.push_back(std::vector<double>({0,-1,  2., 0.,0})) ;
                data.push_back(std::vector<double>({0,-1, -1., 1.,0})) ;
                data.push_back(std::vector<double>({0,-1,  0., 1.,0})) ;
                data.push_back(std::vector<double>({0,-1,  1., 1.,0})) ;
                data.push_back(std::vector<double>({0,-1,  2., 1.,0})) ;
                data.push_back(std::vector<double>({0,-1, -1., 2.,0})) ;
                data.push_back(std::vector<double>({0,-1,  0., 2.,0})) ;
                data.push_back(std::vector<double>({0,-1,  1., 2.,0})) ;
                data.push_back(std::vector<double>({0,-1,  2., 2.,0})) ;
                
                // give all kinetics a A in the middle in order to get an 
                // expected kinetic signal of 1 (normalization is N / 1 :-))
                // -> see kmer map below
                sequence = std::string("AAAAA") ;
                for(const auto& datum : data)
                {   kinetics.push_back(ngsai::KineticSignal(sequence, sequence,
                                                            datum,    datum,
                                                            datum,    datum)) ;
                }
            }

            // setup background model
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


// initialise static members
ngsai::KmerMap*                   DiPositionNormalizedKineticModelTest::kmermap = nullptr ;
std::string                       DiPositionNormalizedKineticModelTest::sequence = std::string() ;
std::vector<ngsai::KineticSignal> DiPositionNormalizedKineticModelTest::kinetics = std::vector<ngsai::KineticSignal>() ;


/*!
 * \brief Tests that the default constructor works properly.
 */
TEST_F(DiPositionNormalizedKineticModelTest, constructor)
{
    DiPositionNormalizedKineticModelToy model ;
    ASSERT_EQ(model.getPairNb(), 0) ;
    ASSERT_EQ(model.getPairIndex().size(), 0) ;
    ASSERT_EQ(model.getHistogramsIPD().size(), 0) ;
    ASSERT_EQ(model.getHistogramsPWD().size(), 0) ;
    ASSERT_EQ(model.isInit(), false) ;
    ASSERT_EQ(model.isDensity(), false) ;
    ASSERT_EQ(model.isLog(), false) ;
    ASSERT_EQ(model.hasKmerMap(), false) ;
    ASSERT_EQ(kmerMapEqual(model.getBckgModel(),
                          ngsai::KmerMap()),
              true) ;

}


/*!
 * \brief Tests that constructor works properly.
 */
TEST_F(DiPositionNormalizedKineticModelTest, constructor_kmermap)
{
    DiPositionNormalizedKineticModelToy model(*kmermap) ;
    ASSERT_EQ(model.getPairNb(), 0) ;
    ASSERT_EQ(model.getPairIndex().size(), 0) ;
    ASSERT_EQ(model.getHistogramsIPD().size(), 0) ;
    ASSERT_EQ(model.getHistogramsPWD().size(), 0) ;
    ASSERT_EQ(kmerMapEqual(model.getBckgModel(), *kmermap), true) ;
    ASSERT_EQ(model.isInit(), false) ;
    ASSERT_EQ(model.isDensity(), false) ;
    ASSERT_EQ(model.isLog(), false) ;
    ASSERT_EQ(model.hasKmerMap(), true) ;
}


/*!
 * \brief Tests that getKineticSignalRequiredSize() method works properly.
 */
TEST_F(DiPositionNormalizedKineticModelTest, getKineticSignalRequiredSize)
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
        DiPositionNormalizedKineticModelToy model(map) ;
        
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
 * \brief Tests that the setParameters() method works properly and that 
 * everything is initialised properly or throws the expected exceptions.
 */
TEST_F(DiPositionNormalizedKineticModelTest, setParameters)
{   
    size_t size(5) ;
    double min(-2.) ;
    double max(2.) ;
    size_t n_bin(4) ;
    double pseudocounts(10.) ;

    // it must work properly
    DiPositionNormalizedKineticModelToy model(*kmermap) ;
    model.setParameters(size,
                        min,
                        max,
                        n_bin,
                        pseudocounts) ;
    // expected pairs of positions to consider
    std::vector<std::pair<size_t,size_t>> pair_index_exp = 
            {std::pair<size_t,size_t>(0,1),
             std::pair<size_t,size_t>(1,2),
             std::pair<size_t,size_t>(2,3),
             std::pair<size_t,size_t>(3,4)} ;

    ASSERT_EQ(model.getPairNb(), pair_index_exp.size()) ;
    ASSERT_THAT(model.getPairIndex(), ContainerEq(pair_index_exp)) ;
    ASSERT_EQ(model.getHistogramsIPD().size(), size) ;
    ASSERT_EQ(model.getHistogramsPWD().size(), size) ;
    ASSERT_EQ(model.isInit(), true) ;
    ASSERT_EQ(model.isDensity(), false) ;
    ASSERT_EQ(model.isLog(), false) ;

    double inf = std::numeric_limits<double>::infinity() ;

    // expected bin boundaries on both axis (lower, upper)
    std::set<std::pair<double,double>> bins_exp ;
    bins_exp.insert(std::pair<double,double>(-inf,-2)) ;
    bins_exp.insert(std::pair<double,double>(-2,  -1)) ;
    bins_exp.insert(std::pair<double,double>(-1,  -0)) ;
    bins_exp.insert(std::pair<double,double>( 0,   1)) ;
    bins_exp.insert(std::pair<double,double>( 1,   2)) ;
    bins_exp.insert(std::pair<double,double>( 2,   inf)) ;
        
    const auto& index_pairs = model.getPairIndex() ;
    auto& histograms_ipd = model.getHistogramsIPD() ;
    for(const auto& index_pair : index_pairs)
    {   
        const size_t& pos1 = index_pair.first ;
        const size_t& pos2 = index_pair.second ;

        auto& histogram = histograms_ipd[pos1][pos2] ;
        // check bin content and boundaries
        std::set<std::pair<double,double>> bins_axis1, bins_axis2 ;
        for(const auto& x : bh::indexed(histogram, bh::coverage::all))
        {   // bin along 1st axis
            bins_axis1.insert(std::make_pair(x.bin(0).lower(),
                                             x.bin(0).upper())) ;
            // bin along 2nd axis
            bins_axis2.insert(std::make_pair(x.bin(1).lower(),
                                             x.bin(1).upper())) ;
            ASSERT_EQ(*x, pseudocounts) ;
        }
        EXPECT_THAT(bins_axis1, ContainerEq(bins_exp)) ;
        EXPECT_THAT(bins_axis2, ContainerEq(bins_exp)) ;
    }
    auto& histograms_pwd = model.getHistogramsPWD() ;
    for(const auto& index_pair : index_pairs)
    {   
        const size_t& pos1 = index_pair.first ;
        const size_t& pos2 = index_pair.second ;

        auto& histogram = histograms_pwd[pos1][pos2] ;
        // check bin content and boundaries
        std::set<std::pair<double,double>> bins_axis1, bins_axis2 ;
        for(const auto& x : bh::indexed(histogram, bh::coverage::all))
        {   // bin along 1st axis
            bins_axis1.insert(std::make_pair(x.bin(0).lower(),
                                             x.bin(0).upper())) ;
            // bin along 2nd axis
            bins_axis2.insert(std::make_pair(x.bin(1).lower(),
                                             x.bin(1).upper())) ;
            ASSERT_EQ(*x, pseudocounts) ;
        }
        EXPECT_THAT(bins_axis1, ContainerEq(bins_exp)) ;
        EXPECT_THAT(bins_axis2, ContainerEq(bins_exp)) ;
    }
    
    // it must fail, size < 2
    EXPECT_THROW(model.setParameters(0, min, max, n_bin ,pseudocounts),
                 std::invalid_argument) ;
    EXPECT_THROW(model.setParameters(1, min, max, n_bin ,pseudocounts),
                 std::invalid_argument) ;
    
    // it must fail, x_min > x_max
    EXPECT_THROW(model.setParameters(size, 10., 2., n_bin ,pseudocounts),
                 std::invalid_argument) ;
    
    // it must fail, x_min == x_max
    EXPECT_THROW(model.setParameters(size, 2., 2., n_bin ,pseudocounts),
                 std::invalid_argument) ;
    
    // it must fail, n_bins == 0
    EXPECT_THROW(model.setParameters(size, min, max, 0 ,pseudocounts),
                 std::invalid_argument) ;
}


/*!
 * \brief Tests that the copy constructor works properly.
 */
TEST_F(DiPositionNormalizedKineticModelTest, constructor_copy)
{
    // copy empty model
    DiPositionNormalizedKineticModelToy model1(*kmermap) ;
    DiPositionNormalizedKineticModelToy model2(model1) ;
    ASSERT_EQ(model1, model2) ;

    
    // 2+2 bins on both axes (-inf,0), [0,1), [0,2), [2,inf)
    model1.setParameters(3, 0., 2., 2, 0.) ;
    model1.add(kinetics[0]) ; // add +1 to bin (-inf,0),(-inf,0) at each pos
    // at this moment the model contains, for IPD and PWD, for all 
    // POSn vs POSn+1 histogram
    //          (-inf,0), [0,1), [0,2), [2,inf)
    // [2,inf)      0       0      0       0
    // [0,2)        0       0      0       0
    // [0,1)        0       0      0       0
    // (-inf,0)     1       0      0       0
    std::vector<std::vector<double>> m_expected({{1.,0.,0.,0.},
                                                 {0.,0.,0.,0.},
                                                 {0.,0.,0.,0.},
                                                 {0.,0.,0.,0.}}) ;
    // check histograms
    auto& index_pairs = model1.getPairIndex() ;
    auto& histograms1_ipd = model1.getHistogramsIPD() ;
    auto& histograms1_pwd = model1.getHistogramsPWD() ;
    for(const auto& index_pair : index_pairs)
    {   const size_t& pos1 = index_pair.first ;
        const size_t& pos2 = index_pair.second ;

        std::vector<std::vector<double>> m_ipd = 
                ngsai::hist2d_to_matrix(histograms1_ipd[pos1][pos2]) ;
        EXPECT_THAT(m_ipd, m_expected) ;

        std::vector<std::vector<double>> m_pwd = 
                ngsai::hist2d_to_matrix(histograms1_pwd[pos1][pos2]) ;
        EXPECT_THAT(m_pwd, m_expected) ;
    }

    // copy model with data
    DiPositionNormalizedKineticModelToy model3(model1) ;
    ASSERT_EQ(model1, model3) ;

}


/*!
 * \brief Tests that the move constructor works properly.
 */
TEST_F(DiPositionNormalizedKineticModelTest, constructor_move)
{
    // assign empty model
    DiPositionNormalizedKineticModelToy model1(*kmermap) ;
    DiPositionNormalizedKineticModelToy model2(model1) ;
    DiPositionNormalizedKineticModelToy model3(std::move(model2)) ;
    ASSERT_EQ(model1, model3) ;


    // 2+2 bins on both axes (-inf,0), [0,1), [0,2), [2,inf)
    model1.setParameters(3, 0., 2., 2, 0.) ;
    model1.add(kinetics[0]) ; // add +1 to bin (-inf,0),(-inf,0) at each pos
    // at this moment the model contains, for IPD and PWD, for all 
    // POSn vs POSn+1 histogram
    //          (-inf,0), [0,1), [0,2), [2,inf)
    // [2,inf)      0       0      0       0
    // [0,2)        0       0      0       0
    // [0,1)        0       0      0       0
    // (-inf,0)     1       0      0       0
    std::vector<std::vector<double>> m_expected({{1.,0.,0.,0.},
                                                 {0.,0.,0.,0.},
                                                 {0.,0.,0.,0.},
                                                 {0.,0.,0.,0.}}) ;
    // check histograms
    auto& index_pairs = model1.getPairIndex() ;
    auto& histograms1_ipd = model1.getHistogramsIPD() ;
    auto& histograms1_pwd = model1.getHistogramsPWD() ;
    for(const auto& index_pair : index_pairs)
    {   const size_t& pos1 = index_pair.first ;
        const size_t& pos2 = index_pair.second ;

        std::vector<std::vector<double>> m_ipd = 
                ngsai::hist2d_to_matrix(histograms1_ipd[pos1][pos2]) ;
        EXPECT_THAT(m_ipd, m_expected) ;

        std::vector<std::vector<double>> m_pwd = 
                ngsai::hist2d_to_matrix(histograms1_pwd[pos1][pos2]) ;
        EXPECT_THAT(m_pwd, m_expected) ;
    }

    // move construct
    DiPositionNormalizedKineticModelToy model4(model1) ;
    DiPositionNormalizedKineticModelToy model5(std::move(model4)) ;
    ASSERT_EQ(model1, model5) ;
}


/*!
 * \brief Tests that the assignment operator works properly.
 */
TEST_F(DiPositionNormalizedKineticModelTest, assignment_operator)
{
    // assign empty model
    DiPositionNormalizedKineticModelToy model1(*kmermap) ;
    DiPositionNormalizedKineticModelToy model2 = model1 ;
    ASSERT_EQ(model1, model2) ;

    // 2+2 bins on both axes (-inf,0), [0,1), [0,2), [2,inf)
    model1.setParameters(3, 0., 2., 2, 0.) ;
    model1.add(kinetics[0]) ; // add +1 to bin (-inf,0),(-inf,0) at each pos
    // at this moment the model contains, for IPD and PWD, for all 
    // POSn vs POSn+1 histogram
    //          (-inf,0), [0,1), [0,2), [2,inf)
    // [2,inf)      0       0      0       0
    // [0,2)        0       0      0       0
    // [0,1)        0       0      0       0
    // (-inf,0)     1       0      0       0
    std::vector<std::vector<double>> m_expected({{1.,0.,0.,0.},
                                                 {0.,0.,0.,0.},
                                                 {0.,0.,0.,0.},
                                                 {0.,0.,0.,0.}}) ;
    // check histograms
    auto& index_pairs = model1.getPairIndex() ;
    auto& histograms1_ipd = model1.getHistogramsIPD() ;
    auto& histograms1_pwd = model1.getHistogramsPWD() ;
    for(const auto& index_pair : index_pairs)
    {   const size_t& pos1 = index_pair.first ;
        const size_t& pos2 = index_pair.second ;

        std::vector<std::vector<double>> m_ipd = 
                ngsai::hist2d_to_matrix(histograms1_ipd[pos1][pos2]) ;
        EXPECT_THAT(m_ipd, m_expected) ;

        std::vector<std::vector<double>> m_pwd = 
                ngsai::hist2d_to_matrix(histograms1_pwd[pos1][pos2]) ;
        EXPECT_THAT(m_pwd, m_expected) ;
    }

    // assign
    model2 = model1 ;
    ASSERT_EQ(model1, model2) ;
}


/*!
 * \brief Tests that the move assignment operator works properly.
 */
TEST_F(DiPositionNormalizedKineticModelTest, assignment_move_operator)
{
    // assign empty model
    DiPositionNormalizedKineticModelToy model1(*kmermap) ;
    DiPositionNormalizedKineticModelToy model2 = model1 ;
    DiPositionNormalizedKineticModelToy model3(std::move(model2)) ;
    ASSERT_EQ(model1, model3) ;

    // 2+2 bins on both axes (-inf,0), [0,1), [0,2), [2,inf)
    model1.setParameters(3, 0., 2., 2, 0.) ;
    model1.add(kinetics[0]) ; // add +1 to bin (-inf,0),(-inf,0) at each pos
    // at this moment the model contains, for IPD and PWD, for all 
    // POSn vs POSn+1 histogram
    //          (-inf,0), [0,1), [0,2), [2,inf)
    // [2,inf)      0       0      0       0
    // [0,2)        0       0      0       0
    // [0,1)        0       0      0       0
    // (-inf,0)     1       0      0       0
    std::vector<std::vector<double>> m_expected({{1.,0.,0.,0.},
                                                 {0.,0.,0.,0.},
                                                 {0.,0.,0.,0.},
                                                 {0.,0.,0.,0.}}) ;
    // check histograms
    auto& index_pairs = model1.getPairIndex() ;
    auto& histograms1_ipd = model1.getHistogramsIPD() ;
    auto& histograms1_pwd = model1.getHistogramsPWD() ;
    for(const auto& index_pair : index_pairs)
    {   const size_t& pos1 = index_pair.first ;
        const size_t& pos2 = index_pair.second ;

        std::vector<std::vector<double>> m_ipd = 
                ngsai::hist2d_to_matrix(histograms1_ipd[pos1][pos2]) ;
        EXPECT_THAT(m_ipd, m_expected) ;

        std::vector<std::vector<double>> m_pwd = 
                ngsai::hist2d_to_matrix(histograms1_pwd[pos1][pos2]) ;
        EXPECT_THAT(m_pwd, m_expected) ;
    }

    // assign move
    model3 = model1 ;
    DiPositionNormalizedKineticModelToy model4(std::move(model3)) ;
    ASSERT_EQ(model1, model4) ;
}


/*!
 * \brief Tests that the copy() method works properly.
 */
TEST_F(DiPositionNormalizedKineticModelTest, copy)
{
    // copy empty model
    DiPositionNormalizedKineticModelToy model1(*kmermap) ;
    DiPositionNormalizedKineticModelToy* copy(nullptr);
    copy = model1.copy() ;
    ASSERT_EQ(model1, *copy) ;
    delete copy ;
    copy = nullptr ;

    // 2+2 bins on both axes (-inf,0), [0,1), [0,2), [2,inf)
    model1.setParameters(3, 0., 2., 2, 0.) ;
    model1.add(kinetics[0]) ; // add +1 to bin (-inf,0),(-inf,0) at each pos
    // at this moment the model contains, for IPD and PWD, for all 
    // POSn vs POSn+1 histogram
    //          (-inf,0), [0,1), [0,2), [2,inf)
    // [2,inf)      0       0      0       0
    // [0,2)        0       0      0       0
    // [0,1)        0       0      0       0
    // (-inf,0)     1       0      0       0
    std::vector<std::vector<double>> m_expected({{1.,0.,0.,0.},
                                                 {0.,0.,0.,0.},
                                                 {0.,0.,0.,0.},
                                                 {0.,0.,0.,0.}}) ;
    // check histograms
    auto& index_pairs = model1.getPairIndex() ;
    auto& histograms1_ipd = model1.getHistogramsIPD() ;
    auto& histograms1_pwd = model1.getHistogramsPWD() ;
    for(const auto& index_pair : index_pairs)
    {   const size_t& pos1 = index_pair.first ;
        const size_t& pos2 = index_pair.second ;

        std::vector<std::vector<double>> m_ipd = 
                ngsai::hist2d_to_matrix(histograms1_ipd[pos1][pos2]) ;
        EXPECT_THAT(m_ipd, m_expected) ;

        std::vector<std::vector<double>> m_pwd = 
                ngsai::hist2d_to_matrix(histograms1_pwd[pos1][pos2]) ;
        EXPECT_THAT(m_pwd, m_expected) ;
    }

    // copy
    copy = model1.copy() ;
    ASSERT_EQ(model1, *copy) ;
    
    delete copy ;
    copy = nullptr ;
}


/*!
 * \brief Tests that the add() method with KineticModel
 * works properly.
 */
TEST_F(DiPositionNormalizedKineticModelTest, plusequal_operator)
{
    DiPositionNormalizedKineticModelToy tmp ;
    DiPositionNormalizedKineticModelToy model_empty_uninit ;
    DiPositionNormalizedKineticModelToy model_empty ;
    model_empty.setParameters(3, 0., 2., 2, 0.) ;
    DiPositionNormalizedKineticModelToy model_uninit(*kmermap) ;
    DiPositionNormalizedKineticModelToy model(*kmermap) ;
    model.setParameters(3, 0., 2., 2, 0.) ;


    // must fail, both operandes are not initialised
    EXPECT_THROW(model_empty_uninit.add(model_empty_uninit), std::invalid_argument) ;
    EXPECT_THROW(model_empty_uninit.add(model_uninit),       std::invalid_argument) ;
    EXPECT_THROW(model_uninit.add(model_empty_uninit),       std::invalid_argument) ;
    EXPECT_THROW(model_uninit.add(model_uninit),             std::invalid_argument) ;

    // must fail, one operandes is not initialised
    EXPECT_THROW(model_empty_uninit.add(model_empty), std::invalid_argument) ;
    EXPECT_THROW(model_empty.add(model_empty_uninit), std::invalid_argument) ;
    EXPECT_THROW(model_empty.add(model_uninit),       std::invalid_argument) ; 
    EXPECT_THROW(model_uninit.add(model_empty),       std::invalid_argument) ; 

    DiPositionNormalizedKineticModelToy model1(*kmermap) ;
    DiPositionNormalizedKineticModelToy model2(*kmermap) ;
    DiPositionNormalizedKineticModelToy model3 ;


    // 2+2 bins on both axes (-inf,0), [0,1), [0,2), [2,inf)
    model1.setParameters(3, 0., 2., 2, 0.) ;
    // 2+2 bins (-inf,0), [0,1), [0,2), [2,inf)
    model2.setParameters(6, 0., 2., 2, 0.) ;

    // must fail, not same size
    EXPECT_THROW(model1.add(model2), std::invalid_argument) ;

    // 1+2 bins (-inf,-1), [-1,1), [1,inf)
    model2.setParameters(3, -1, 1, 1, 0.) ;

    // must fail, not same nb of bins
    EXPECT_THROW(model1.add(model2), std::runtime_error) ;

    // train model
    for(const auto& kinetic : kinetics)
    {   model1.add(kinetic) ; }
    // at this moment the model contains, for IPD and PWD, for all 
    // POSn vs POSn+1 histogram
    //          (-inf,0), [0,1), [0,2), [2,inf)
    // [2,inf)      5       5      5       5
    // [0,2)        5       1      1       1
    // [0,1)        5       1      1       1
    // (-inf,0)     9       1      1       1
    std::vector<std::vector<double>> m_expected({{9.,5.,5.,5.},
                                                 {5.,1.,1.,1.},
                                                 {5.,1.,1.,1.},
                                                 {5.,1.,1.,1.}}) ;
    std::vector<std::pair<size_t,size_t>> pair_exp = {std::make_pair(0,1),
                                                      std::make_pair(0,2),
                                                      std::make_pair(1,2)} ;
    
    // reset models
    // 2+2 bins on both axes (-inf,0), [0,1), [0,2), [2,inf)
    model2.setParameters(3, 0., 2., 2, 0.) ;
    model3.setParameters(3, 0., 2., 2, 0.) ;


    // sum empty models
    tmp = model_empty ;
    model_empty.add(model_empty) ;
    ASSERT_EQ(model_empty, tmp) ;


    // sum trained model and empty model
    // 2+2 bins on both axes (-inf,0), [0,1), [0,2), [2,inf)
    model2.setParameters(3, 0., 2., 2, 0.) ;
    tmp = model1 ;
    model1.add(model2) ;
    ASSERT_EQ(tmp, model1) ;
    // model2 is 2x model1
    tmp = model2 ;
    tmp.add(model1) ;
    tmp.add(model1) ;
    model2.add(model1) ;
    model2.add(model1) ;
    ASSERT_EQ(model2, tmp) ;
    // model3 is 3x model1
    tmp = model3 ;
    tmp.add(model1) ;
    tmp.add(model1) ;
    tmp.add(model1) ;
    model3.add(model1) ;
    model3.add(model1) ;
    model3.add(model1) ;
    ASSERT_EQ(model3, tmp) ;
    

    // reset models
    // 2+2 bins on both axes (-inf,0), [0,1), [0,2), [2,inf)
    model2.setParameters(3, 0., 2., 2, 0.) ;
    model3.setParameters(3, 0., 2., 2, 0.) ;


    // sum density models
    model1.density() ;
    // model2 is 2x model1
    tmp = model2 ;
    tmp.add(model1) ;
    tmp.add(model1) ;
    model2.add(model1) ;
    model2.add(model1) ;
    ASSERT_EQ(model2, tmp) ;
    // model3 is 3x model1
    tmp = model3 ;
    tmp.add(model1) ;
    tmp.add(model1) ;
    tmp.add(model1) ;
    model3.add(model1) ;
    model3.add(model1) ;
    model3.add(model1) ;
    ASSERT_EQ(model3, tmp) ;
    

    // reset models
    // 2+2 bins on both axes (-inf,0), [0,1), [0,2), [2,inf)
    model2.setParameters(3, 0., 2., 2, 0.) ;
    model3.setParameters(3, 0., 2., 2, 0.) ;


    // sum log density models
    model1.log() ;
   
    // model2 is 2x model1
    tmp = model2 ;
    tmp.add(model1) ;
    tmp.add(model1) ;
    model2.add(model1) ;
    model2.add(model1) ;
    ASSERT_EQ(model2, tmp) ;
    // model3 is 3x model1
    tmp = model3 ;
    tmp.add(model1) ;
    tmp.add(model1) ;
    tmp.add(model1) ;
    model3.add(model1) ;
    model3.add(model1) ;
    model3.add(model1) ;
    ASSERT_EQ(model3, tmp) ;
}


/*!
 * \brief Tests that the logLikelihood() method works properly.
 */
TEST_F(DiPositionNormalizedKineticModelTest, logLikelihood)
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

    DiPositionNormalizedKineticModelToy model ;

    // no kmermap
    EXPECT_THROW(model.logLikelihood(s1_a), std::runtime_error) ;

    model = DiPositionNormalizedKineticModelToy(*kmermap) ;

    // not initialized
    EXPECT_THROW(model.logLikelihood(s1_a), std::runtime_error) ;

    // 2+2 bins on both axes (-inf,0), [0,1), [0,2), [2,inf)
    model.setParameters(3, 0., 2., 2, 0.) ;

    // must fail not log density
    EXPECT_THROW(model.logLikelihood(kinetics[0]), std::runtime_error) ;

    // train model
    for(const auto& kinetic : kinetics)
    {   model.add(kinetic) ; }
    // at this moment the model contains, for IPD and PWD, for all 
    // POSn vs POSn+1 histogram
    //          (-inf,0), [0,1), [1,2), [2,inf)
    // [2,inf)      5       1      1       1
    // [1,2)        5       1      1       1
    // [0,1)        5       1      1       1
    // (-inf,0)     9       5      5       5

    model.density() ;
    model.log() ;

    double n = 48 ;
    ASSERT_DOUBLE_EQ(model.logLikelihood(kinetics[0]),
                     2.*(log(9/n) + log(9/n))) ;
    ASSERT_DOUBLE_EQ(model.logLikelihood(kinetics[1]),
                     2.*(log(5/n) + log(9/n))) ;
    ASSERT_DOUBLE_EQ(model.logLikelihood(kinetics[2]),
                     2.*(log(5/n) + log(9/n))) ;
    ASSERT_DOUBLE_EQ(model.logLikelihood(kinetics[3]),
                     2.*(log(5/n) + log(9/n))) ;
    ASSERT_DOUBLE_EQ(model.logLikelihood(kinetics[4]),
                     2.*(log(5/n) + log(5/n))) ;
    ASSERT_DOUBLE_EQ(model.logLikelihood(kinetics[5]),
                     2.*(log(1/n) + log(5/n))) ;
    ASSERT_DOUBLE_EQ(model.logLikelihood(kinetics[6]),
                     2.*(log(1/n) + log(5/n))) ;
    ASSERT_DOUBLE_EQ(model.logLikelihood(kinetics[7]),
                     2.*(log(1/n) + log(5/n))) ;
    ASSERT_DOUBLE_EQ(model.logLikelihood(kinetics[8]),
                     2.*(log(5/n) + log(5/n))) ;
    ASSERT_DOUBLE_EQ(model.logLikelihood(kinetics[9]),
                     2.*(log(1/n) + log(5/n))) ;
    ASSERT_DOUBLE_EQ(model.logLikelihood(kinetics[10]),
                     2.*(log(1/n) + log(5/n))) ;
    ASSERT_DOUBLE_EQ(model.logLikelihood(kinetics[11]),
                     2.*(log(1/n) + log(5/n))) ;
    ASSERT_DOUBLE_EQ(model.logLikelihood(kinetics[12]),
                     2.*(log(5/n) + log(5/n))) ;
    ASSERT_DOUBLE_EQ(model.logLikelihood(kinetics[13]),
                     2.*(log(1/n) + log(5/n))) ;
    ASSERT_DOUBLE_EQ(model.logLikelihood(kinetics[14]),
                     2.*(log(1/n) + log(5/n))) ;
    ASSERT_DOUBLE_EQ(model.logLikelihood(kinetics[15]),
                     2.*(log(1/n) + log(5/n))) ;

    ASSERT_DOUBLE_EQ(model.logLikelihood(kinetics[16]),
                     2.*(log(9/n) + log(9/n))) ;
    ASSERT_DOUBLE_EQ(model.logLikelihood(kinetics[17]),
                     2.*(log(5/n) + log(9/n))) ;
    ASSERT_DOUBLE_EQ(model.logLikelihood(kinetics[18]),
                     2.*(log(5/n) + log(9/n))) ;
    ASSERT_DOUBLE_EQ(model.logLikelihood(kinetics[19]),
                     2.*(log(5/n) + log(9/n))) ;
    ASSERT_DOUBLE_EQ(model.logLikelihood(kinetics[20]),
                     2.*(log(9/n) + log(5/n))) ;
    ASSERT_DOUBLE_EQ(model.logLikelihood(kinetics[21]),
                     2.*(log(5/n) + log(5/n))) ;
    ASSERT_DOUBLE_EQ(model.logLikelihood(kinetics[22]),
                     2.*(log(5/n) + log(5/n))) ;
    ASSERT_DOUBLE_EQ(model.logLikelihood(kinetics[23]),
                     2.*(log(5/n) + log(5/n))) ;
    ASSERT_DOUBLE_EQ(model.logLikelihood(kinetics[24]),
                     2.*(log(9/n) + log(5/n))) ;
    ASSERT_DOUBLE_EQ(model.logLikelihood(kinetics[25]),
                     2.*(log(5/n) + log(5/n))) ;
    ASSERT_DOUBLE_EQ(model.logLikelihood(kinetics[26]),
                     2.*(log(5/n) + log(5/n))) ;
    ASSERT_DOUBLE_EQ(model.logLikelihood(kinetics[27]),
                     2.*(log(5/n) + log(5/n))) ;
    ASSERT_DOUBLE_EQ(model.logLikelihood(kinetics[28]),
                     2.*(log(9/n) + log(5/n))) ;
    ASSERT_DOUBLE_EQ(model.logLikelihood(kinetics[29]),
                     2.*(log(5/n) + log(5/n))) ;
    ASSERT_DOUBLE_EQ(model.logLikelihood(kinetics[30]),
                     2.*(log(5/n) + log(5/n))) ;
    ASSERT_DOUBLE_EQ(model.logLikelihood(kinetics[31]),
                     2.*(log(5/n) + log(5/n))) ;

    ASSERT_DOUBLE_EQ(model.logLikelihood(kinetics[32]),
                     2.*(log(9/n) + log(9/n))) ;
    ASSERT_DOUBLE_EQ(model.logLikelihood(kinetics[33]),
                     2.*(log(5/n) + log(5/n))) ;
    ASSERT_DOUBLE_EQ(model.logLikelihood(kinetics[34]),
                     2.*(log(5/n) + log(5/n))) ;
    ASSERT_DOUBLE_EQ(model.logLikelihood(kinetics[35]),
                     2.*(log(5/n) + log(5/n))) ;
    ASSERT_DOUBLE_EQ(model.logLikelihood(kinetics[36]),
                     2.*(log(9/n) + log(5/n))) ;
    ASSERT_DOUBLE_EQ(model.logLikelihood(kinetics[37]),
                     2.*(log(5/n) + log(1/n))) ;
    ASSERT_DOUBLE_EQ(model.logLikelihood(kinetics[38]),
                     2.*(log(5/n) + log(1/n))) ;
    ASSERT_DOUBLE_EQ(model.logLikelihood(kinetics[39]),
                     2.*(log(5/n) + log(1/n))) ;
    ASSERT_DOUBLE_EQ(model.logLikelihood(kinetics[40]),
                     2.*(log(9/n) + log(5/n))) ;
    ASSERT_DOUBLE_EQ(model.logLikelihood(kinetics[41]),
                     2.*(log(5/n) + log(1/n))) ;
    ASSERT_DOUBLE_EQ(model.logLikelihood(kinetics[42]),
                     2.*(log(5/n) + log(1/n))) ;
    ASSERT_DOUBLE_EQ(model.logLikelihood(kinetics[43]),
                     2.*(log(5/n) + log(1/n))) ;
    ASSERT_DOUBLE_EQ(model.logLikelihood(kinetics[44]),
                     2.*(log(9/n) + log(5/n))) ;
    ASSERT_DOUBLE_EQ(model.logLikelihood(kinetics[45]),
                     2.*(log(5/n) + log(1/n))) ;
    ASSERT_DOUBLE_EQ(model.logLikelihood(kinetics[46]),
                     2.*(log(5/n) + log(1/n))) ;
    ASSERT_DOUBLE_EQ(model.logLikelihood(kinetics[47]),
                     2.*(log(5/n) + log(1/n))) ;

    // must fail, size of kinetic signal does not match size of model
    std::string seq_e1("A") ;
    std::string seq_e2("AA") ;
    std::string seq_e4("AAAA") ;
    std::vector<double> val_e1({1}) ;
    std::vector<double> val_e2({1,1}) ;
    std::vector<double> val_e4({1,1,1,1}) ;
    ngsai::KineticSignal s1_err(seq_e1,seq_e1,val_e1,val_e1,val_e1,val_e1) ;
    ngsai::KineticSignal s2_err(seq_e2,seq_e2,val_e2,val_e2,val_e2,val_e2) ;
    ngsai::KineticSignal s4_err(seq_e4,seq_e4,val_e4,val_e4,val_e4,val_e4) ;
    EXPECT_THROW(model.logLikelihood(s1_err), std::invalid_argument) ; // tpo short
    EXPECT_THROW(model.logLikelihood(s2_err), std::invalid_argument) ; // too short
    EXPECT_THROW(model.logLikelihood(s4_err), std::invalid_argument) ; // too long
}