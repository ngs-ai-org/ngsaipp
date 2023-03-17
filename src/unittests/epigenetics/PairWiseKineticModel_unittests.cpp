#include <gtest/gtest.h>
#include <gmock/gmock.h>
#include <string>
#include <vector>
#include <set>
#include <stdexcept>
#include <limits>
#include <cmath>

#include <ngsaipp/epigenetics/PairWiseKineticModel.hpp>
#include "utility.hpp"


using testing::ElementsAre ;
using testing::ContainerEq ;


/*
 * This file ontains tests for the ngsai::PairWiseKineticModel class from 
 * src/epigenetics/PairWiseKineticModel.cpp
 */



/*!
 * \brief A child class of RawKineticModel allowing to access the protected 
 * fields of PairWiseKineticModel for testing purposes.
 */
class PairWiseKineticModelToy : public ngsai::PairWiseKineticModel
{
    public:

        PairWiseKineticModelToy()
            : PairWiseKineticModel()
        { ; }

        PairWiseKineticModelToy(const ngsai::PairWiseKineticModel& other)
            : PairWiseKineticModel(other)
        { ; }

        PairWiseKineticModelToy(ngsai::PairWiseKineticModel&& other)
            : PairWiseKineticModel(other)
        { ; }

        ~PairWiseKineticModelToy()
        { ; }

        // ugly workaround to test RawKineticModel.copy() using a
        // RawKineticModelTest instance
        PairWiseKineticModelToy* copy() const
        {   PairWiseKineticModel* copy = PairWiseKineticModel::copy() ;
            PairWiseKineticModelToy* copy2 = new PairWiseKineticModelToy(*copy) ;
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

        bool& getIsInit()
        {   return m_is_init ; }

        bool& getIsDensity()
        {   return m_is_density ; }

        bool& getIsLog()
        {   return m_is_log ; }

} ;


class PairWiseKineticModelTest : public ::testing::Test
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

        static void
        SetUpTestSuite()
        {
            if(kinetics.size() == 0)
            {   
                std::vector<std::vector<double>> data ;
                data.push_back(std::vector<double>({-1.,-1.,-1.})) ;
                data.push_back(std::vector<double>({ 0.,-1.,-1.})) ;
                data.push_back(std::vector<double>({ 1.,-1.,-1.})) ;
                data.push_back(std::vector<double>({ 2.,-1.,-1.})) ;
                data.push_back(std::vector<double>({-1., 0.,-1.})) ;
                data.push_back(std::vector<double>({ 0., 0.,-1.})) ;
                data.push_back(std::vector<double>({ 1., 0.,-1.})) ;
                data.push_back(std::vector<double>({ 2., 0.,-1.})) ;
                data.push_back(std::vector<double>({-1., 1.,-1.})) ;
                data.push_back(std::vector<double>({ 0., 1.,-1.})) ;
                data.push_back(std::vector<double>({ 1., 1.,-1.})) ;
                data.push_back(std::vector<double>({ 2., 1.,-1.})) ;
                data.push_back(std::vector<double>({-1., 2.,-1.})) ;
                data.push_back(std::vector<double>({ 0., 2.,-1.})) ;
                data.push_back(std::vector<double>({ 1., 2.,-1.})) ;
                data.push_back(std::vector<double>({ 2., 2.,-1.})) ;
                
                data.push_back(std::vector<double>({-1, -1.,-1.})) ;
                data.push_back(std::vector<double>({ 0, -1.,-1.})) ;
                data.push_back(std::vector<double>({ 1, -1.,-1.})) ;
                data.push_back(std::vector<double>({ 2, -1.,-1.})) ;
                data.push_back(std::vector<double>({-1, -1., 0.})) ;
                data.push_back(std::vector<double>({ 0, -1., 0.})) ;
                data.push_back(std::vector<double>({ 1, -1., 0.})) ;
                data.push_back(std::vector<double>({ 2, -1., 0.})) ;
                data.push_back(std::vector<double>({-1, -1., 1.})) ;
                data.push_back(std::vector<double>({ 0, -1., 1.})) ;
                data.push_back(std::vector<double>({ 1, -1., 1.})) ;
                data.push_back(std::vector<double>({ 2, -1., 1.})) ;
                data.push_back(std::vector<double>({-1, -1., 2.})) ;
                data.push_back(std::vector<double>({ 0, -1., 2.})) ;
                data.push_back(std::vector<double>({ 1, -1., 2.})) ;
                data.push_back(std::vector<double>({ 2, -1., 2.})) ;
                
                data.push_back(std::vector<double>({-1, -1.,-1.})) ;
                data.push_back(std::vector<double>({-1,  0.,-1.})) ;
                data.push_back(std::vector<double>({-1,  1.,-1.})) ;
                data.push_back(std::vector<double>({-1,  2.,-1.})) ;
                data.push_back(std::vector<double>({-1, -1., 0.})) ;
                data.push_back(std::vector<double>({-1,  0., 0.})) ;
                data.push_back(std::vector<double>({-1,  1., 0.})) ;
                data.push_back(std::vector<double>({-1,  2., 0.})) ;
                data.push_back(std::vector<double>({-1, -1., 1.})) ;
                data.push_back(std::vector<double>({-1,  0., 1.})) ;
                data.push_back(std::vector<double>({-1,  1., 1.})) ;
                data.push_back(std::vector<double>({-1,  2., 1.})) ;
                data.push_back(std::vector<double>({-1, -1., 2.})) ;
                data.push_back(std::vector<double>({-1,  0., 2.})) ;
                data.push_back(std::vector<double>({-1,  1., 2.})) ;
                data.push_back(std::vector<double>({-1,  2., 2.})) ;
                
                sequence = std::string("AAA") ;
                for(const auto& datum : data)
                {   kinetics.push_back(ngsai::KineticSignal(sequence, sequence,
                                                            datum,    datum,
                                                            datum,    datum)) ;
                }
            }
        }
} ;


// initialise static members
std::string                       PairWiseKineticModelTest::sequence = std::string() ;
std::vector<ngsai::KineticSignal> PairWiseKineticModelTest::kinetics = std::vector<ngsai::KineticSignal>() ;


/*!
 * \brief Tests that constructor works properly.
 */
TEST_F(PairWiseKineticModelTest, constructor)
{
    PairWiseKineticModelToy model ;

    ASSERT_EQ(model.getPairNb(), 0) ;
    ASSERT_EQ(model.getPairIndex().size(), 0) ;
    ASSERT_EQ(model.getHistogramsIPD().size(), 0) ;
    ASSERT_EQ(model.getHistogramsPWD().size(), 0) ;
    ASSERT_EQ(model.isInit(), false) ;
    ASSERT_EQ(model.isDensity(), false) ;
    ASSERT_EQ(model.isLog(), false) ;
}


/*!
 * \brief Tests that isInit() method works properly.
 */
TEST_F(PairWiseKineticModelTest, isInit)
{
    PairWiseKineticModelToy model ;

    ASSERT_EQ(model.isInit(), false) ;
    bool& init = model.getIsInit() ;
    init = true ;
    ASSERT_EQ(model.isInit(), true) ;
}


/*!
 * \brief Tests that isDensity() method works properly.
 */
TEST_F(PairWiseKineticModelTest, isDensity)
{
    PairWiseKineticModelToy model ;

    ASSERT_EQ(model.isDensity(), false) ;
    bool& density = model.getIsDensity() ;
    density = true ;
    ASSERT_EQ(model.isDensity(), true) ;
}


/*!
 * \brief Tests that isLog() method works properly.
 */
TEST_F(PairWiseKineticModelTest, isLog)
{
    PairWiseKineticModelToy model ;

    ASSERT_EQ(model.isLog(), false) ;
    bool& islog = model.getIsLog() ;
    islog = true ;
    ASSERT_EQ(model.isLog(), true) ;
}


/*!
 * \brief Tests that the setParameters() method works properly and that 
 * everything is initialised properly or throws the expected exceptions.
 */
TEST_F(PairWiseKineticModelTest, setParameters)
{   
    size_t size(5) ;
    double min(-2.) ;
    double max(2.) ;
    size_t n_bin(4) ;
    double pseudocounts(10.) ;

    // it must work properly
    PairWiseKineticModelToy model ;
    model.setParameters(size,
                        min,
                        max,
                        n_bin,
                        pseudocounts) ;
    // expected pairs of positions to consider
    std::vector<std::pair<size_t,size_t>> pair_index_exp = 
            {std::pair<size_t,size_t>(0,1), std::pair<size_t,size_t>(0,2),
             std::pair<size_t,size_t>(0,3), std::pair<size_t,size_t>(0,4),
             std::pair<size_t,size_t>(1,2), std::pair<size_t,size_t>(1,3),
             std::pair<size_t,size_t>(1,4), std::pair<size_t,size_t>(2,3),
             std::pair<size_t,size_t>(2,4), std::pair<size_t,size_t>(3,4)} ;

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
 * \brief Tests that the size() method works properly.
 */
TEST_F(PairWiseKineticModelTest, size)
{   // before initialisation
    PairWiseKineticModelToy model ;
    ASSERT_EQ(model.size(), 0) ;
    
    // 2bp
    model.setParameters(2, -2, 2, 4, 10) ;
    ASSERT_EQ(model.size(), 2) ;

    // 5bp
    model.setParameters(5, -2, 2, 4, 10) ;
    ASSERT_EQ(model.size(), 5) ;

    // 10bp
    model.setParameters(10, -2, 2, 4, 10) ;
    ASSERT_EQ(model.size(), 10) ;
}


/*!
 * \brief Tests that the getBinBoundaries() method works properly.
 */
TEST_F(PairWiseKineticModelTest, getBinBoundaries)
{   
    PairWiseKineticModelToy model ;

    // 4+2 bins [-inf,-2), [-2,-1), [-1,0), [0,1), [1,2), [2,inf)
    model.setParameters(5, -2, 2, 4, 10) ;
    double inf = std::numeric_limits<double>::infinity() ;
    std::vector<std::pair<double,double>> bins_exp = 
                                            {std::make_pair(-inf, -2.),
                                             std::make_pair(-2., -1.),
                                             std::make_pair(-1.,  0.),
                                             std::make_pair( 0.,  1.),
                                             std::make_pair( 1.,  2.),
                                             std::make_pair( 2.,  inf)} ;
    std::vector<std::pair<double,double>> bins = model.getBinBoundaries() ;
    EXPECT_THAT(bins, ContainerEq(bins_exp)) ;

    // 1+2 bins [-inf,-2), [-10,10), [2,inf)
    model.setParameters(5, -10, 10, 1, 10) ;
    bins_exp = {std::make_pair(-inf, -10.),
                std::make_pair(-10., 10.),
                std::make_pair( 10.,  inf)} ;
    bins = model.getBinBoundaries() ;
    EXPECT_THAT(bins, ContainerEq(bins_exp)) ;
}


/*!
 * \brief Tests that the add() method with KineticSignal works properly.
 */
TEST_F(PairWiseKineticModelTest, add_kineticsignal)
{   
    PairWiseKineticModelToy model ;

    // not yet initialised
    EXPECT_THROW(model.add(kinetics[0]), std::runtime_error) ;

    // 2+2 bins on both axes (-inf,0), [0,1), [0,2), [2,inf)
    model.setParameters(3, 0., 2., 2, 0.) ;

    // train model
    for(const auto& kinetic : kinetics)
    {   model.add(kinetic) ; }
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
    
    // check histograms
    auto& index_pairs = model.getPairIndex() ;

    auto& histograms_ipd = model.getHistogramsIPD() ;
    auto& histograms_pwd = model.getHistogramsPWD() ;
    for(const auto& index_pair : index_pairs)
    {   const size_t& pos1 = index_pair.first ;
        const size_t& pos2 = index_pair.second ;

        std::vector<std::vector<double>> m_ipd = 
                ngsai::hist2d_to_matrix(histograms_ipd[pos1][pos2]) ;
        EXPECT_THAT(m_ipd, m_expected) ;

        std::vector<std::vector<double>> m_pwd = 
                ngsai::hist2d_to_matrix(histograms_pwd[pos1][pos2]) ;
        EXPECT_THAT(m_pwd, m_expected) ;
    }

    std::vector<double> val_e1({1.}) ;
    std::vector<double> val_e2({1.,2.}) ;
    std::vector<double> val_e4({1.,2.,3.,4.}) ;
    std::string seq_e1("A") ;
    std::string seq_e2("AA") ;
    std::string seq_e4("AAAA") ;
    ngsai::KineticSignal k_err0 ;
    ngsai::KineticSignal k_err1(seq_e1, seq_e1, val_e1, val_e1, val_e1, val_e1) ;
    ngsai::KineticSignal k_err2(seq_e2, seq_e2, val_e2, val_e2, val_e2, val_e2) ;
    ngsai::KineticSignal k_err4(seq_e4, seq_e4, val_e4, val_e4, val_e4, val_e4) ;

    // must fail because not same size as model
    EXPECT_THROW(model.add(k_err0), std::invalid_argument) ; // empty
    EXPECT_THROW(model.add(k_err1), std::invalid_argument) ; // too short
    EXPECT_THROW(model.add(k_err2), std::invalid_argument) ; // too short
    EXPECT_THROW(model.add(k_err4), std::invalid_argument) ; // too long

}


/*!
 * \brief Tests that the copy constructor works properly.
 */
TEST_F(PairWiseKineticModelTest, constructor_copy)
{
    // copy empty model
    PairWiseKineticModelToy model1 ;
    PairWiseKineticModelToy model2(model1) ;
    ASSERT_EQ(model1.isInit(),    model2.isInit()) ;
    ASSERT_EQ(model1.isDensity(), model2.isDensity()) ;
    ASSERT_EQ(model1.isLog(),     model2.isLog()) ;
    EXPECT_THAT(model1.getPairIndex(), 
                ContainerEq(model2.getPairIndex())) ;
    EXPECT_THAT(model1.getHistogramsIPD(), 
                ContainerEq(model2.getHistogramsIPD())) ;
    EXPECT_THAT(model1.getHistogramsPWD(), 
                ContainerEq(model2.getHistogramsPWD())) ;

    
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
    PairWiseKineticModelToy model3(model1) ;
    auto& histograms3_ipd = model3.getHistogramsIPD() ;
    auto& histograms3_pwd = model3.getHistogramsPWD() ;

    // compare them
    ASSERT_EQ(model1.isInit(),    model3.isInit()) ;
    ASSERT_EQ(model1.isDensity(), model3.isDensity()) ;
    ASSERT_EQ(model1.isLog(),     model3.isLog()) ;
    EXPECT_THAT(model1.getPairIndex(), ContainerEq(model3.getPairIndex())) ;
    EXPECT_THAT(histograms1_ipd, ContainerEq(histograms3_ipd)) ;
    EXPECT_THAT(histograms1_pwd, ContainerEq(histograms3_pwd)) ;

}


/*!
 * \brief Tests that the move constructor works properly.
 */
TEST_F(PairWiseKineticModelTest, constructor_move)
{
    // assign empty model
    PairWiseKineticModelToy model1 ;
    PairWiseKineticModelToy model2(model1) ;
    PairWiseKineticModelToy model3(std::move(model2)) ;
    ASSERT_EQ(model1.isInit(),    model3.isInit()) ;
    ASSERT_EQ(model1.isDensity(), model3.isDensity()) ;
    ASSERT_EQ(model1.isLog(),     model3.isLog()) ;
    EXPECT_THAT(model1.getPairIndex(), 
                ContainerEq(model3.getPairIndex())) ;
    EXPECT_THAT(model1.getHistogramsIPD(), 
                ContainerEq(model3.getHistogramsIPD())) ;
    EXPECT_THAT(model1.getHistogramsPWD(), 
                ContainerEq(model3.getHistogramsPWD())) ;


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
    PairWiseKineticModelToy model4(model1) ;
    PairWiseKineticModelToy model5(std::move(model4)) ;
    auto& histograms5_ipd = model5.getHistogramsIPD() ;
    auto& histograms5_pwd = model5.getHistogramsPWD() ;

    // compare
    ASSERT_EQ(model1.isInit(),    model5.isInit()) ;
    ASSERT_EQ(model1.isDensity(), model5.isDensity()) ;
    ASSERT_EQ(model1.isLog(),     model5.isLog()) ;
    EXPECT_THAT(model1.getPairIndex(), ContainerEq(model5.getPairIndex())) ;
    EXPECT_THAT(histograms1_ipd, ContainerEq(histograms5_ipd)) ;
    EXPECT_THAT(histograms1_pwd, ContainerEq(histograms5_pwd)) ;
}


/*!
 * \brief Tests that the assignment operator works properly.
 */
TEST_F(PairWiseKineticModelTest, assignment_operator)
{
    // assign empty model
    PairWiseKineticModelToy model1 ;
    PairWiseKineticModelToy model2 = model1 ;
    ASSERT_EQ(model1.isInit(),    model2.isInit()) ;
    ASSERT_EQ(model1.isDensity(), model2.isDensity()) ;
    ASSERT_EQ(model1.isLog(),     model2.isLog()) ;
    EXPECT_THAT(model1.getPairIndex(), 
                ContainerEq(model2.getPairIndex())) ;
    EXPECT_THAT(model1.getHistogramsIPD(), 
                ContainerEq(model2.getHistogramsIPD())) ;
    EXPECT_THAT(model1.getHistogramsPWD(), 
                ContainerEq(model2.getHistogramsPWD())) ;

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
    auto& histograms2_ipd = model2.getHistogramsIPD() ;
    auto& histograms2_pwd = model2.getHistogramsPWD() ;

    // compare
    ASSERT_EQ(model1.isInit(),    model2.isInit()) ;
    ASSERT_EQ(model1.isDensity(), model2.isDensity()) ;
    ASSERT_EQ(model1.isLog(),     model2.isLog()) ;
    EXPECT_THAT(model1.getPairIndex(), ContainerEq(model2.getPairIndex())) ;
    EXPECT_THAT(histograms1_ipd, ContainerEq(histograms2_ipd)) ;
    EXPECT_THAT(histograms1_pwd, ContainerEq(histograms2_pwd)) ;
}


/*!
 * \brief Tests that the move assignment operator works properly.
 */
TEST_F(PairWiseKineticModelTest, assignment_move_operator)
{
    // assign empty model
    PairWiseKineticModelToy model1 ;
    PairWiseKineticModelToy model2 = model1 ;
    PairWiseKineticModelToy model3(std::move(model2)) ;
    ASSERT_EQ(model1.isInit(),    model3.isInit()) ;
    ASSERT_EQ(model1.isDensity(), model3.isDensity()) ;
    ASSERT_EQ(model1.isLog(),     model3.isLog()) ;
    EXPECT_THAT(model1.getPairIndex(), 
                ContainerEq(model3.getPairIndex())) ;
    EXPECT_THAT(model1.getHistogramsIPD(), 
                ContainerEq(model3.getHistogramsIPD())) ;
    EXPECT_THAT(model1.getHistogramsPWD(), 
                ContainerEq(model3.getHistogramsPWD())) ;

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
    PairWiseKineticModelToy model4(std::move(model3)) ;
    auto& histograms4_ipd = model4.getHistogramsIPD() ;
    auto& histograms4_pwd = model4.getHistogramsPWD() ;

    // compare
    ASSERT_EQ(model1.isInit(),    model4.isInit()) ;
    ASSERT_EQ(model1.isDensity(), model4.isDensity()) ;
    ASSERT_EQ(model1.isLog(),     model4.isLog()) ;
    EXPECT_THAT(model1.getPairIndex(), ContainerEq(model4.getPairIndex())) ;
    EXPECT_THAT(histograms1_ipd, ContainerEq(histograms4_ipd)) ;
    EXPECT_THAT(histograms1_pwd, ContainerEq(histograms4_pwd)) ;
}


/*!
 * \brief Tests that the copy() method works properly.
 */
TEST_F(PairWiseKineticModelTest, copy)
{
    // copy empty model
    PairWiseKineticModelToy model1 ;
    PairWiseKineticModelToy* copy(nullptr);
    copy = model1.copy() ;
    ASSERT_EQ(model1.isInit(),    copy->isInit()) ;
    ASSERT_EQ(model1.isDensity(), copy->isDensity()) ;
    ASSERT_EQ(model1.isLog(),     copy->isLog()) ;
    EXPECT_THAT(model1.getPairIndex(), 
                ContainerEq(copy->getPairIndex())) ;
    EXPECT_THAT(model1.getHistogramsIPD(), 
                ContainerEq(copy->getHistogramsIPD())) ;
    EXPECT_THAT(model1.getHistogramsPWD(), 
                ContainerEq(copy->getHistogramsPWD())) ;
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
    auto& histogramscp_ipd = copy->getHistogramsIPD() ;
    auto& histogramscp_pwd = copy->getHistogramsPWD() ;

    // compare
    ASSERT_EQ(model1.isInit(),    copy->isInit()) ;
    ASSERT_EQ(model1.isDensity(), copy->isDensity()) ;
    ASSERT_EQ(model1.isLog(),     copy->isLog()) ;
    EXPECT_THAT(model1.getPairIndex(), ContainerEq(copy->getPairIndex())) ;
    EXPECT_THAT(histograms1_ipd, ContainerEq(histogramscp_ipd)) ;
    EXPECT_THAT(histograms1_pwd, ContainerEq(histogramscp_pwd)) ;
    
    delete copy ;
    copy = nullptr ;
}


/*!
 * \brief Tests that the density() method works properly. 
 */
TEST_F(PairWiseKineticModelTest, density)
{   
    PairWiseKineticModelToy model ;

    /// 2+2 bins on both axes (-inf,0), [0,1), [0,2), [2,inf)
    model.setParameters(3, 0., 2., 2, 0.) ;

    // train model
    for(const auto& kinetic : kinetics)
    {   model.add(kinetic) ; }
    // at this moment the model contains, for IPD and PWD, for all 
    // POSn vs POSn+1 histogram
    //          (-inf,0), [0,1), [0,2), [2,inf)
    // [2,inf)      5       5      5       5
    // [0,2)        5       1      1       1
    // [0,1)        5       1      1       1
    // (-inf,0)     9       1      1       1
   
    ASSERT_EQ(model.isDensity(), false) ;

    // to density
    model.density() ;
    ASSERT_EQ(model.isDensity(), true) ;
    double n = 48. ;
    std::vector<std::vector<double>> m_expected({{9./n,5./n,5./n,5./n},
                                                 {5./n,1./n,1./n,1./n},
                                                 {5./n,1./n,1./n,1./n},
                                                 {5./n,1./n,1./n,1./n}}) ;
    // check histograms
    auto& index_pairs = model.getPairIndex() ;
    auto& histograms_ipd = model.getHistogramsIPD() ;
    auto& histograms_pwd = model.getHistogramsPWD() ;
    for(const auto& index_pair : index_pairs)
    {   const size_t& pos1 = index_pair.first ;
        const size_t& pos2 = index_pair.second ;

        std::vector<std::vector<double>> m_ipd = 
                ngsai::hist2d_to_matrix(histograms_ipd[pos1][pos2]) ;
        std::vector<std::vector<double>> m_pwd = 
                ngsai::hist2d_to_matrix(histograms_pwd[pos1][pos2]) ;
        
        for(size_t i=0; i<m_ipd.size(); i++)
        {   for(size_t j=0; j<m_ipd[i].size(); j++)
            {   ASSERT_DOUBLE_EQ(m_ipd[i][j], m_expected[i][j]) ;
                ASSERT_DOUBLE_EQ(m_pwd[i][j], m_expected[i][j]) ;
            }
        }
    }
}


/*!
 * \brief Tests that the log() method works properly. 
 */
TEST_F(PairWiseKineticModelTest, log)
{   
    PairWiseKineticModelToy model ;

    /// 2+2 bins on both axes (-inf,0), [0,1), [0,2), [2,inf)
    model.setParameters(3, 0., 2., 2, 0.) ;

    // train model
    for(const auto& kinetic : kinetics)
    {   model.add(kinetic) ; }
    // at this moment the model contains, for IPD and PWD, for all 
    // POSn vs POSn+1 histogram
    //          (-inf,0), [0,1), [0,2), [2,inf)
    // [2,inf)      5       5      5       5
    // [0,2)        5       1      1       1
    // [0,1)        5       1      1       1
    // (-inf,0)     9       1      1       1
   
    ASSERT_EQ(model.isLog(), false) ;

    // log
    model.log() ;
    ASSERT_EQ(model.isLog(), true) ;
    std::vector<std::vector<double>> m_expected({{log(9.),log(5.),log(5.),log(5.)},
                                                 {log(5.),log(1.),log(1.),log(1.)},
                                                 {log(5.),log(1.),log(1.),log(1.)},
                                                 {log(5.),log(1.),log(1.),log(1.)}}) ;
    // check histograms
    auto& index_pairs = model.getPairIndex() ;
    auto& histograms_ipd = model.getHistogramsIPD() ;
    auto& histograms_pwd = model.getHistogramsPWD() ;
    for(const auto& index_pair : index_pairs)
    {   const size_t& pos1 = index_pair.first ;
        const size_t& pos2 = index_pair.second ;

        std::vector<std::vector<double>> m_ipd = 
                ngsai::hist2d_to_matrix(histograms_ipd[pos1][pos2]) ;
        std::vector<std::vector<double>> m_pwd = 
                ngsai::hist2d_to_matrix(histograms_pwd[pos1][pos2]) ;
        
        for(size_t i=0; i<m_ipd.size(); i++)
        {   for(size_t j=0; j<m_ipd[i].size(); j++)
            {   ASSERT_DOUBLE_EQ(m_ipd[i][j], m_expected[i][j]) ;
                ASSERT_DOUBLE_EQ(m_pwd[i][j], m_expected[i][j]) ;
            }
        }
    }
}


/*!
 * \brief Tests that the exp() method works properly. 
 */
TEST_F(PairWiseKineticModelTest, exp)
{   
    PairWiseKineticModelToy model ;

    /// 2+2 bins on both axes (-inf,0), [0,1), [0,2), [2,inf)
    model.setParameters(3, 0., 2., 2, 0.) ;

    // train model
    for(const auto& kinetic : kinetics)
    {   model.add(kinetic) ; }
    // at this moment the model contains, for IPD and PWD, for all 
    // POSn vs POSn+1 histogram
    //          (-inf,0), [0,1), [0,2), [2,inf)
    // [2,inf)      5       5      5       5
    // [0,2)        5       1      1       1
    // [0,1)        5       1      1       1
    // (-inf,0)     9       1      1       1
   
    ASSERT_EQ(model.isLog(), false) ;
    model.log() ; 
    ASSERT_EQ(model.isLog(), true) ;
    model.exp() ;
    ASSERT_EQ(model.isLog(), false) ;

    // exp
    model.exp() ;
    ASSERT_EQ(model.isLog(), false) ;
    std::vector<std::vector<double>> m_expected({{exp(exp(log(9.))),exp(exp(log(5.))),exp(exp(log(5.))),exp(exp(log(5.)))},
                                                 {exp(exp(log(5.))),exp(exp(log(1.))),exp(exp(log(1.))),exp(exp(log(1.)))},
                                                 {exp(exp(log(5.))),exp(exp(log(1.))),exp(exp(log(1.))),exp(exp(log(1.)))},
                                                 {exp(exp(log(5.))),exp(exp(log(1.))),exp(exp(log(1.))),exp(exp(log(1.)))}}) ;
    // check histograms
    auto& index_pairs = model.getPairIndex() ;
    auto& histograms_ipd = model.getHistogramsIPD() ;
    auto& histograms_pwd = model.getHistogramsPWD() ;
    for(const auto& index_pair : index_pairs)
    {   const size_t& pos1 = index_pair.first ;
        const size_t& pos2 = index_pair.second ;

        std::vector<std::vector<double>> m_ipd = 
                ngsai::hist2d_to_matrix(histograms_ipd[pos1][pos2]) ;
        std::vector<std::vector<double>> m_pwd = 
                ngsai::hist2d_to_matrix(histograms_pwd[pos1][pos2]) ;
        
        for(size_t i=0; i<m_ipd.size(); i++)
        {   for(size_t j=0; j<m_ipd[i].size(); j++)
            {   ASSERT_DOUBLE_EQ(m_ipd[i][j], m_expected[i][j]) ;
                ASSERT_DOUBLE_EQ(m_pwd[i][j], m_expected[i][j]) ;
            }
        }
    }
}


/*!
 * \brief Tests that the += operator works properly. 
 */
TEST_F(PairWiseKineticModelTest, add_kineticmodel)
{   
    PairWiseKineticModelToy model1, model2, model3, model4 ;

    // must fail, both operande are not initialised
    EXPECT_THROW(model1.add(model1), std::invalid_argument) ;

    // 2+2 bins on both axes (-inf,0), [0,1), [0,2), [2,inf)
    model1.setParameters(3, 0., 2., 2, 0.) ;
    model2.setParameters(3, 0., 2., 2, 0.) ;
    model3.setParameters(3, 0., 2., 2, 0.) ;

    // must fail, model4 not initialised
    EXPECT_THROW(model1.add(model4), std::invalid_argument) ;

    // 2+2 bins (-inf,0), [0,1), [0,2), [2,inf)
    model4.setParameters(6, 0., 2., 2, 0.) ;

    // must fail, not same size
    EXPECT_THROW(model1.add(model4), std::invalid_argument) ;

    // 1+2 bins (-inf,-1), [-1,1), [1,inf)
    model4.setParameters(3, -1, 1, 1, 0.) ;

    // must fail, not same nb of bins
    EXPECT_THROW(model1.add(model4), std::runtime_error) ;


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
    // model2 += (model1 + model1) ;
    model2.add(model1) ;
    model2.add(model1) ;
    // model3 += (model1 + model1 + model1) ;
    model3.add(model1) ;
    model3.add(model1) ;
    model3.add(model1) ;


    // model2 is 2x model1
    auto& index_pairs = model2.getPairIndex() ;
    auto& histograms_ipd = model2.getHistogramsIPD() ;
    auto& histograms_pwd = model2.getHistogramsPWD() ;
    for(const auto& index_pair : index_pairs)
    {   const size_t& pos1 = index_pair.first ;
        const size_t& pos2 = index_pair.second ;
        // ipd
        std::vector<std::vector<double>> m_ipd = 
                ngsai::hist2d_to_matrix(histograms_ipd[pos1][pos2]) ;
        EXPECT_THAT(m_ipd, ContainerEq(m_expected*2.)) ;
        // pwd
        std::vector<std::vector<double>> m_pwd = 
                ngsai::hist2d_to_matrix(histograms_pwd[pos1][pos2]) ;
        EXPECT_THAT(m_ipd, ContainerEq(m_expected*2.)) ;
    }
    ASSERT_EQ(model2.getIsInit(), true) ;
    ASSERT_EQ(model2.getIsDensity(), false) ;
    ASSERT_EQ(model2.getIsLog(), false) ;
    // ASSERT_EQ(model2.getPairNb(), pair_exp.size()) ;
    EXPECT_THAT(model2.getPairIndex(), ContainerEq(pair_exp)) ;
    // model3 is 3x model1
    index_pairs = model3.getPairIndex() ;
    histograms_ipd = model3.getHistogramsIPD() ;
    histograms_pwd = model3.getHistogramsPWD() ;
    for(const auto& index_pair : index_pairs)
    {   const size_t& pos1 = index_pair.first ;
        const size_t& pos2 = index_pair.second ;
        // ipd
        std::vector<std::vector<double>> m_ipd = 
                ngsai::hist2d_to_matrix(histograms_ipd[pos1][pos2]) ;
        EXPECT_THAT(m_ipd, ContainerEq(m_expected*3.)) ;
        // pwd
        std::vector<std::vector<double>> m_pwd = 
                ngsai::hist2d_to_matrix(histograms_pwd[pos1][pos2]) ;
        EXPECT_THAT(m_ipd, ContainerEq(m_expected*3.)) ;
    }
    ASSERT_EQ(model3.getIsInit(), true) ;
    ASSERT_EQ(model3.getIsDensity(), false) ;
    ASSERT_EQ(model3.getIsLog(), false) ;
    ASSERT_EQ(model3.getPairNb(), pair_exp.size()) ;
    EXPECT_THAT(model3.getPairIndex(), ContainerEq(pair_exp)) ;


    
    // sum density models
    PairWiseKineticModelToy model5, model6 ;
    model5.setParameters(3, 0., 2., 2, 0.) ;
    model6.setParameters(3, 0., 2., 2, 0.) ;

    model1.density() ;
    // model5 += (model1 + model1) ;
    model5.add(model1) ;
    model5.add(model1) ;

    // model6 += (model1 + model1 + model1) ;
    model6.add(model1) ;
    model6.add(model1) ;
    model6.add(model1) ;

    // model5 is 2x model1
    index_pairs = model5.getPairIndex() ;
    histograms_ipd = model5.getHistogramsIPD() ;
    histograms_pwd = model5.getHistogramsPWD() ;
    for(const auto& index_pair : index_pairs)
    {   const size_t& pos1 = index_pair.first ;
        const size_t& pos2 = index_pair.second ;
        // ipd
        std::vector<std::vector<double>> m_ipd = 
                ngsai::hist2d_to_matrix(histograms_ipd[pos1][pos2]) ;
        EXPECT_THAT(m_ipd, ContainerEq((m_expected/48.)*2.)) ;
        // pwd
        std::vector<std::vector<double>> m_pwd = 
                ngsai::hist2d_to_matrix(histograms_pwd[pos1][pos2]) ;
        EXPECT_THAT(m_ipd, ContainerEq((m_expected/48.)*2.)) ;
    }
    ASSERT_EQ(model5.getIsInit(), true) ;
    ASSERT_EQ(model5.getIsDensity(), false) ;
    ASSERT_EQ(model5.getIsLog(), false) ;
    ASSERT_EQ(model5.getPairNb(), pair_exp.size()) ;
    EXPECT_THAT(model5.getPairIndex(), ContainerEq(pair_exp)) ;

    // model6 is 3x model1
    index_pairs = model6.getPairIndex() ;
    histograms_ipd = model6.getHistogramsIPD() ;
    histograms_pwd = model6.getHistogramsPWD() ;
    for(const auto& index_pair : index_pairs)
    {   const size_t& pos1 = index_pair.first ;
        const size_t& pos2 = index_pair.second ;
        // ipd
        std::vector<std::vector<double>> m_ipd = 
                ngsai::hist2d_to_matrix(histograms_ipd[pos1][pos2]) ;
        EXPECT_THAT(m_ipd, ContainerEq((m_expected/48.)*3.)) ;
        // pwd
        std::vector<std::vector<double>> m_pwd = 
                ngsai::hist2d_to_matrix(histograms_pwd[pos1][pos2]) ;
        EXPECT_THAT(m_ipd, ContainerEq((m_expected/48.)*3.)) ;
    }
    ASSERT_EQ(model6.getIsInit(), true) ;
    ASSERT_EQ(model6.getIsDensity(), false) ;
    ASSERT_EQ(model6.getIsLog(), false) ;
    ASSERT_EQ(model6.getPairNb(), pair_exp.size()) ;
    EXPECT_THAT(model6.getPairIndex(), ContainerEq(pair_exp)) ;
    
    
    // sum log density models
    PairWiseKineticModelToy model7, model8 ;
    model7.setParameters(3, 0., 2., 2, 0.) ;
    model8.setParameters(3, 0., 2., 2, 0.) ;
    
    model1.log() ;

    // model7 += (model1 + model1) ;
    model7.add(model1) ;
    model7.add(model1) ;

    // model8 += (model1 + model1 + model1) ;
    model8.add(model1) ;
    model8.add(model1) ;
    model8.add(model1) ;

    // model7 is 2x model1
    index_pairs = model7.getPairIndex() ;
    histograms_ipd = model7.getHistogramsIPD() ;
    histograms_pwd = model7.getHistogramsPWD() ;
    for(const auto& index_pair : index_pairs)
    {   const size_t& pos1 = index_pair.first ;
        const size_t& pos2 = index_pair.second ;
        // ipd
        std::vector<std::vector<double>> m_ipd = 
                ngsai::hist2d_to_matrix(histograms_ipd[pos1][pos2]) ;
        EXPECT_THAT(m_ipd, ContainerEq(log(m_expected/48.)*2.)) ;
        // pwd
        std::vector<std::vector<double>> m_pwd = 
                ngsai::hist2d_to_matrix(histograms_pwd[pos1][pos2]) ;
        EXPECT_THAT(m_ipd, ContainerEq(log(m_expected/48.)*2.)) ;
    }
    ASSERT_EQ(model7.getIsInit(), true) ;
    ASSERT_EQ(model7.getIsDensity(), false) ;
    ASSERT_EQ(model7.getIsLog(), false) ;
    ASSERT_EQ(model7.getPairNb(), pair_exp.size()) ;
    EXPECT_THAT(model7.getPairIndex(), ContainerEq(pair_exp)) ;

    // model8 is 3x model1
    index_pairs = model8.getPairIndex() ;
    histograms_ipd = model8.getHistogramsIPD() ;
    histograms_pwd = model8.getHistogramsPWD() ;
    for(const auto& index_pair : index_pairs)
    {   const size_t& pos1 = index_pair.first ;
        const size_t& pos2 = index_pair.second ;
        // ipd
        std::vector<std::vector<double>> m_ipd = 
                ngsai::hist2d_to_matrix(histograms_ipd[pos1][pos2]) ;
        EXPECT_THAT(m_ipd, ContainerEq(log(m_expected/48.)*3.)) ;
        // pwd
        std::vector<std::vector<double>> m_pwd = 
                ngsai::hist2d_to_matrix(histograms_pwd[pos1][pos2]) ;
        EXPECT_THAT(m_ipd, ContainerEq(log(m_expected/48.)*3.)) ;
    }
    ASSERT_EQ(model8.getIsInit(), true) ;
    ASSERT_EQ(model8.getIsDensity(), false) ;
    ASSERT_EQ(model8.getIsLog(), false) ;
    ASSERT_EQ(model8.getPairNb(), pair_exp.size()) ;
    EXPECT_THAT(model8.getPairIndex(), ContainerEq(pair_exp)) ;
}


/*!
 * \brief Tests that the logLikelihood() method works properly.
 */
TEST_F(PairWiseKineticModelTest, logLikelihood)
{   
    PairWiseKineticModelToy model ;

    // must fail not initialised
    EXPECT_THROW(model.logLikelihood(kinetics[0]), std::runtime_error) ;

    // 2+2 bins on both axes (-inf,0), [0,1), [0,2), [2,inf)
    model.setParameters(3, 0., 2., 2, 0.) ;

    // must fail not log density
    EXPECT_THROW(model.logLikelihood(kinetics[0]), std::runtime_error) ;

    // train model
    for(const auto& kinetic : kinetics)
    {   model.add(kinetic) ; }
    // at this moment the model contains, for IPD and PWD, for all 
    // POSn vs POSn+1 histogram
    //          (-inf,0), [0,1), [0,2), [2,inf)
    // [2,inf)      5       5      5       5
    // [0,2)        5       1      1       1
    // [0,1)        5       1      1       1
    // (-inf,0)     9       1      1       1

    model.density() ;
    model.log() ;

    double n = 48 ;
    ASSERT_DOUBLE_EQ(model.logLikelihood(kinetics[0]),
                     2.*(log(9/n) + log(9/n) + log(9/n))) ;
    ASSERT_DOUBLE_EQ(model.logLikelihood(kinetics[1]),
                     2.*(log(5/n) + log(5/n) + log(9/n))) ;
    ASSERT_DOUBLE_EQ(model.logLikelihood(kinetics[2]),
                     2.*(log(5/n) + log(5/n) + log(9/n))) ;
    ASSERT_DOUBLE_EQ(model.logLikelihood(kinetics[3]),
                     2.*(log(5/n) + log(5/n) + log(9/n))) ;
    ASSERT_DOUBLE_EQ(model.logLikelihood(kinetics[4]),
                     2.*(log(5/n) + log(9/n) + log(5/n))) ;
    ASSERT_DOUBLE_EQ(model.logLikelihood(kinetics[5]),
                     2.*(log(1/n) + log(5/n) + log(5/n))) ;
    ASSERT_DOUBLE_EQ(model.logLikelihood(kinetics[6]),
                     2.*(log(1/n) + log(5/n) + log(5/n))) ;
    ASSERT_DOUBLE_EQ(model.logLikelihood(kinetics[7]),
                     2.*(log(1/n) + log(5/n) + log(5/n))) ;
    ASSERT_DOUBLE_EQ(model.logLikelihood(kinetics[8]),
                     2.*(log(5/n) + log(9/n) + log(5/n))) ;
    ASSERT_DOUBLE_EQ(model.logLikelihood(kinetics[9]),
                     2.*(log(1/n) + log(5/n) + log(5/n))) ;
    ASSERT_DOUBLE_EQ(model.logLikelihood(kinetics[10]),
                     2.*(log(1/n) + log(5/n) + log(5/n))) ;
    ASSERT_DOUBLE_EQ(model.logLikelihood(kinetics[11]),
                     2.*(log(1/n) + log(5/n) + log(5/n))) ;
    ASSERT_DOUBLE_EQ(model.logLikelihood(kinetics[12]),
                     2.*(log(5/n) + log(9/n) + log(5/n))) ;
    ASSERT_DOUBLE_EQ(model.logLikelihood(kinetics[13]),
                     2.*(log(1/n) + log(5/n) + log(5/n))) ;
    ASSERT_DOUBLE_EQ(model.logLikelihood(kinetics[14]),
                     2.*(log(1/n) + log(5/n) + log(5/n))) ;


    ASSERT_DOUBLE_EQ(model.logLikelihood(kinetics[16]),
                     2.*(log(9/n) + log(9/n) + log(9/n))) ;
    ASSERT_DOUBLE_EQ(model.logLikelihood(kinetics[17]),
                     2.*(log(5/n) + log(5/n) + log(9/n))) ;
    ASSERT_DOUBLE_EQ(model.logLikelihood(kinetics[18]),
                     2.*(log(5/n) + log(5/n) + log(9/n))) ;
    ASSERT_DOUBLE_EQ(model.logLikelihood(kinetics[19]),
                     2.*(log(5/n) + log(5/n) + log(9/n))) ;
    ASSERT_DOUBLE_EQ(model.logLikelihood(kinetics[20]),
                     2.*(log(9/n) + log(5/n) + log(5/n))) ;
    ASSERT_DOUBLE_EQ(model.logLikelihood(kinetics[21]),
                     2.*(log(5/n) + log(1/n) + log(5/n))) ;
    ASSERT_DOUBLE_EQ(model.logLikelihood(kinetics[22]),
                     2.*(log(5/n) + log(1/n) + log(5/n))) ;
    ASSERT_DOUBLE_EQ(model.logLikelihood(kinetics[23]),
                     2.*(log(5/n) + log(1/n) + log(5/n))) ;
    ASSERT_DOUBLE_EQ(model.logLikelihood(kinetics[24]),
                     2.*(log(9/n) + log(5/n) + log(5/n))) ;
    ASSERT_DOUBLE_EQ(model.logLikelihood(kinetics[25]),
                     2.*(log(5/n) + log(1/n) + log(5/n))) ;
    ASSERT_DOUBLE_EQ(model.logLikelihood(kinetics[26]),
                     2.*(log(5/n) + log(1/n) + log(5/n))) ;
    ASSERT_DOUBLE_EQ(model.logLikelihood(kinetics[27]),
                     2.*(log(5/n) + log(1/n) + log(5/n))) ;
    ASSERT_DOUBLE_EQ(model.logLikelihood(kinetics[28]),
                     2.*(log(9/n) + log(5/n) + log(5/n))) ;
    ASSERT_DOUBLE_EQ(model.logLikelihood(kinetics[29]),
                     2.*(log(5/n) + log(1/n) + log(5/n))) ;
    ASSERT_DOUBLE_EQ(model.logLikelihood(kinetics[30]),
                     2.*(log(5/n) + log(1/n) + log(5/n))) ;
    ASSERT_DOUBLE_EQ(model.logLikelihood(kinetics[31]),
                     2.*(log(5/n) + log(1/n) + log(5/n))) ;

    ASSERT_DOUBLE_EQ(model.logLikelihood(kinetics[32]),
                     2.*(log(9/n) + log(9/n) + log(9/n))) ;
    ASSERT_DOUBLE_EQ(model.logLikelihood(kinetics[33]),
                     2.*(log(5/n) + log(9/n) + log(5/n))) ;
    ASSERT_DOUBLE_EQ(model.logLikelihood(kinetics[34]),
                     2.*(log(5/n) + log(9/n) + log(5/n))) ;
    ASSERT_DOUBLE_EQ(model.logLikelihood(kinetics[35]),
                     2.*(log(5/n) + log(9/n) + log(5/n))) ;
    ASSERT_DOUBLE_EQ(model.logLikelihood(kinetics[36]),
                     2.*(log(9/n) + log(5/n) + log(5/n))) ;
    ASSERT_DOUBLE_EQ(model.logLikelihood(kinetics[37]),
                     2.*(log(5/n) + log(5/n) + log(1/n))) ;
    ASSERT_DOUBLE_EQ(model.logLikelihood(kinetics[38]),
                     2.*(log(5/n) + log(5/n) + log(1/n))) ;
    ASSERT_DOUBLE_EQ(model.logLikelihood(kinetics[39]),
                     2.*(log(5/n) + log(5/n) + log(1/n))) ;
    ASSERT_DOUBLE_EQ(model.logLikelihood(kinetics[40]),
                     2.*(log(9/n) + log(5/n) + log(5/n))) ;
    ASSERT_DOUBLE_EQ(model.logLikelihood(kinetics[41]),
                     2.*(log(5/n) + log(5/n) + log(1/n))) ;
    ASSERT_DOUBLE_EQ(model.logLikelihood(kinetics[42]),
                     2.*(log(5/n) + log(5/n) + log(1/n))) ;
    ASSERT_DOUBLE_EQ(model.logLikelihood(kinetics[43]),
                     2.*(log(5/n) + log(5/n) + log(1/n))) ;
    ASSERT_DOUBLE_EQ(model.logLikelihood(kinetics[44]),
                     2.*(log(9/n) + log(5/n) + log(5/n))) ;
    ASSERT_DOUBLE_EQ(model.logLikelihood(kinetics[45]),
                     2.*(log(5/n) + log(5/n) + log(1/n))) ;
    ASSERT_DOUBLE_EQ(model.logLikelihood(kinetics[46]),
                     2.*(log(5/n) + log(5/n) + log(1/n))) ;
    ASSERT_DOUBLE_EQ(model.logLikelihood(kinetics[47]),
                     2.*(log(5/n) + log(5/n) + log(1/n))) ;

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