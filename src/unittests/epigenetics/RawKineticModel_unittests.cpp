#include <gtest/gtest.h>
#include <gmock/gmock.h>
#include <string>
#include <vector>
#include <stdexcept>
#include <limits>
#include <cmath>

#include <ngsaipp/epigenetics/RawKineticModel.hpp>
#include <ngsaipp/epigenetics/KineticModel.hpp>
#include <ngsaipp/epigenetics/NormalizedKineticModel.hpp>
#include <ngsaipp/epigenetics/PairWiseKineticModel.hpp>
#include <ngsaipp/epigenetics/PairWiseNormalizedKineticModel.hpp>
#include <ngsaipp/epigenetics/DiPositionKineticModel.hpp>
#include <ngsaipp/epigenetics/DiPositionNormalizedKineticModel.hpp>


using testing::ElementsAre ;
using testing::ContainerEq ;


/*
 * This file ontains tests for the ngsai::RawKineticModel class from 
 * src/epigenetics/RawKineticModel.cpp
 */


/*!
 * \brief A child class of RawKineticModel allowing to access the protected 
 * fields of RawKineticModel for testing purposes.
 */
class RawKineticModelToy : public ngsai::RawKineticModel
{
    public:

        RawKineticModelToy()
            : RawKineticModel()
        { ; }

        RawKineticModelToy(const ngsai::RawKineticModel& other)
            : RawKineticModel(other)
        { ; }

        RawKineticModelToy(ngsai::RawKineticModel&& other)
            : RawKineticModel(other)
        { ; }

        ~RawKineticModelToy()
        { ; }

        // ugly workaround to test RawKineticModel.copy() using a
        // RawKineticModelToy instance
        RawKineticModelToy* copy() const
        {   RawKineticModel* copy = RawKineticModel::copy() ;
            RawKineticModelToy* copy2 = new RawKineticModelToy(*copy) ;
            delete copy ;
            return copy2 ;
        }

        std::vector<hist_1d_double>& getHistogramsIPD()
        {   return m_histograms_ipd ; }

        std::vector<hist_1d_double>& getHistogramsPWD()
        {   return m_histograms_pwd ; }

        bool& getIsInit()
        {   return m_is_init ; }

        bool& getIsDensity()
        {   return m_is_density ; }

        bool& getIsLog()
        {   return m_is_log ; }

} ;



/*!
 * \brief Tests that constructor works properly.
 */
TEST(RawKineticModelTest, constructor)
{
    RawKineticModelToy model ;

    ASSERT_EQ(model.getHistogramsIPD().size(), 0) ;
    ASSERT_EQ(model.getHistogramsPWD().size(), 0) ;
    ASSERT_EQ(model.isInit(), false) ;
    ASSERT_EQ(model.isDensity(), false) ;
    ASSERT_EQ(model.isLog(), false) ;
}


/*!
 * \brief Tests that isInit() method works properly.
 */
TEST(RawKineticModelTest, isInit)
{
    RawKineticModelToy model ;

    ASSERT_EQ(model.isInit(), false) ;
    bool& init = model.getIsInit() ;
    init = true ;
    ASSERT_EQ(model.isInit(), true) ;
}


/*!
 * \brief Tests that isDensity() method works properly.
 */
TEST(RawKineticModelTest, isDensity)
{
    RawKineticModelToy model ;

    ASSERT_EQ(model.isDensity(), false) ;
    bool& density = model.getIsDensity() ;
    density = true ;
    ASSERT_EQ(model.isDensity(), true) ;
}


/*!
 * \brief Tests that the setParameters() method works properly and that 
 * everything is initialised properly or throws the expected exceptions.
 */
TEST(RawKineticModelTest, setParameters)
{   
    size_t size(5) ;
    double min(-2.) ;
    double max(2.) ;
    size_t n_bin(4) ;
    double pseudocounts(10.) ;

    // it must work properly
    RawKineticModelToy model ;
    model.setParameters(size,
                        min,
                        max,
                        n_bin,
                        pseudocounts) ;
    ASSERT_EQ(model.getHistogramsIPD().size(), size) ;
    ASSERT_EQ(model.getHistogramsPWD().size(), size) ;
    ASSERT_EQ(model.isInit(), true) ;
    ASSERT_EQ(model.isDensity(), false) ;
    ASSERT_EQ(model.isLog(), false) ;
    double inf = std::numeric_limits<double>::infinity() ;
    std::vector<std::pair<double,double>> bins_exp = 
                                            {std::make_pair(-inf, -2),
                                             std::make_pair(-2., -1.),
                                             std::make_pair(-1.,  0.),
                                             std::make_pair( 0.,  1.),
                                             std::make_pair( 1.,  2.),
                                             std::make_pair( 2.,  inf)} ;
    auto& histograms_ipd = model.getHistogramsIPD() ;
    for(auto& histogram : histograms_ipd)
    {   // check bin content and boundaries
        std::vector<std::pair<double,double>> bins ;
        for(const auto& x : bh::indexed(histogram, bh::coverage::all))
        {   bins.push_back(std::make_pair(x.bin().lower(),
                                          x.bin().upper())) ;
            ASSERT_EQ(*x, pseudocounts) ;
        }
        EXPECT_THAT(bins, ContainerEq(bins_exp)) ;
    }
    auto& histograms_pwd = model.getHistogramsPWD() ;
    for(auto& histogram : histograms_pwd)
    {   // check bin content and boundaries
        std::vector<std::pair<double,double>> bins ;
        for(const auto& x : bh::indexed(histogram, bh::coverage::all))
        {   bins.push_back(std::make_pair(x.bin().lower(),
                                          x.bin().upper())) ;
            ASSERT_EQ(*x, pseudocounts) ;
        }
        EXPECT_THAT(bins, ContainerEq(bins_exp)) ;
    }

    // it must fail, size == 0
    EXPECT_THROW(model.setParameters(0, min, max, n_bin ,pseudocounts),
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
TEST(RawKineticModelTest, size)
{   // before initialisation
    RawKineticModelToy model ;
    ASSERT_EQ(model.size(), 0) ;

    // 5bp
    model.setParameters(5, -2, 2, 4, 10) ;
    ASSERT_EQ(model.size(), 5) ;

    // 10bp
    model.setParameters(10, -2, 2, 4, 10) ;
    ASSERT_EQ(model.size(), 10) ;

    // 1bp
    model.setParameters(1, -2, 2, 4, 10) ;
    ASSERT_EQ(model.size(), 1) ;
}


/*!
 * \brief Tests that getKineticSignalRequiredSize() method works properly.
 */
TEST(RawKineticModelTest, getKineticSignalRequiredSize)
{   
    size_t size(5) ;
    double min(-2.) ;
    double max(2.) ;
    size_t n_bin(4) ;
    double pseudocounts(10.) ;

    RawKineticModelToy model ;
    
    // uninitialised
    ASSERT_EQ(model.getKineticSignalRequiredSize(), 0) ;
    
    // initialised
    model.setParameters(size,
                        min,
                        max,
                        n_bin,
                        pseudocounts) ;
    ASSERT_EQ(model.getKineticSignalRequiredSize(), size) ;
}


/*!
 * \brief Tests that the getBinBoundaries() method works properly.
 */
TEST(RawKineticModelTest, getBinBoundaries)
{
    RawKineticModelToy model ;

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
 * \brief Tests that the add(double, size_t) method works properly.
 */
TEST(RawKineticModelTest, add_value)
{
    double pseudocounts(2.) ;

    RawKineticModelToy model ;

    // model is not initialised yet, must fail
    EXPECT_THROW(model.add(1., 1., 0), std::runtime_error) ;
    EXPECT_THROW(model.add(1., 1., 1), std::runtime_error) ;
    EXPECT_THROW(model.add(1., 1., 2), std::runtime_error) ;
    EXPECT_THROW(model.add(1., 1., 3), std::runtime_error) ;
    EXPECT_THROW(model.add(1., 1., 4), std::runtime_error) ;

    // 4+2 bins [-inf,-2), [-2,-1), [-1,0), [0,1), [1,2), [2,inf)
    model.setParameters(5, -2, 2, 4, pseudocounts) ;

    // each value goes to a different bin
    double inf = std::numeric_limits<double>::infinity() ;
    std::vector<double> values({-inf, -2., -1., 0, 1., inf}) ;

    // add +1 count to each bin at pos 0
    model.add(values[0], values[0], 0) ;
    model.add(values[1], values[1], 0) ;
    model.add(values[2], values[2], 0) ;
    model.add(values[3], values[3], 0) ;
    model.add(values[4], values[4], 0) ;
    model.add(values[5], values[5], 0) ;

    // add +1 count to each bin at pos 1
    model.add(values[0], values[0], 1) ;
    model.add(values[1], values[1], 1) ;
    model.add(values[2], values[2], 1) ;
    model.add(values[3], values[3], 1) ;
    model.add(values[4], values[4], 1) ;
    model.add(values[5], values[5], 1) ;

    // add +1 count to each bin at pos 2
    model.add(values[0], values[0], 2) ;
    model.add(values[1], values[1], 2) ;
    model.add(values[2], values[2], 2) ;
    model.add(values[3], values[3], 2) ;
    model.add(values[4], values[4], 2) ;
    model.add(values[5], values[5], 2) ;

    // add +1 count to each bin at pos 3
    model.add(values[0], values[0], 3) ;
    model.add(values[1], values[1], 3) ;
    model.add(values[2], values[2], 3) ;
    model.add(values[3], values[3], 3) ;
    model.add(values[4], values[4], 3) ;
    model.add(values[5], values[5], 3) ;
    
    // add +1 count to each bin at pos 4
    model.add(values[0], values[0], 4) ;
    model.add(values[1], values[1], 4) ;
    model.add(values[2], values[2], 4) ;
    model.add(values[3], values[3], 4) ;
    model.add(values[4], values[4], 4) ;
    model.add(values[5], values[5], 4) ;

    auto& histograms_ipd = model.getHistogramsIPD() ;
    for(auto& histogram : histograms_ipd)
    {   for(const auto& x : bh::indexed(histogram, bh::coverage::all))
        {   ASSERT_EQ(*x, pseudocounts + 1) ; }
    }

    auto& histograms_pwd = model.getHistogramsPWD() ;
    for(auto& histogram : histograms_pwd)
    {   for(const auto& x : bh::indexed(histogram, bh::coverage::all))
        {   ASSERT_EQ(*x, pseudocounts + 1) ; }
    }
}


/*!
 * \brief Tests that the add() method with KineticSignal works properly.
 */
TEST(RawKineticModelTest, add_kineticsignal)
{
    double pseudocounts(2.) ;

    RawKineticModelToy model ;

    double inf = std::numeric_limits<double>::infinity() ;

    std::string seq("AAAAA") ;

    // add +1 count in [-inf,-2) bin at each position
    std::vector<double> val_1({-inf, -inf, -inf, -inf, -inf}) ;

    // add +1 count in [-2,-1) bin at each position
    std::vector<double> val_2({ -2.,  -2.,  -2.,  -2.,  -2.}) ;

    // add +1 count in [-1,0) bin at each position
    std::vector<double> val_3({ -1.,  -1.,  -1.,  -1.,  -1.}) ;

    // add +1 count in [0,1) bin at each position
    std::vector<double> val_4({  0.,   0.,   0.,   0.,   0.}) ;

    // add +1 count in [1,2) bin at each position
    std::vector<double> val_5({  1.,   1.,   1.,   1.,   1.}) ;

    // add +1 count in [2, inf) bin at each position
    std::vector<double> val_6({inf, inf, inf, inf, inf}) ;

    ngsai::KineticSignal k1(seq, seq, val_1, val_1, val_1, val_1) ;
    ngsai::KineticSignal k2(seq, seq, val_2, val_2, val_2, val_2) ;
    ngsai::KineticSignal k3(seq, seq, val_3, val_3, val_3, val_3) ;
    ngsai::KineticSignal k4(seq, seq, val_4, val_4, val_4, val_4) ;
    ngsai::KineticSignal k5(seq, seq, val_5, val_5, val_5, val_5) ;
    ngsai::KineticSignal k6(seq, seq, val_6, val_6, val_6, val_6) ;

    // model is not initialised yet, must fail
    EXPECT_THROW(model.add(k2), std::runtime_error) ;

    // 4+2 bins (-inf,-2), [-2,-1), [-1,0), [0,1), [1,2), [2,inf)
    model.setParameters(5, -2, 2, 4, pseudocounts) ;

    // add +1 count to each bin at each pos
    model.add(k1) ;
    model.add(k2) ;
    model.add(k3) ;
    model.add(k4) ;
    model.add(k5) ;
    model.add(k6) ;

    auto& histograms_ipd = model.getHistogramsIPD() ;
    for(auto& histogram : histograms_ipd)
    {   for(const auto& x : bh::indexed(histogram, bh::coverage::all))
        {   ASSERT_EQ(*x, pseudocounts + 1) ; }
    }

    auto& histograms_pwd = model.getHistogramsPWD() ;
    for(auto& histogram : histograms_pwd)
    {   for(const auto& x : bh::indexed(histogram, bh::coverage::all))
        {   ASSERT_EQ(*x, pseudocounts + 1) ; }
    }

    std::vector<double> val_e1({1.}) ;
    std::vector<double> val_e2({1.,2.}) ;
    std::vector<double> val_e3({1.,2.,3.}) ;
    std::vector<double> val_e4({1.,2.,3.,4.}) ;
    std::vector<double> val_e6({1.,2.,3.,4.,5.,6.}) ;
    std::string seq_e1("A") ;
    std::string seq_e2("AA") ;
    std::string seq_e3("AAA") ;
    std::string seq_e4("AAAA") ;
    std::string seq_e6("AAAAAA") ;
    ngsai::KineticSignal k_err0 ;
    ngsai::KineticSignal k_err1(seq_e1, seq_e1, val_e1, val_e1, val_e1, val_e1) ;
    ngsai::KineticSignal k_err2(seq_e2, seq_e2, val_e2, val_e2, val_e2, val_e2) ;
    ngsai::KineticSignal k_err3(seq_e3, seq_e3, val_e3, val_e3, val_e3, val_e3) ;
    ngsai::KineticSignal k_err4(seq_e4, seq_e4, val_e4, val_e4, val_e4, val_e4) ;
    ngsai::KineticSignal k_err6(seq_e6, seq_e6, val_e6, val_e6, val_e6, val_e6) ;

    // must fail because not same size as model
    EXPECT_THROW(model.add(k_err0), std::invalid_argument) ; // empty
    EXPECT_THROW(model.add(k_err1), std::invalid_argument) ; // too short
    EXPECT_THROW(model.add(k_err2), std::invalid_argument) ; // too short
    EXPECT_THROW(model.add(k_err3), std::invalid_argument) ; // too short
    EXPECT_THROW(model.add(k_err4), std::invalid_argument) ; // too short
    EXPECT_THROW(model.add(k_err6), std::invalid_argument) ; // too long
}


/*!
 * \brief Tests that the copy constructor works properly.
 */
TEST(RawKineticModelTest, constructor_copy)
{
    // copy empty model
    RawKineticModelToy model1 ;
    RawKineticModelToy model2(model1) ;
    ASSERT_EQ(model1.getHistogramsIPD().size(), 
              model2.getHistogramsIPD().size()) ;
    ASSERT_EQ(model1.getHistogramsPWD().size(), 
              model2.getHistogramsPWD().size()) ;
    ASSERT_EQ(model1.isInit(),    model2.isInit()) ;
    ASSERT_EQ(model1.isDensity(), model2.isDensity()) ;
    ASSERT_EQ(model1.isLog(),     model2.isLog()) ;

    // copy model with data
    // [-inf,-1), [-1,0),  [0,1), [1,inf)
    model1.setParameters(1, -1, 1, 2, 1) ;
    model1.add(0., 0., 0) ; // add +1 to bin [0,1) at each pos

    RawKineticModelToy model3(model1) ;

    ASSERT_EQ(model1.getHistogramsIPD().size(), 
              model3.getHistogramsIPD().size()) ;
    ASSERT_EQ(model1.getHistogramsPWD().size(), 
              model3.getHistogramsPWD().size()) ;
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
TEST(RawKineticModelTest, constructor_move)
{
    // assign empty model
    RawKineticModelToy model1 ;
    RawKineticModelToy model2(model1) ;
    RawKineticModelToy model3(std::move(model2)) ;
    ASSERT_EQ(model1.getHistogramsIPD().size(), 
              model3.getHistogramsPWD().size()) ;
    ASSERT_EQ(model1.isInit(),    model3.isInit()) ;
    ASSERT_EQ(model1.isDensity(), model3.isDensity()) ;
    ASSERT_EQ(model1.isLog(),     model3.isLog()) ;

    // copy model with data
    // [-inf,-1), [-1,0),  [0,1), [1,inf)
    model1.setParameters(1, -1, 1, 2, 1) ;
    model1.add(0., 0., 0) ; // add +1 to bin [0,1)

    RawKineticModelToy model4(model1) ;
    RawKineticModelToy model5(std::move(model4)) ;

    ASSERT_EQ(model1.getHistogramsIPD().size(),
              model5.getHistogramsIPD().size()) ;
    ASSERT_EQ(model1.getHistogramsPWD().size(),
              model5.getHistogramsPWD().size()) ;
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
TEST(RawKineticModelTest, assignment_operator)
{
    // assign empty model
    RawKineticModelToy model1 ;
    RawKineticModelToy model2 = model1 ;
    ASSERT_EQ(model1.getHistogramsIPD().size(),
              model2.getHistogramsIPD().size()) ;
    ASSERT_EQ(model1.getHistogramsPWD().size(),
              model2.getHistogramsPWD().size()) ;
    ASSERT_EQ(model1.isInit(),    model2.isInit()) ;
    ASSERT_EQ(model1.isDensity(), model2.isDensity()) ;
    ASSERT_EQ(model1.isLog(),     model2.isLog()) ;

    // assign model with data
    // [-inf,-1), [-1,0),  [0,1), [1,inf)
    model1.setParameters(1, -1, 1, 2, 1) ;
    model1.add(0., 0., 0) ; // add +1 to bin [0,1)

    model2 = model1 ;

    ASSERT_EQ(model1.getHistogramsIPD().size(),
              model2.getHistogramsIPD().size()) ;
    ASSERT_EQ(model1.getHistogramsPWD().size(),
              model2.getHistogramsPWD().size()) ;
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
TEST(RawKineticModelTest, assignment_move_operator)
{
    // assign empty model
    RawKineticModelToy model1 ;
    RawKineticModelToy model2 = model1 ;
    RawKineticModelToy model3(std::move(model2)) ;
    ASSERT_EQ(model1.getHistogramsIPD().size(), 
              model3.getHistogramsIPD().size()) ;
    ASSERT_EQ(model1.getHistogramsPWD().size(), 
              model3.getHistogramsPWD().size()) ;
    ASSERT_EQ(model1.isInit(),    model3.isInit()) ;
    ASSERT_EQ(model1.isDensity(), model3.isDensity()) ;
    ASSERT_EQ(model1.isLog(),     model3.isLog()) ;

    // assign model with data
    // [-inf,-1), [-1,0),  [0,1), [1,inf)
    model1.setParameters(1, -1, 1, 2, 1) ;
    model1.add(0., 0., 0) ; // add +1 to bin [0,1)

    model3 = model1 ;

    RawKineticModelToy model4(std::move(model3)) ;

    ASSERT_EQ(model1.getHistogramsIPD().size(), 
              model4.getHistogramsIPD().size()) ;
    ASSERT_EQ(model1.getHistogramsPWD().size(), 
              model4.getHistogramsPWD().size()) ;
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
TEST(RawKineticModelTest, copy)
{
    // copy empty model
    RawKineticModelToy model1 ;
    RawKineticModelToy* copy(nullptr);
    copy = model1.copy() ;
    ASSERT_EQ(model1.getHistogramsIPD().size(), 
              copy->getHistogramsIPD().size()) ;
    ASSERT_EQ(model1.getHistogramsPWD().size(), 
              copy->getHistogramsPWD().size()) ;
    ASSERT_EQ(model1.isInit(),    copy->isInit()) ;
    ASSERT_EQ(model1.isDensity(), copy->isDensity()) ;
    ASSERT_EQ(model1.isLog(),     copy->isLog()) ;
    delete copy ;
    copy = nullptr ;

    // copy model with data
    // [-inf,-1), [-1,0),  [0,1), [1,inf)
    model1.setParameters(1, -1, 1, 2, 1) ;
    model1.add(0., 0., 0) ; // add +1 to bin [0,1)

    copy = model1.copy() ;

    ASSERT_EQ(model1.getHistogramsIPD().size(), 
              copy->getHistogramsIPD().size()) ;
    ASSERT_EQ(model1.getHistogramsPWD().size(), 
              copy->getHistogramsPWD().size()) ;
    ASSERT_EQ(model1.isInit(),    copy->isInit()) ;
    ASSERT_EQ(model1.isDensity(), copy->isDensity()) ;
    ASSERT_EQ(model1.isLog(),     copy->isLog()) ;
    auto& histograms1_ipd = model1.getHistogramsIPD() ;
    auto& histograms1_pwd = model1.getHistogramsPWD() ;
    auto& histogramsc_ipd = copy->getHistogramsIPD() ;
    auto& histogramsc_pwd = copy->getHistogramsPWD() ;
    EXPECT_THAT(histograms1_ipd, ContainerEq(histogramsc_ipd)) ;
    EXPECT_THAT(histograms1_pwd, ContainerEq(histogramsc_pwd)) ;
}


/*!
 * \brief Tests taht the density() method works properly. 
 */
TEST(RawKineticModelTest, density)
{   
    RawKineticModelToy model ;

    // 3+2 bins [-inf,0), [0,1), [1,2), [2,3), [3,inf)
    model.setParameters(3, 0., 3., 3, 0.) ;

    std::string seq("AAA") ;

    std::vector<double> val_1({-1., -1., -1.}) ;
    std::vector<double> val_2({ 0.,  0.,  0.}) ;
    std::vector<double> val_3({ 1.,  1.,  1.}) ;
    std::vector<double> val_4({ 1.,  1.,  1.}) ;
    std::vector<double> val_5({ 2.,  2.,  2.}) ;
    std::vector<double> val_6({ 4.,  4.,  4.}) ;

    ngsai::KineticSignal k1(seq, seq, val_1, val_1, val_1, val_1) ;
    ngsai::KineticSignal k2(seq, seq, val_2, val_2, val_2, val_2) ;
    ngsai::KineticSignal k3(seq, seq, val_3, val_3, val_3, val_3) ;
    ngsai::KineticSignal k4(seq, seq, val_4, val_4, val_4, val_4) ;
    ngsai::KineticSignal k5(seq, seq, val_5, val_5, val_5, val_5) ;
    ngsai::KineticSignal k6(seq, seq, val_6, val_6, val_6, val_6) ;

    model.add(k1) ;
    model.add(k2) ;
    model.add(k3) ;
    model.add(k4) ;
    model.add(k5) ;
    model.add(k6) ;

    ASSERT_EQ(model.isDensity(), false) ;

    model.density() ;
    std::vector<double> density_exp({1./6.,1./6.,2./6.,1./6.,1./6.}) ;

    // IPD
    ASSERT_EQ(model.isDensity(), true) ;
    auto& histograms_ipd = model.getHistogramsIPD() ;
    for(auto& histogram : histograms_ipd)
    {   size_t i=0 ;
        for(const auto& x : bh::indexed(histogram, bh::coverage::all))
        {   ASSERT_EQ(*x, density_exp[i]) ;
            i++ ;
        }
    }

    // PWD
    ASSERT_EQ(model.isDensity(), true) ;
    auto& histograms_pwd = model.getHistogramsPWD() ;
    for(auto& histogram : histograms_pwd)
    {   size_t i=0 ;
        for(const auto& x : bh::indexed(histogram, bh::coverage::all))
        {   ASSERT_EQ(*x, density_exp[i]) ;
            i++ ;
        }
    }
}


/*!
 * \brief Tests that the log() method works properly. 
 */
TEST(RawKineticModelTest, log)
{
    RawKineticModelToy model ;

    // 3+2 bins [-inf,0), [0,1), [1,2), [2,3), [3,inf)
    model.setParameters(3, 0., 3., 3, 0.) ;

    std::string seq("AAA") ;

    std::vector<double> val_1({-1., -1., -1.}) ;
    std::vector<double> val_2({ 0.,  0.,  0.}) ;
    std::vector<double> val_3({ 1.,  1.,  1.}) ;
    std::vector<double> val_4({ 1.,  1.,  1.}) ;
    std::vector<double> val_5({ 2.,  2.,  2.}) ;
    std::vector<double> val_6({ 4.,  4.,  4.}) ;

    ngsai::KineticSignal k1(seq, seq, val_1, val_1, val_1, val_1) ;
    ngsai::KineticSignal k2(seq, seq, val_2, val_2, val_2, val_2) ;
    ngsai::KineticSignal k3(seq, seq, val_3, val_3, val_3, val_3) ;
    ngsai::KineticSignal k4(seq, seq, val_4, val_4, val_4, val_4) ;
    ngsai::KineticSignal k5(seq, seq, val_5, val_5, val_5, val_5) ;
    ngsai::KineticSignal k6(seq, seq, val_6, val_6, val_6, val_6) ;

    model.add(k1) ;
    model.add(k2) ;
    model.add(k3) ;
    model.add(k4) ;
    model.add(k5) ;
    model.add(k6) ;

    ASSERT_EQ(model.isLog(), false) ;

    model.log() ;
    std::vector<double> log_exp({log(1.),
                                 log(1.),
                                 log(2.),
                                 log(1.),
                                 log(1.)}) ;

    // IPD
    ASSERT_EQ(model.isLog(), true) ;
    auto& histograms_ipd = model.getHistogramsIPD() ;
    for(auto& histogram : histograms_ipd)
    {   size_t i=0 ;
        for(const auto& x : bh::indexed(histogram, bh::coverage::all))
        {   ASSERT_EQ(*x, log_exp[i]) ;
            i++ ;
        }
    }

    // PWD
    ASSERT_EQ(model.isLog(), true) ;
    auto& histograms_pwd = model.getHistogramsPWD() ;
    for(auto& histogram : histograms_pwd)
    {   size_t i=0 ;
        for(const auto& x : bh::indexed(histogram, bh::coverage::all))
        {   ASSERT_EQ(*x, log_exp[i]) ;
            i++ ;
        }
    }
}


/*!
 * \brief Tests that the exp() method works properly. 
 */
TEST(RawKineticModelTest, exp)
{   
    RawKineticModelToy model ;

    // 3+2 bins [-inf,0), [0,1), [1,2), [2,3), [3,inf)
    model.setParameters(3, 0., 3., 3, 0.) ;

    std::string seq("AAA") ;

    std::vector<double> val_1({-1., -1., -1.}) ;
    std::vector<double> val_2({ 0.,  0.,  0.}) ;
    std::vector<double> val_3({ 1.,  1.,  1.}) ;
    std::vector<double> val_4({ 1.,  1.,  1.}) ;
    std::vector<double> val_5({ 2.,  2.,  2.}) ;
    std::vector<double> val_6({ 4.,  4.,  4.}) ;

    ngsai::KineticSignal k1(seq, seq, val_1, val_1, val_1, val_1) ;
    ngsai::KineticSignal k2(seq, seq, val_2, val_2, val_2, val_2) ;
    ngsai::KineticSignal k3(seq, seq, val_3, val_3, val_3, val_3) ;
    ngsai::KineticSignal k4(seq, seq, val_4, val_4, val_4, val_4) ;
    ngsai::KineticSignal k5(seq, seq, val_5, val_5, val_5, val_5) ;
    ngsai::KineticSignal k6(seq, seq, val_6, val_6, val_6, val_6) ;

    model.add(k1) ;
    model.add(k2) ;
    model.add(k3) ;
    model.add(k4) ;
    model.add(k5) ;
    model.add(k6) ;

    // converts to log and back to exp to check isLog
    // this is equal to do nothing
    ASSERT_EQ(model.isLog(), false) ;
    model.log() ;
    ASSERT_EQ(model.isLog(), true) ;
    model.exp() ;
    ASSERT_EQ(model.isLog(), false) ;

    // exp
    model.exp() ;
    ASSERT_EQ(model.isLog(), false) ;
    std::vector<double> exp_exp({exp(1.),
                                 exp(1.),
                                 exp(2.),
                                 exp(1.),
                                 exp(1.)}) ;

    // IPD
    auto& histograms_ipd = model.getHistogramsIPD() ;
    for(auto& histogram : histograms_ipd)
    {   size_t i=0 ;
        for(const auto& x : bh::indexed(histogram, bh::coverage::all))
        {   ASSERT_EQ(*x, exp_exp[i]) ;
            i++ ;
        }
    }

    // PWD
    auto& histograms_pwd = model.getHistogramsPWD() ;
    for(auto& histogram : histograms_pwd)
    {   size_t i=0 ;
        for(const auto& x : bh::indexed(histogram, bh::coverage::all))
        {   ASSERT_EQ(*x, exp_exp[i]) ;
            i++ ;
        }
    }
}


/*!
 * \brief Tests that add() method with KineticSignal works 
 * properly.
 */
TEST(RawKineticModelTest, add_kineticmodel)
{
    double pseudocounts(2.) ;

    RawKineticModelToy model1, model2, model3, model4 ;

    double inf = std::numeric_limits<double>::infinity() ;

    std::string seq("AAAAA") ;

    // add +1 count in [-inf,-2) bin at each position
    std::vector<double> val_1({-inf, -inf, -inf, -inf, -inf}) ;

    // add +1 count in [-2,-1) bin at each position
    std::vector<double> val_2({ -2.,  -2.,  -2.,  -2.,  -2.}) ;

    // add +1 count in [-1,0) bin at each position
    std::vector<double> val_3({ -1.,  -1.,  -1.,  -1.,  -1.}) ;

    // add +1 count in [0,1) bin at each position
    std::vector<double> val_4({  0.,   0.,   0.,   0.,   0.}) ;

    // add +1 count in [1,2) bin at each position
    std::vector<double> val_5({  1.,   1.,   1.,   1.,   1.}) ;

    // add +1 count in [2, inf) bin at each position
    std::vector<double> val_6({inf, inf, inf, inf, inf}) ;

    ngsai::KineticSignal k1(seq, seq, val_1, val_1, val_1, val_1) ;
    ngsai::KineticSignal k2(seq, seq, val_2, val_2, val_2, val_2) ;
    ngsai::KineticSignal k3(seq, seq, val_3, val_3, val_3, val_3) ;
    ngsai::KineticSignal k4(seq, seq, val_4, val_4, val_4, val_4) ;
    ngsai::KineticSignal k5(seq, seq, val_5, val_5, val_5, val_5) ;
    ngsai::KineticSignal k6(seq, seq, val_6, val_6, val_6, val_6) ;

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

    // add +1 count to each bin at each pos
    model1.add(k1) ;
    model1.add(k2) ;
    model1.add(k3) ;
    model1.add(k4) ;
    model1.add(k5) ;
    model1.add(k6) ;

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
    {   for(const auto& x : bh::indexed(histogram, bh::coverage::all))
        {   ASSERT_DOUBLE_EQ(*x, 2.*(pseudocounts + 1)) ; }
    }
    auto& histograms_pwd = model2.getHistogramsPWD() ;
    for(auto& histogram : histograms_pwd)
    {   for(const auto& x : bh::indexed(histogram, bh::coverage::all))
        {   ASSERT_DOUBLE_EQ(*x, 2.*(pseudocounts + 1)) ; }
    }
    ASSERT_EQ(model2.size(), 5) ;
    ASSERT_EQ(model2.getIsInit(), true) ;
    ASSERT_EQ(model2.getIsDensity(), false) ;
    ASSERT_EQ(model2.getIsLog(), false) ;
    // model3 is 3x model1
    histograms_ipd = model3.getHistogramsIPD() ;
    for(auto& histogram : histograms_ipd)
    {   for(const auto& x : bh::indexed(histogram, bh::coverage::all))
        {   ASSERT_DOUBLE_EQ(*x, 3.*(pseudocounts + 1)) ; }
    }
    histograms_pwd = model3.getHistogramsPWD() ;
    for(auto& histogram : histograms_pwd)
    {   for(const auto& x : bh::indexed(histogram, bh::coverage::all))
        {   ASSERT_DOUBLE_EQ(*x, 3.*(pseudocounts + 1)) ; }
    }
    ASSERT_EQ(model3.size(), 5) ;
    ASSERT_EQ(model3.getIsInit(), true) ;
    ASSERT_EQ(model3.getIsDensity(), false) ;
    ASSERT_EQ(model3.getIsLog(), false) ;

    
    // sum density models
    model1.density() ;

    RawKineticModelToy model5, model6 ;
    model5.setParameters(5, -2, 2, 4, 0.) ;
    model6.setParameters(5, -2, 2, 4, 0.) ;
    // model5 and model6 are not density but it is allowed to count with 
    // density, log with density, etc
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
    {   for(const auto& x : bh::indexed(histogram, bh::coverage::all))
        {   ASSERT_DOUBLE_EQ(*x, (2.*(pseudocounts + 1) / 18.)) ; }
    }
    histograms_pwd = model5.getHistogramsPWD() ;
    for(auto& histogram : histograms_pwd)
    {   for(const auto& x : bh::indexed(histogram, bh::coverage::all))
        {   ASSERT_DOUBLE_EQ(*x, (2.*(pseudocounts + 1) / 18.)) ; }
    }
    ASSERT_EQ(model2.size(), 5) ;
    ASSERT_EQ(model2.getIsInit(), true) ;
    ASSERT_EQ(model2.getIsDensity(), false) ;
    ASSERT_EQ(model2.getIsLog(), false) ;
    // model6 is 3x model1
    histograms_ipd = model6.getHistogramsIPD() ;
    for(auto& histogram : histograms_ipd)
    {   for(const auto& x : bh::indexed(histogram, bh::coverage::all))
        {   ASSERT_DOUBLE_EQ(*x,(3.*(pseudocounts + 1) / 18.)) ; }
    }
    histograms_pwd = model6.getHistogramsPWD() ;
    for(auto& histogram : histograms_pwd)
    {   for(const auto& x : bh::indexed(histogram, bh::coverage::all))
        {   ASSERT_DOUBLE_EQ(*x, (3.*(pseudocounts + 1) / 18.)) ; }
    }
    ASSERT_EQ(model3.size(), 5) ;
    ASSERT_EQ(model3.getIsInit(), true) ;
    ASSERT_EQ(model3.getIsDensity(), false) ;
    ASSERT_EQ(model3.getIsLog(), false) ;

    
    // sum log density models
    model1.log() ;

    RawKineticModelToy model7, model8 ;
    // model7 and model8 are not log but it is allowed to count with 
    // density, log with density, etc
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
    {   for(const auto& x : bh::indexed(histogram, bh::coverage::all))
        {   ASSERT_DOUBLE_EQ(*x, 2.*log((pseudocounts + 1) / 18.)) ; }
    }
    histograms_pwd = model7.getHistogramsPWD() ;
    for(auto& histogram : histograms_pwd)
    {   for(const auto& x : bh::indexed(histogram, bh::coverage::all))
        {   ASSERT_DOUBLE_EQ(*x, 2.*log((pseudocounts + 1) / 18.)) ; }
    }
    ASSERT_EQ(model2.size(), 5) ;
    ASSERT_EQ(model2.getIsInit(), true) ;
    ASSERT_EQ(model2.getIsDensity(), false) ;
    ASSERT_EQ(model2.getIsLog(), false) ;
    // model8 is 3x model1
    histograms_ipd = model8.getHistogramsIPD() ;
    for(auto& histogram : histograms_ipd)
    {   for(const auto& x : bh::indexed(histogram, bh::coverage::all))
        {   ASSERT_DOUBLE_EQ(*x, 3.*log((pseudocounts + 1) / 18.)) ; }
    }
    histograms_pwd = model8.getHistogramsPWD() ;
    for(auto& histogram : histograms_pwd)
    {   for(const auto& x : bh::indexed(histogram, bh::coverage::all))
        {   ASSERT_DOUBLE_EQ(*x, 3.*log((pseudocounts + 1) / 18.)) ; }
    }
    ASSERT_EQ(model3.size(), 5) ;
    ASSERT_EQ(model3.getIsInit(), true) ;
    ASSERT_EQ(model3.getIsDensity(), false) ;
    ASSERT_EQ(model3.getIsLog(), false) ;
}


/*!
 * \brief Tests that the logLikelihood() method works properly.
 */
TEST(RawKineticModelTest, logLikelihood)
{      
    std::string seq("AAA") ;

    std::vector<double> val_1({-1., -1., -1.}) ;
    std::vector<double> val_2({ 0.,  0.,  0.}) ;
    std::vector<double> val_3({ 1.,  1.,  1.}) ;
    std::vector<double> val_4({ 1.,  1.,  1.}) ;
    std::vector<double> val_5({ 2.,  2.,  2.}) ;
    std::vector<double> val_6({ 4.,  4.,  4.}) ;

    ngsai::KineticSignal k1(seq, seq, val_1, val_1, val_1, val_1) ;
    ngsai::KineticSignal k2(seq, seq, val_2, val_2, val_2, val_2) ;
    ngsai::KineticSignal k3(seq, seq, val_3, val_3, val_3, val_3) ;
    ngsai::KineticSignal k4(seq, seq, val_4, val_4, val_4, val_4) ;
    ngsai::KineticSignal k5(seq, seq, val_5, val_5, val_5, val_5) ;
    ngsai::KineticSignal k6(seq, seq, val_6, val_6, val_6, val_6) ;

    RawKineticModelToy model ;

    // not initialized
    EXPECT_THROW(model.logLikelihood(k3), std::runtime_error) ;

    // 3+2 bins (-inf,0), [0,1), [1,2), [2,3), [3,inf)
    model.setParameters(3, 0., 3., 3, 0.) ;

    model.add(k1) ;
    model.add(k2) ;
    model.add(k3) ;
    model.add(k4) ;
    model.add(k5) ;
    model.add(k6) ;


    std::string seq2("AAA") ;
    std::vector<double> val_s1( {1, 1, 1}) ;
    std::vector<double> val_s2( {-1, 0, 1}) ;
    std::vector<double> val_s3(  {0, 1, 2}) ;
    std::vector<double> val_s4(  {1, 2, 3}) ;
    std::vector<double> val_s5({-1,-1,-1}) ;
    std::vector<double> val_s6({ 0, 0, 0}) ;
    std::vector<double> val_s7({ 1, 1, 1}) ;
    std::vector<double> val_s8({ 2, 2, 2}) ;
    std::vector<double> val_s9({ 3, 3, 3}) ;

    // not log density
    ngsai::KineticSignal s1(seq, seq, val_s1, val_s1, val_s1, val_s1) ;
    EXPECT_THROW(model.logLikelihood(s1), std::runtime_error) ;

    model.density() ;
    
    // density but not log
    EXPECT_THROW(model.logLikelihood(s1), std::runtime_error) ;

    model.log() ;

    ngsai::KineticSignal s2(seq, seq, val_s2, val_s2, val_s2, val_s2) ;
    ngsai::KineticSignal s3(seq, seq, val_s3, val_s3, val_s3, val_s3) ;
    ngsai::KineticSignal s4(seq, seq, val_s4, val_s4, val_s4, val_s4) ;
    ngsai::KineticSignal s5(seq, seq, val_s5, val_s5, val_s5, val_s5) ;
    ngsai::KineticSignal s6(seq, seq, val_s6, val_s6, val_s6, val_s6) ;
    ngsai::KineticSignal s7(seq, seq, val_s7, val_s7, val_s7, val_s7) ;
    ngsai::KineticSignal s8(seq, seq, val_s8, val_s8, val_s8, val_s8) ;
    ngsai::KineticSignal s9(seq, seq, val_s9, val_s9, val_s9, val_s9) ;

    // must work          
    ASSERT_DOUBLE_EQ(model.logLikelihood(s2),
                     log(1./6.) + log(1./6.) + log(2./6.) + 
                     log(1./6.) + log(1./6.) + log(2./6.)) ;
    
    ASSERT_DOUBLE_EQ(model.logLikelihood(s3),
                     log(1./6.) + log(2./6.) + log(1./6.) +
                     log(1./6.) + log(2./6.) + log(1./6.)) ;
    
    ASSERT_DOUBLE_EQ(model.logLikelihood(s4),
                     log(2./6.) + log(1./6.) + log(1./6.) +
                     log(2./6.) + log(1./6.) + log(1./6.)) ;

    ASSERT_DOUBLE_EQ(model.logLikelihood(s5),
                     log(1./6.) + log(1./6.) + log(1./6.) +
                     log(1./6.) + log(1./6.) + log(1./6.)) ;

    ASSERT_DOUBLE_EQ(model.logLikelihood(s6),
                     log(1./6.) + log(1./6.) + log(1./6.) +
                     log(1./6.) + log(1./6.) + log(1./6.)) ;
    
    ASSERT_DOUBLE_EQ(model.logLikelihood(s7),
                     log(2./6.) + log(2./6.) + log(2./6.) +
                     log(2./6.) + log(2./6.) + log(2./6.)) ;

    ASSERT_DOUBLE_EQ(model.logLikelihood(s8),
                     log(1./6.) + log(1./6.) + log(1./6.) + 
                     log(1./6.) + log(1./6.) + log(1./6.)) ;
    
    ASSERT_DOUBLE_EQ(model.logLikelihood(s9),
                     log(1./6.) + log(1./6.) + log(1./6.) +
                     log(1./6.) + log(1./6.) + log(1./6.)) ;
    

    // must return nan because data have same size as model but missing data
    ngsai::KineticSignal s1_inc, s2_inc, s3_inc, s4_inc, s5_inc, s6_inc ;
    // misses seq fw
                                 s1_inc.setSequenceRv(seq2) ;
    s1_inc.setIPDFw(val_s2)    ; s1_inc.setIPDRv(val_s2) ;
    s1_inc.setPWDFw(val_s2)    ; s1_inc.setPWDRv(val_s2) ;
    // misses seq rv
    s2_inc.setSequenceFw(seq2) ;
    s2_inc.setIPDFw(val_s2)    ; s2_inc.setIPDRv(val_s2) ;
    s2_inc.setPWDFw(val_s2)    ; s2_inc.setPWDRv(val_s2) ;
    // misses IPD fw
    s3_inc.setSequenceFw(seq2) ; s3_inc.setSequenceRv(seq2) ;
                                 s3_inc.setIPDRv(val_s2)    ;
    s3_inc.setPWDFw(val_s2)    ; s3_inc.setPWDFw(val_s2) ;
    // misses IPD rv
    s4_inc.setSequenceFw(seq2) ; s4_inc.setSequenceRv(seq2) ;
    s4_inc.setIPDFw(val_s2)    ;
    s4_inc.setPWDFw(val_s2)    ; s4_inc.setPWDRv(val_s2) ;
    // misses PWD fw
    s5_inc.setSequenceFw(seq2) ; s5_inc.setSequenceRv(seq2) ;
    s5_inc.setIPDFw(val_s2)    ; s5_inc.setIPDRv(val_s2)    ;
                                 s5_inc.setPWDRv(val_s2) ;
    // misses PWD rv
    s6_inc.setSequenceFw(seq2) ; s6_inc.setSequenceRv(seq2) ;
    s6_inc.setIPDFw(val_s2)    ; s6_inc.setIPDRv(val_s2)    ;
    s6_inc.setPWDFw(val_s2)    ;
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
    std::vector<double> val_e1({1}) ;
    std::vector<double> val_e2({1,1}) ;
    std::vector<double> val_e4({1,1,1,1}) ;
    ngsai::KineticSignal s1_err(seq_e1,seq_e1,val_e1,val_e1,val_e1,val_e1) ;
    ngsai::KineticSignal s2_err(seq_e2,seq_e2,val_e2,val_e2,val_e2,val_e2) ;
    ngsai::KineticSignal s4_err(seq_e4,seq_e4,val_e4,val_e4,val_e4,val_e4) ;
    EXPECT_THROW(model.logLikelihood(s1_err), std::invalid_argument) ; // tpo short
    EXPECT_THROW(model.logLikelihood(s2_err), std::invalid_argument) ; // tpo short
    EXPECT_THROW(model.logLikelihood(s4_err), std::invalid_argument) ; // too long
}

