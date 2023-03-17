#include <ngsaipp/epigenetics/PairWiseKineticModel.hpp>

#include <utility>
#include <stdexcept>                   // std::invalid_argument, std::runtime_error, std::out_of_range
#include <cmath>                       // log(), exp(), nan()
#include <sstream>                     // std::stringstream
#include <fstream>
#include <set>                         // std::set
#include <map>                         // std::map
#include <pbbam/CompositeBamReader.h>  // GenomicIntervalCompositeBamReader
#include <pbbam/BamRecord.h>      
#include <boost/histogram.hpp>
#include <boost/archive/text_iarchive.hpp>      // boost::archive::text_iarchive
#include <boost/archive/text_oarchive.hpp>      // boost::archive::text_oarchive
#include <boost/serialization/utility.hpp>      // stl data structures serialization

#include <ngsaipp/io/bed_io.hpp>                         // ngsai::BedReader, ngsai::BedRecord
#include <ngsaipp/genome/constants.hpp>                  // ngsai::genome::strand
#include <ngsaipp/epigenetics/CcsKineticExtractor.hpp>   // ngsai::CcsKineticExtractor
#include <ngsaipp/epigenetics/KineticSignal.hpp>         // ngsai::KineticSignal


namespace bh = boost::histogram ;



ngsai::PairWiseKineticModel::PairWiseKineticModel()
    : KineticModel(),
      m_pair_nb(0),
      m_pair_index(0),
      m_histograms_ipd(0),
      m_histograms_pwd(0)
{ ; }


ngsai::PairWiseKineticModel::PairWiseKineticModel(
                    const ngsai::PairWiseKineticModel& other)
    : KineticModel(other),
      m_pair_nb(other.m_pair_nb),
      m_pair_index(other.m_pair_index),
      m_histograms_ipd(other.m_histograms_ipd),
      m_histograms_pwd(other.m_histograms_pwd)
{ ; }


ngsai::PairWiseKineticModel::PairWiseKineticModel(
                    ngsai::PairWiseKineticModel&& other)
    : KineticModel(other),
      m_pair_nb(std::move(other.m_pair_nb)),
      m_pair_index(std::move(other.m_pair_index)),
      m_histograms_ipd(std::move(other.m_histograms_ipd)),
      m_histograms_pwd(std::move(other.m_histograms_pwd))

{ ; }


ngsai::PairWiseKineticModel::~PairWiseKineticModel()
{ ; }


ngsai::PairWiseKineticModel&
ngsai::PairWiseKineticModel::operator = (const ngsai::PairWiseKineticModel& other)
{   m_size           = other.m_size ;
    m_pair_nb        = other.m_pair_nb ;
    m_pair_index     = other.m_pair_index ;
    m_histograms_ipd = other.m_histograms_ipd ;
    m_histograms_pwd = other.m_histograms_pwd ;
    m_is_init        = other.m_is_init ;
    m_is_density     = other.m_is_density ;
    m_is_log         = other.m_is_log ;
    return *this ;
}


ngsai::PairWiseKineticModel&
ngsai::PairWiseKineticModel::operator = (ngsai::PairWiseKineticModel&& other)
{   m_size           = std::move(other.m_size) ;
    m_pair_nb        = std::move(other.m_pair_nb) ;
    m_pair_index     = std::move(other.m_pair_index) ;
    m_histograms_ipd = std::move(other.m_histograms_ipd) ;
    m_histograms_pwd = std::move(other.m_histograms_pwd) ;
    m_is_init        = std::move(other.m_is_init) ;
    m_is_density     = std::move(other.m_is_density) ;
    m_is_log         = std::move(other.m_is_log) ;
    return *this ;
}


bool
ngsai::PairWiseKineticModel::operator == (
                                const ngsai::PairWiseKineticModel& other) const
{   if(ngsai::KineticModel::operator!=(other))
    {   return false ; }

    if(m_pair_nb != other.m_pair_nb)
    {   return false ; }

    for(size_t i=0; i<m_pair_nb; i++)
    {   const auto& p1 = m_pair_index[i] ;
        const auto& p2 = other.m_pair_index[i] ;
        if(p1 != p2)
        {   return false ; }
    }

    for(const auto& p : m_pair_index)
    {   auto pos1 = p.first ;
        auto pos2 = p.second ;
        if(m_histograms_ipd[pos1][pos2] != other.m_histograms_ipd[pos1][pos2])
        {   return false ;}
        if(m_histograms_pwd[pos1][pos2] != other.m_histograms_pwd[pos1][pos2])
        {   return false ;}
    }

    return true ;
}

bool
ngsai::PairWiseKineticModel::operator != (
                                const ngsai::PairWiseKineticModel& other) const
{   return not ((*this) == other) ; }


void 
ngsai::PairWiseKineticModel::setParameters(size_t size,
                                             double x_min,
                                             double x_max,
                                             size_t n_bins,
                                             double pseudo_counts)
{   if(size < 2)
    {   char msg[4096] ;
        sprintf(msg, 
                "PairWiseKineticModel error : model size must be >1 (%zu)", 
                size) ;
        throw std::invalid_argument(msg) ;
    }
    else if(n_bins == 0)
    {   char msg[4096] ;
        sprintf(msg, 
                "PairWiseKineticModel error : histograms must have >0 bins "
                "(%zu)", 
                n_bins) ;
        throw std::invalid_argument(msg) ;
    }
    if(not (x_min < x_max))
    {   char msg[4096] ;
        sprintf(msg, 
                "PairWiseKineticModel error : histogram minimum (%lf) must "
                "be bigger than histogram maximum (%lf)", 
                x_min, 
                x_max) ;
        throw std::invalid_argument(msg) ;
    }
    m_size = size ;
    // sets m_pair_index and m_pair_nb
    this->computeMeaningfulPairs() ;
    m_histograms_ipd = std::vector<std::vector<hist_2d_double>>(
                            m_size, std::vector<hist_2d_double>(m_size)) ; 
    m_histograms_pwd = std::vector<std::vector<hist_2d_double>>(
                            m_size, std::vector<hist_2d_double>(m_size)) ; 
   
    // histogram with uniform pseudo count 
    hist_2d_double pseudo_count = 
        bh::make_histogram(bh::axis::regular<>(n_bins,
                                               x_min, 
                                               x_max, 
                                               "pos_1"),
                           bh::axis::regular<>(n_bins,
                                               x_min, 
                                               x_max, 
                                               "pos_2")) ;
    for (const auto& x : bh::indexed(pseudo_count, bh::coverage::all))
    {   (*x) += pseudo_counts ; }

    // add pseudo counts to models
    for(const auto& pair_index : m_pair_index)
    {   const size_t& pos1 = pair_index.first ;
        const size_t& pos2 = pair_index.second ;
        // IPDs
        m_histograms_ipd[pos1][pos2] = 
                    bh::make_histogram(bh::axis::regular<>(n_bins,
                                                           x_min,
                                                           x_max, 
                                                           "pos_1"),
                                       bh::axis::regular<>(n_bins,
                                                           x_min,
                                                           x_max, 
                                                           "pos_2")) ;
        m_histograms_ipd[pos1][pos2] += pseudo_count ;

        // PWDs
        m_histograms_pwd[pos1][pos2] = 
                    bh::make_histogram(bh::axis::regular<>(n_bins,
                                                           x_min,
                                                           x_max, 
                                                           "pos_1"),
                                       bh::axis::regular<>(n_bins,
                                                           x_min,
                                                           x_max, 
                                                           "pos_2")) ;
        m_histograms_pwd[pos1][pos2] += pseudo_count ;
    }
    m_is_init = true ;
}


std::vector<std::pair<double,double>>
ngsai::PairWiseKineticModel::getBinBoundaries() const
{   
    if(not m_is_init)
    {   throw std::runtime_error("PairWiseKineticModel error : cannot get "
                                 "the bin boundaries of a model that has not "
                                 "been initialized") ;
    }

    // the 2D histogram is squared -> same bins on 1st and 2nd axes
    // use a set because the enumerator list all bins one after the other
    // -> on a given axis, the bin boundaries are listed several times  
    // (and this happens for both axes)
    std::set<std::pair<double,double>> bins ;

    // only consider one of the histograms since they are all the same
    const size_t& pos1 = m_pair_index[0].first ;
    const size_t& pos2 = m_pair_index[0].second ;
    for(const auto& x : bh::indexed(m_histograms_ipd[pos1][pos2],
                                    bh::coverage::all))
    {    bins.emplace(std::make_pair(x.bin(0).lower(),
                                     x.bin(0).upper())) ;
    }
    return std::vector<std::pair<double,double>>(bins.begin(), bins.end()) ;
}


void 
ngsai::PairWiseKineticModel::add(const ngsai::KineticSignal& kinetics)
{   
    if(not m_is_init)
    {   throw std::runtime_error("PairWiseKineticModel error : cannot add "
                                 "data to a model that has not been "
                                 "initialized") ;
    }
    if(kinetics.size() != m_size)
    {   char msg[4096]  ;
        sprintf(msg,
               "PairWiseKineticModel error : number of IPDs / PWDs (%zu) does"
               "not match model size (%zu)",
               kinetics.size(),
               m_size) ;
        throw std::invalid_argument(msg) ;
    }
    if(not kinetics.isComplete())
    {   throw std::invalid_argument("PairWiseKineticModel error : "
                                    "KineticSignal does not have "
                                    "complete data") ;
    }
    
    // compute average (fw and rv) of IPD / PWD
    std::vector<double> ipds = kinetics.getIPDMean() ;
    std::vector<double> pwds = kinetics.getPWDMean() ;
    
    for(const auto& pair_index : m_pair_index)
    {   const size_t& pos1 = pair_index.first ;
        const size_t& pos2 = pair_index.second ;
        m_histograms_ipd[pos1][pos2](ipds[pos1], ipds[pos2]) ;
        m_histograms_pwd[pos1][pos2](pwds[pos1], pwds[pos2]) ;
    }
}


void
ngsai::PairWiseKineticModel::add(
                            const KineticModel& model)
{   model.addTo(*this) ; }


ngsai::PairWiseKineticModel*
ngsai::PairWiseKineticModel::copy() const
{
    return new ngsai::PairWiseKineticModel(*this) ;
}


void
ngsai::PairWiseKineticModel::density()
{   if(not m_is_init)
    {   throw std::invalid_argument("PairWiseKineticModel error : cannot "
                                    "compute densities from a model that has "
                                    "not been initialized") ;
    }
    for(const auto& pair_index : m_pair_index)
    {   const size_t& pos1 = pair_index.first ;
        const size_t& pos2 = pair_index.second ;
        m_histograms_ipd[pos1][pos2] /= 
                bh::algorithm::sum(m_histograms_ipd[pos1][pos2]) ;
        m_histograms_pwd[pos1][pos2] /= 
                bh::algorithm::sum(m_histograms_pwd[pos1][pos2]) ;
    }
    m_is_density = true ;
}


void 
ngsai::PairWiseKineticModel::log()
{   if(not m_is_init)
    {   throw std::invalid_argument("PairWiseKineticModel error : cannot "
                                    "log transform a model that has not been "
                                    "initialized") ;
    }
    // transform to log
    for(const auto& pair_index : m_pair_index)
    {   const size_t& pos1 = pair_index.first ;
        const size_t& pos2 = pair_index.second ;
        for(const auto& x : 
            bh::indexed(m_histograms_ipd[pos1][pos2], bh::coverage::all))
        {   *x = std::log(*x) ; }
        for(const auto& x : 
            bh::indexed(m_histograms_pwd[pos1][pos2], bh::coverage::all))
        {   *x = std::log(*x) ; }
    }
    m_is_log = true ;
}


void 
ngsai::PairWiseKineticModel::exp()
{   if(not m_is_init)
    {   throw std::invalid_argument("PairWiseKineticModel error : cannot "
                                    "exponantialize a model that has not been "
                                    "initialized") ;
    }
    // transform to exp
    for(const auto& pair_index : m_pair_index)
    {   const size_t& pos1 = pair_index.first ;
        const size_t& pos2 = pair_index.second ;
        for(const auto& x : 
            bh::indexed(m_histograms_ipd[pos1][pos2], bh::coverage::all))
        {   *x = std::exp(*x) ; }
        for(const auto& x : 
            bh::indexed(m_histograms_pwd[pos1][pos2], bh::coverage::all))
        {   *x = std::exp(*x) ; }
    }
    m_is_log = false ;
}


double
ngsai::PairWiseKineticModel::logLikelihood(
                                    const ngsai::KineticSignal& kinetics) const
{   // model side errors
    if(not m_is_init)
    {   throw std::runtime_error("PairWiseKineticModel error : cannot "
                                 "compute a loglikelihood from a model that "
                                 "has not been initialized") ;
    }
    if(not m_is_density)
    {   throw std::runtime_error("PairWiseKineticModel error :cannot "
                                 "compute a loglikelihood from a model that "
                                 "does not contain signal densities") ;
    }
    if(not m_is_log)
    {   throw std::runtime_error("PairWiseKineticModel error : cannot "
                                 "compute a loglikelihood from a model that "
                                 "does not contain signal log densities") ;

    }

    // cannot compute from something incomplete
    if(not kinetics.isComplete())
    {   return std::nan("") ; }
    
    // kinetic side errors
    if(kinetics.size() != m_size)
    {   char msg[4096]  ;
        sprintf(msg,
               "PairWiseKineticModel error : number of IPDs / PWDs (%zu) does"
               "not match model size (%zu)",
               kinetics.size(),
               m_size) ;
        throw std::invalid_argument(msg) ;
    }
    
   
    // compute average (fw and rv) of IPD / PWD
    std::vector<double> ipds = kinetics.getIPDMean() ;
    std::vector<double> pwds = kinetics.getPWDMean() ;

    // log likelihood
    double ll = 0. ;
    for(const auto& pair_index : m_pair_index)
    {   const size_t& pos1 = pair_index.first ;
        const size_t& pos2 = pair_index.second ;
        size_t ipds_idx_1 = m_histograms_ipd[pos1][pos2].axis(0).index(ipds[pos1]) ;
        size_t ipds_idx_2 = m_histograms_ipd[pos1][pos2].axis(1).index(ipds[pos2]) ;
        size_t pwds_idx_1 = m_histograms_pwd[pos1][pos2].axis(0).index(pwds[pos1]) ;
        size_t pwds_idx_2 = m_histograms_pwd[pos1][pos2].axis(1).index(pwds[pos2]) ;

        // std::cerr <<   "IDP " << ipds[pos1] << " " << ipds[pos2] 
        //           << "\tPWD " << pwds[pos1] << " " << pwds[pos2] << std::endl ; 
        // std::cerr << ll << " + " << m_histograms_ipd[pos1][pos2].at(ipds_idx_1, ipds_idx_2)
        //                 << " + " << m_histograms_pwd[pos1][pos2].at(pwds_idx_1, pwds_idx_2) << std::endl ;
        // std::cerr << "------------------------------------------------------" << std::endl ;

        ll += m_histograms_ipd[pos1][pos2].at(ipds_idx_1, ipds_idx_2) +
              m_histograms_pwd[pos1][pos2].at(pwds_idx_1, pwds_idx_2) ;
    }
    // std::cerr << "======================================================" << std::endl ;
    return ll ;
}


void
ngsai::PairWiseKineticModel::computeMeaningfulPairs()
{   m_pair_nb = ((m_size*m_size) - m_size) / 2 ;
    // m_pair_index.resize(m_pair_nb) ;
    // m_pair_index.empty() ;
    m_pair_index.clear() ;
    for(size_t i=0; i<m_size-1; i++)
    {   for(size_t j=i+1; j<m_size; j++)
        {   m_pair_index.emplace_back(std::make_pair(i,j)) ; }
    }
}


std::string
ngsai::PairWiseKineticModel::toString() const
{   
    if(not this->isInit())
    {   return std::string() ; }

    std::stringstream stream ;

    std::vector<std::pair<double,double>> bins = this->getBinBoundaries() ;

    // lower bin bound
    stream << "bin_lower" << '\t' ;
    for(size_t i=0; i<bins.size(); i++)
    {   stream << bins[i].first  ;
        if(i < (bins.size() - 1))
        {   stream << '\t' ; }
    }
    stream << std::endl ;

    // upper bin bound
    stream << "bin_upper" << '\t' ;
    for(size_t i=0; i<bins.size(); i++)
    {   stream << bins[i].second  ;
        if(i < (bins.size() - 1))
        {   stream << '\t' ; }
    }
    stream << std::endl ;

    // print each IPD histogram
    for(const auto& pair_index : m_pair_index)
    {   const size_t& pos1 = pair_index.first ;
        const size_t& pos2 = pair_index.second ;

        // convert to matrix
        auto matrix = ngsai::hist2d_to_matrix(m_histograms_ipd[pos1][pos2]) ;
        // print it with row names
        for(const auto& row : matrix)
        {   stream << "ipd_" << pos1+1 << "_" << pos2+1 << '\t' ;
            size_t j = 0 ;
            for(const auto& x : row)
            {   stream << x ;
                if(j < (bins.size() - 1))
                {   stream << '\t' ; }
                j++ ;
            }
            stream << std::endl ;
        }
    }

    // print each PWD histogram
    for(const auto& pair_index : m_pair_index)
    {   const size_t& pos1 = pair_index.first ;
        const size_t& pos2 = pair_index.second ;

        // convert to matrix
        auto matrix = ngsai::hist2d_to_matrix(m_histograms_pwd[pos1][pos2]) ;
        // print it with row names
        for(const auto& row : matrix)
        {   stream << "pwd_" << pos1+1 << "_" << pos2+1  << '\t' ;
            size_t j = 0 ;
            for(const auto& x : row)
            {   stream << x ;
                if(j < (bins.size() - 1))
                {   stream << '\t' ; }
                j++ ;
            }
            stream << std::endl ;
        }
    }
    return stream.str() ;
}


void
ngsai::PairWiseKineticModel::save(const std::string& path) const
{   
    std::ofstream f_out(path) ;
    boost::archive::text_oarchive arch_out(f_out) ;
    arch_out << *this ;
    f_out.close() ;
}


void
ngsai::PairWiseKineticModel::load(const std::string& path)
{   std::ifstream f_in(path) ;
    boost::archive::text_iarchive arch_in(f_in) ;
    arch_in >> *this ;
    f_in.close() ;
}


std::vector<std::vector<double>>
ngsai::hist2d_to_matrix(const hist_2d_double& h)
{   
    // get the bin boundaries and assign each bin with an index
    // cannot user x.bin().index(0/1) because for -inf/inf bins the index 
    // is std::numeric_limits<size_t>::max()...
    // need to construct it -> map use dichotomic search/insert
    // this means the container is always sorted \o/ (a bit slower than 
    // hasing using unordered_map but useful)
    size_t index_1 = 0 ;
    size_t index_2 = 0 ;
    std::map<std::pair<double,double>, size_t> bins_1 ;
    std::map<std::pair<double,double>, size_t> bins_2 ;
    for(const auto& x : bh::indexed(h,  bh::coverage::all))
    {      
        std::pair<double,double> p1(x.bin(0).lower(), x.bin(0).upper()) ;
        std::pair<double,double> p2(x.bin(1).lower(), x.bin(1).upper()) ;

        if(bins_1.find(p1) == bins_1.end())
        {   bins_1[p1] = index_1 ;
            index_1++ ;
        }

        if(bins_2.find(p2) == bins_2.end())
        {   bins_2[p2] = index_2 ;
            index_2++ ;
        }
    }

    // transform into matrix
    std::vector<std::vector<double>> m(bins_1.size(), 
                                       std::vector<double>(bins_2.size(), 
                                                           -1.)) ;
    for(const auto& x : bh::indexed(h,  bh::coverage::all))
    {      
        std::pair<double,double> p1(x.bin(0).lower(), x.bin(0).upper()) ;
        std::pair<double,double> p2(x.bin(1).lower(), x.bin(1).upper()) ;
        size_t i = bins_1.find(p1)->second ;
        size_t j = bins_2.find(p2)->second ;
        m[i][j] = *x ;
    }

    return  m ;
}


void 
ngsai::PairWiseKineticModel::addTo(
            ngsai::RawKineticModel&) const
{
   throw std::invalid_argument(
        "PairWiseKineticModel error : cannot add a "
        "PairWiseKineticModel to a RawKineticModel") ;
}


void
ngsai::PairWiseKineticModel::addTo(
            ngsai::NormalizedKineticModel&) const
{   throw std::invalid_argument(
        "PairWiseKineticModel error : cannot add a "
        "PairWiseKineticModel to a "
        "NormalizedKineticModel") ;
}


void 
ngsai::PairWiseKineticModel::addTo(
                ngsai::PairWiseKineticModel& model) const
{
    if(not m_is_init)
    {   throw std::invalid_argument(
            "PairWiseKineticModel error : this model has "
            "not been initialized") ;
    }
     if(not model.isInit())
    {   throw std::invalid_argument(
        "PairWiseKineticModel error : given model has "
        "not been initialized") ;
    }
    if(this->size() != model.size())
    {   char msg[4096] ; 
        sprintf(msg,
                "PairWiseKineticModel error : cannot "
                "add models having different sizes "
                "(%zu, %zu)",
                this->size(),
                model.size()) ;
        throw(std::invalid_argument(msg)) ;
    }

    try
    {   // for(size_t i=0; i<this->size(); i++)
        // {   for(size_t j=0; j<m_histograms_ipd[i].size(); j++)
        //     {   m_histograms_ipd[i][j] += other.m_histograms_ipd[i][j] ;
        //         m_histograms_pwd[i][j] += other.m_histograms_pwd[i][j] ;
        //     }
        // }
        for(const auto& [i,j] : model.m_pair_index)
        {   model.m_histograms_ipd[i][j] += 
                            m_histograms_ipd[i][j] ;
            model.m_histograms_pwd[i][j] +=
                            m_histograms_pwd[i][j] ;
        }
    }
    catch(const std::exception& e)
    {   throw std::runtime_error(e.what()) ; }
}


void
ngsai::PairWiseKineticModel::addTo(
            ngsai::PairWiseNormalizedKineticModel&) const
{   throw std::invalid_argument(
        "PairWiseKineticModel error : cannot add a "
        "PairWiseKineticModel to a "
        "PairWiseNormalizedKineticModel") ;
}


void
ngsai::PairWiseKineticModel::addTo(
            ngsai::DiPositionKineticModel&) const
{   throw std::invalid_argument(
        "PairWiseKineticModel error : cannot add a "
        "PairWiseKineticModel to a "
        "DiPositionKineticModel") ;
}


void
ngsai::PairWiseKineticModel::addTo(
            ngsai::DiPositionNormalizedKineticModel&) const
{   throw std::invalid_argument(
        "PairWiseKineticModel error : cannot add a "
        "PairWiseKineticModel to a "
        "DiPositionNormalizedKineticModel") ;
}


