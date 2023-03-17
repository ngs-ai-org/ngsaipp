#include <ngsaipp/epigenetics/RawKineticModel.hpp>

#include <utility>
#include <stdexcept>                   // std::invalid_argument, std::runtime_error, std::out_of_range
#include <cmath>                       // log(), exp(), nan()
#include <sstream>                     // std::stringstream
#include <fstream>
#include <pbbam/CompositeBamReader.h>  // GenomicIntervalCompositeBamReader
#include <pbbam/BamRecord.h>      
#include <boost/histogram.hpp>
#include <boost/archive/text_iarchive.hpp>      // boost::archive::text_iarchive
#include <boost/archive/text_oarchive.hpp>      // boost::archive::text_oarchive
#include <boost/serialization/utility.hpp>      // stl data structures serialization

#include <ngsaipp/io/bed_io.hpp>                        // ngsai::BedReader, ngsai::BedRecord
#include <ngsaipp/genome/constants.hpp>                 // ngsai::genome::strand
#include <ngsaipp/epigenetics/CcsKineticExtractor.hpp>  // ngsai::CcsKineticExtractor
#include <ngsaipp/epigenetics/KineticSignal.hpp>

template<class T>
std::ostream&
ngsai::operator << (std::ostream& stream, const std::vector<T>& v)
{   for(const auto& x : v)
    {   stream << x << " " ; }
    return stream ;
}


namespace bh = boost::histogram ;


ngsai::RawKineticModel::RawKineticModel()
    : KineticModel(),
      m_histograms_ipd(0),
      m_histograms_pwd(0)
{ ; }


ngsai::RawKineticModel::RawKineticModel(const ngsai::RawKineticModel& other)
    : KineticModel(other),
      m_histograms_ipd(other.m_histograms_ipd),
      m_histograms_pwd(other.m_histograms_pwd)
{ ; }


ngsai::RawKineticModel::RawKineticModel(ngsai::RawKineticModel&& other)
    : KineticModel(other),
      m_histograms_ipd(std::move(other.m_histograms_ipd)),
      m_histograms_pwd(std::move(other.m_histograms_pwd))

{ ; }


ngsai::RawKineticModel::~RawKineticModel()
{ ; }



ngsai::RawKineticModel&
ngsai::RawKineticModel::operator = (const ngsai::RawKineticModel& other)
{   m_size           = other.m_size ;
    m_histograms_ipd = other.m_histograms_ipd ;
    m_histograms_pwd = other.m_histograms_pwd ;
    m_is_init        = other.m_is_init ;
    m_is_density     = other.m_is_density ;
    m_is_log         = other.m_is_log ;
    return *this ;
}


ngsai::RawKineticModel&
ngsai::RawKineticModel::operator = (ngsai::RawKineticModel&& other)
{   m_size           = std::move(other.m_size) ;
    m_histograms_ipd = std::move(other.m_histograms_ipd) ;
    m_histograms_pwd = std::move(other.m_histograms_pwd) ;
    m_is_init        = std::move(other.m_is_init) ;
    m_is_density     = std::move(other.m_is_density) ;
    m_is_log         = std::move(other.m_is_log) ;
    return *this ;
}


bool
ngsai::RawKineticModel::operator == (const ngsai::RawKineticModel& other) const
{   
    if(ngsai::KineticModel::operator!=(other))
    {   return false ; }

    if((m_histograms_ipd.size() != other.m_histograms_ipd.size()) or
       (m_histograms_pwd.size() != other.m_histograms_pwd.size()))
    {   return false ; }

    for(size_t i=0; i<m_histograms_ipd.size(); i++)
    {   if(m_histograms_ipd[i] != other.m_histograms_ipd[i])
        {   return false ; }
    }
    for(size_t i=0; i<m_histograms_pwd.size(); i++)
    {   if(m_histograms_pwd[i] != other.m_histograms_pwd[i])
        {   return false ; }
    }

    return true ;
}

bool
ngsai::RawKineticModel::operator != (const ngsai::RawKineticModel& other) const
{   return not ((*this) == other) ; }


void 
ngsai::RawKineticModel::setParameters(size_t size,
                                      double x_min,
                                      double x_max,
                                      size_t n_bins,
                                      double pseudo_counts)
{   if(size == 0)
    {   char msg[4096] ;
        sprintf(msg, 
                "RawKineticModel error : model size must be >0 (%zu)", 
                size) ;
        throw std::invalid_argument(msg) ;
    }
    else if(n_bins == 0)
    {   char msg[4096] ;
        sprintf(msg, 
                "RawKineticModel error : histograms must have >0 bins (%zu)", 
                n_bins) ;
        throw std::invalid_argument(msg) ;
    }
    if(not (x_min < x_max))
    {   char msg[4096] ;
        sprintf(msg, 
                "RawKineticModel error : histogram minimum (%lf) must be "
                "bigger than histogram maximum (%lf)", 
                x_min, 
                x_max) ;
        throw std::invalid_argument(msg) ;
    }
    m_size = size ;
    m_histograms_ipd = std::vector<hist_1d_double>(m_size) ;
    m_histograms_pwd = std::vector<hist_1d_double>(m_size) ;

    // histogram with uniform pseudo count 
    hist_1d_double pseudo_count_ipd = 
        bh::make_histogram(bh::axis::regular<>(n_bins,
                                               x_min, 
                                               x_max, 
                                               "ipd")) ;
    for (const auto& x : bh::indexed(pseudo_count_ipd, bh::coverage::all))
    {   (*x) += pseudo_counts ; }
    hist_1d_double pseudo_count_pwd = 
        bh::make_histogram(bh::axis::regular<>(n_bins,
                                               x_min, 
                                               x_max, 
                                               "pwd")) ;
    for (const auto& x : bh::indexed(pseudo_count_pwd, bh::coverage::all))
    {   (*x) += pseudo_counts ; }

    // add pseudo counts to models
    for(size_t i=0; i<m_size; i++)
    {   // IPDs
        m_histograms_ipd[i] = 
                    bh::make_histogram(bh::axis::regular<>(n_bins,
                                                           x_min,
                                                           x_max, 
                                                           "ipd")) ;
        m_histograms_ipd[i] += pseudo_count_ipd ;

        // PWDs
        m_histograms_pwd[i] = 
                    bh::make_histogram(bh::axis::regular<>(n_bins,
                                                           x_min,
                                                           x_max, 
                                                           "pwd")) ;
        m_histograms_pwd[i] += pseudo_count_pwd ;
    }
    m_is_init = true ;
}


std::vector<std::pair<double,double>>
ngsai::RawKineticModel::getBinBoundaries() const
{   if(not m_is_init)
    {   throw std::runtime_error("RawKineticModel error : cannot get the bin "
                                 "boundaries of a model that has not been "
                                 "initialized") ;
    }
    std::vector<std::pair<double,double>> bins ;
    for (const auto& x : bh::indexed(m_histograms_ipd[0], bh::coverage::all))
    {   bins.push_back(std::make_pair(x.bin().lower(),
                                      x.bin().upper())) ;
    }
    return bins ;
}


void 
ngsai::RawKineticModel::add(double ipd, double pwd, size_t position)
{   
    if(not m_is_init)
    {   throw std::runtime_error("RawKineticModel error : cannot add data to "
                                 "a model that has not been initialized") ;
    }
    if(position > m_size - 1)
    {   char msg[4096]  ;
        sprintf(msg,
               "RawKineticModel error : position (%zu) is out of range for "
               "model size (%zu)",
                position,
                m_size) ;
        throw std::out_of_range(msg) ;
    }
    m_histograms_ipd[position](ipd) ;
    m_histograms_pwd[position](pwd) ;
}


void 
ngsai::RawKineticModel::add(const ngsai::KineticSignal& kinetics)
{   
    if(not m_is_init)
    {   throw std::runtime_error("RawKineticModel error : cannot add data to "
                                 "a model that has not been initialized") ;
    }
    if(kinetics.size() != m_size)
    {   char msg[4096]  ;
        sprintf(msg,
               "RawKineticModel error : number of IPDs/PWDs (%zu) do not " 
               "match model size (%zu)",
               kinetics.size(),
               m_size) ;
        throw std::invalid_argument(msg) ;
    }
    if(not kinetics.isComplete())
    {   throw std::invalid_argument("RawKineticModel error : KineticSignal "
                                   "does not have complete data") ;
    }
    
    // compute average (fw and rv) of IPD / PWD
    std::vector<double> ipds = kinetics.getIPDMean() ;
    std::vector<double> pwds = kinetics.getPWDMean() ;

    for(size_t i=0; i<kinetics.size(); i++)
    {   m_histograms_ipd[i](ipds[i]) ;
        m_histograms_pwd[i](pwds[i]) ; 
    }
}


void
ngsai::RawKineticModel::add(const KineticModel& model)
{   model.addTo(*this) ; }


ngsai::RawKineticModel*
ngsai::RawKineticModel::copy() const
{
    return new ngsai::RawKineticModel(*this) ;
}


void
ngsai::RawKineticModel::density()
{   if(not m_is_init)
    {   throw std::invalid_argument("RawKineticModel error : cannot compute "
                                    "densities from a model that has not been "
                                    "initialized") ;
    }
    for(size_t i=0; i<m_size; i++)
    {   m_histograms_ipd[i] /= 
                bh::algorithm::sum(m_histograms_ipd[i]) ;
        m_histograms_pwd[i] /= 
                bh::algorithm::sum(m_histograms_pwd[i]) ;
    }
    m_is_density = true ;
}


void 
ngsai::RawKineticModel::log()
{   if(not m_is_init)
    {   throw std::invalid_argument("RawKineticModel error : cannot log "
                                    "transform a model that has not been "
                                    "initialized") ;
    }
    // transform to log
    for(size_t i=0; i<m_size; i++)
    {   for(const auto& x : 
            bh::indexed(m_histograms_ipd[i], bh::coverage::all))
        {   *x = std::log(*x) ; }
        for(const auto& x : 
            bh::indexed(m_histograms_pwd[i], bh::coverage::all))
        {   *x = std::log(*x) ; }
    }
    m_is_log = true ;
}


void 
ngsai::RawKineticModel::exp()
{   if(not m_is_init)
    {   throw std::invalid_argument("RawKineticModel error : cannot "
                                    "exponantialize a model that has not been "
                                    "initialized") ;
    }
    // transform to exp
    for(size_t i=0; i<m_size; i++)
    {   for(const auto& x : 
            bh::indexed(m_histograms_ipd[i], bh::coverage::all))
        {   *x = std::exp(*x) ; }
        for(const auto& x : 
            bh::indexed(m_histograms_pwd[i], bh::coverage::all))
        {   *x = std::exp(*x) ; }
    }
    m_is_log = false ;
}


double
ngsai::RawKineticModel::logLikelihood(const ngsai::KineticSignal& kinetics) const
{   
    // model side errors
    if(not m_is_init)
    {   throw std::runtime_error("RawKineticModel error : cannot compute a "
                                 "loglikelihood from a model that has not "
                                 "been initialized") ;
    }
    if(not m_is_density)
    {   throw std::runtime_error("RawKineticModel error :cannot compute a "
                                 "loglikelihood from a model that does not "
                                 "contain signal densities") ;
    }
    if(not m_is_log)
    {   throw std::runtime_error("RawKineticModel error : cannot compute a "
                                 "loglikelihood from a model that does not "
                                 "contain signal log densities") ;

    }

    // cannot compute from something incomplete
    if(not kinetics.isComplete())
    {   return std::nan("") ; }

    // kinetic side error
    if(kinetics.size() != m_size)
    {   char msg[4096]  ;
        sprintf(msg,
               "RawKineticModel error : number of IPDs / PWDs (%zu) do not "
               "match model size (%zu)",
               kinetics.size(),
               m_size) ;
        throw std::invalid_argument(msg) ;
    }

    // compute average (fw and rv) of IPD / PWD
    std::vector<double> ipds = kinetics.getIPDMean() ;
    std::vector<double> pwds = kinetics.getPWDMean() ;    
    
    // log likelihood
    double ll = 0. ;
    for(size_t i=0; i<kinetics.size(); i++)
    {   size_t ipd_idx = m_histograms_ipd[i].axis(0).index(ipds[i]) ;
        size_t pwd_idx = m_histograms_pwd[i].axis(0).index(pwds[i]) ;
        ll += m_histograms_ipd[i].at(ipd_idx) +
              m_histograms_pwd[i].at(pwd_idx) ; 
    }

    return ll ;
}


std::string
ngsai::RawKineticModel::toString() const
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

    // IPDs
    // positions
    for(size_t i=0; i<this->size(); i++)
    {   stream << "ipd_" << i+1 << '\t' ;
        size_t j = 0 ;
        for(const auto& x : 
            bh::indexed(this->m_histograms_ipd[i], bh::coverage::all))
        {   stream << *x ; 
            if(j < (bins.size() - 1))
            {   stream << '\t' ; }
            j++ ;
        }
        stream << std::endl ;

    }

    // PWDs
    // positions
    for(size_t i=0; i<this->size(); i++)
    {   stream << "pwd_" << i+1 << '\t' ;
        size_t j = 0 ;
        for(const auto& x : 
            bh::indexed(this->m_histograms_pwd[i], bh::coverage::all))
        {   stream << *x ; 
            if(j < (bins.size() - 1))
            {   stream << '\t' ; }
            j++ ;
        }
        stream << std::endl ;

    }
    return stream.str() ;
}


void
ngsai::RawKineticModel::save(const std::string& path) const
{   
    std::ofstream f_out(path) ;
    boost::archive::text_oarchive arch_out(f_out) ;
    arch_out << *this ;
    f_out.close() ;
}


void
ngsai::RawKineticModel::load(const std::string& path)
{   std::ifstream f_in(path) ;
    boost::archive::text_iarchive arch_in(f_in) ;
    arch_in >> *this ;
    f_in.close() ;
}


std::ostream& 
ngsai::operator << (
    std::ostream& stream, 
    const ngsai::RawKineticModel& model)
{   stream << model.toString() ;
    return stream ;
}


void 
ngsai::RawKineticModel::addTo(
            ngsai::RawKineticModel& model) const
{   if(not m_is_init)
    {   throw std::invalid_argument(
            "RawKineticModel error : this model has not "
            "been initialized") ;
    }
     if(not model.isInit())
    {   throw std::invalid_argument(
            "RawKineticModel error : given model has not "
            "been initialized") ;
    }
    if(this->size() != model.size())
    {   char msg[4096] ; 
        sprintf(msg,
                "RawKineticModel error : cannot sum models "
                "having different sizes (%zu, %zu)",
                this->size(),
                model.size()) ;
        throw(std::invalid_argument(msg)) ;
    }

    try
    {   for(size_t i=0; i<model.size(); i++)
        {   model.m_histograms_ipd[i] += 
                            m_histograms_ipd[i] ;
            model.m_histograms_pwd[i] += 
                            m_histograms_pwd[i] ;
        }
    }
    catch(const std::exception& e)
    {   throw std::runtime_error(e.what()) ; }
}


void
ngsai::RawKineticModel::addTo(
            ngsai::NormalizedKineticModel&) const
{   throw std::invalid_argument(
        "RawKineticModel error : cannot add a "
        "RawKineticModel to a NormalizedKineticModel") ;
}


void
ngsai::RawKineticModel::addTo(
            ngsai::PairWiseKineticModel&) const
{   throw std::invalid_argument(
        "RawKineticModel error : cannot add a "
        "RawKineticModel to a PairWiseKineticModel") ;
}


void
ngsai::RawKineticModel::addTo(
            ngsai::PairWiseNormalizedKineticModel&) const
{   throw std::invalid_argument(
        "RawKineticModel error : cannot add a "
        "RawKineticModel to a "
        "PairWiseNormalizedKineticModel") ;
}


void
ngsai::RawKineticModel::addTo(
            ngsai::DiPositionKineticModel&) const
{   throw std::invalid_argument(
        "RawKineticModel error : cannot add a "
        "RawKineticModel to a DiPositionKineticModel") ;
}


void
ngsai::RawKineticModel::addTo(
            ngsai::DiPositionNormalizedKineticModel&) const
{   throw std::invalid_argument(
        "RawKineticModel error : cannot add a "
        "RawKineticModel to a "
        "DiPositionNormalizedKineticModel") ;
}