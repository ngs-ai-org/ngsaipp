#include <ngsaipp/epigenetics/NormalizedKineticModel.hpp>

#include <utility>
#include <vector>
#include <fstream>
#include <stdexcept>
#include <boost/archive/text_iarchive.hpp>      // boost::archive::text_iarchive
#include <boost/archive/text_oarchive.hpp>      // boost::archive::text_oarchive
#include <boost/serialization/utility.hpp>      // stl data structures serialization

#include  <ngsaipp/epigenetics/KineticSignal.hpp>
#include  <ngsaipp/epigenetics/model_utility.hpp>   // ngsai::normalize_kinetics()


ngsai::NormalizedKineticModel::NormalizedKineticModel()
    : RawKineticModel(),
      m_has_kmermap(false),
      m_bckg_model()
{ ; }

ngsai::NormalizedKineticModel::NormalizedKineticModel(
                                    const ngsai::KmerMap& bckg_model)
    : RawKineticModel(),
      m_has_kmermap(true),
      m_bckg_model(bckg_model)
{ ; }


ngsai::NormalizedKineticModel::NormalizedKineticModel(
                                    const ngsai::NormalizedKineticModel& other)
    : RawKineticModel(other),
      m_has_kmermap(other.m_has_kmermap),
      m_bckg_model(other.m_bckg_model)
{ ; }


ngsai::NormalizedKineticModel::NormalizedKineticModel(
                                    ngsai::NormalizedKineticModel&& other)
    : RawKineticModel(other),
      m_has_kmermap(std::move(other.m_has_kmermap)),
      m_bckg_model(std::move(other.m_bckg_model))
{ ; }


ngsai::NormalizedKineticModel::~NormalizedKineticModel()
{ ; }


ngsai::NormalizedKineticModel&
ngsai::NormalizedKineticModel::operator = (
                                    const ngsai::NormalizedKineticModel& other)
{   ngsai::RawKineticModel::operator=(other) ;
    m_bckg_model = other.m_bckg_model ;
    return *this ;
}


ngsai::NormalizedKineticModel&
ngsai::NormalizedKineticModel::operator = (
                                    ngsai::NormalizedKineticModel&& other)
{   ngsai::RawKineticModel::operator=(other) ;
    m_bckg_model = std::move(other.m_bckg_model) ;
    return *this ;
}


bool
ngsai::NormalizedKineticModel::operator == (
                            const ngsai::NormalizedKineticModel& other) const
{   if(ngsai::RawKineticModel::operator!=(other))
    {   return false ; }

    if(m_has_kmermap != other.m_has_kmermap)
    {   return false ; }

    if(m_bckg_model != other.m_bckg_model)
    {   return false ; }

    return true ;
}

bool
ngsai::NormalizedKineticModel::operator != (
                            const ngsai::NormalizedKineticModel& other) const
{ return not ((*this) == other) ; }


void 
ngsai::NormalizedKineticModel::add(const ngsai::KineticSignal& kinetics)
{   
    // because of the normalization, the kinetics vector will contain
    // 0's on the first size_kmer_half and last size_kmer_half positions
    // so kinetic signal must contain m_size + 2*size_kmer_half values
    // and only the central values will be added
    size_t size_extended = this->getKineticSignalRequiredSize() ;
    
    if(not m_has_kmermap)
    {   throw std::runtime_error("NormalizedKineticModel error : an instance "
                                 "without KmerMap cannot be trained.") ;
    }
    if(not m_is_init)
    {   throw std::runtime_error("NormalizedKineticModel error : cannot add "
                                 "data to a model that has not been "
                                 "initialized") ;
    }
    if(kinetics.size() != size_extended)
    {   char msg[4096]  ;
        sprintf(msg,
               "NormalizedKineticModel error : number of IPD / PWD (%zu) must "
               "be %zu to match match model size after normalization (%zu)",
               kinetics.size(),
               size_extended,
               m_size) ;
        throw std::invalid_argument(msg) ;
    }
    if(not kinetics.isComplete())
    {   throw std::invalid_argument("NormalizedKineticModel error : "
                                    "KineticSignal does not have "
                                    "complete data") ;
    }

    // normalize
    auto kinetics_norm_fw = ngsai::normalize_kinetics(
                                            kinetics.getSequenceFw(),
                                            kinetics.getIPDFw(),
                                            kinetics.getPWDFw(),
                                            m_bckg_model) ;
    auto kinetics_norm_rv = ngsai::normalize_kinetics(
                                            kinetics.getSequenceRv(),
                                            kinetics.getIPDRv(),
                                            kinetics.getPWDRv(),
                                            m_bckg_model) ;
    ngsai::KineticSignal kinetics_norm(
        std::string(m_size, 'N'),  // exact seq not needed but must have a value :-)
        std::string(m_size, 'N'),  // exact seq not needed but must have a value :-)
        std::move(kinetics_norm_fw.first),     // IPDR fw
        std::move(kinetics_norm_rv.first),     // IPDR rv
        std::move(kinetics_norm_fw.second),    // PWDR fw
        std::move(kinetics_norm_rv.second)) ;  // PWDR rv

    RawKineticModel::add(kinetics_norm) ;
}


void
ngsai::NormalizedKineticModel::add(
                            const KineticModel& model)
{   model.addTo(*this) ; }


ngsai::NormalizedKineticModel*
ngsai::NormalizedKineticModel::copy() const
{
    return new ngsai::NormalizedKineticModel(*this) ;
}


double
ngsai::NormalizedKineticModel::logLikelihood(
                                const ngsai::KineticSignal& kinetics) const
{   
    // because of the normalization, the kinetics vector will contain
    // 0's on the first size_kmer_half and last size_kmer_half positions
    // so kinetic signal must contain m_size + 2*size_kmer_half values
    // and only the central values will be added
    size_t size_extended  = this->getKineticSignalRequiredSize() ;
    
    // model side errors
    if(not m_has_kmermap)
    {   throw std::runtime_error("NormalizedKineticModel error : an instance "
                                 "without KmerMap cannot compute "
                                 "loglikelihood.") ;
    }
    if(not m_is_init)
    {   throw std::runtime_error("NormalizedKineticModel error : cannot "
                                 "compute a loglikelihood from a model that "
                                 "has not been initialized") ;
    }
    if(not m_is_density)
    {   throw std::runtime_error("NormalizedKineticModel error : cannot "
                                 "compute a loglikelihood from a model that "
                                 "does not contain signal densities") ;
    }
    if(not m_is_log)
    {   throw std::runtime_error("NormalizedKineticModel error : cannot "
                                 "compute a loglikelihood from a model that "
                                 "does not contain signal log densities") ;

    }

    // cannot compute from something incomplete
    if(not kinetics.isComplete())
    {   return std::nan("") ; }

    // kinetic side errors
    if(kinetics.size() != size_extended)
    {   char msg[4096]  ;
        sprintf(msg,
               "NormalizedKineticModel error : size of IPD / PWD / sequence " 
               "(%zu) must be %zu to match match model size after "
               "normalization (%zu)",
               kinetics.size(),
               size_extended, 
               m_size) ;
        throw std::invalid_argument(msg) ;
    }

    // normalize
    auto kinetics_norm_fw = ngsai::normalize_kinetics(
                                            kinetics.getSequenceFw(),
                                            kinetics.getIPDFw(),
                                            kinetics.getPWDFw(),
                                            m_bckg_model) ;
    auto kinetics_norm_rv = ngsai::normalize_kinetics(
                                            kinetics.getSequenceRv(),
                                            kinetics.getIPDRv(),
                                            kinetics.getPWDRv(),
                                            m_bckg_model) ;
    ngsai::KineticSignal kinetics_norm(
        std::string(m_size, 'N'),  // exact seq not needed but must have a value :-)
        std::string(m_size, 'N'),  // exact seq not needed but must have a value :-)
        std::move(kinetics_norm_fw.first),
        std::move(kinetics_norm_rv.first),
        std::move(kinetics_norm_fw.second),
        std::move(kinetics_norm_rv.second)) ;

    return RawKineticModel::logLikelihood(kinetics_norm) ;
}


size_t
ngsai::NormalizedKineticModel::getKineticSignalRequiredSize() const
{   if(not this->isInit())
    {   return 0 ; }
    
    size_t size_kmer_half = m_bckg_model.getKmerSize() / 2 ;
    size_t size_extended  = m_size + 2*size_kmer_half ;
    return size_extended ;    
}


std::string
ngsai::NormalizedKineticModel::toString() const
{   return RawKineticModel::toString() ; }


void
ngsai::NormalizedKineticModel::save(const std::string& path) const
{   std::ofstream f_out(path) ;
    boost::archive::text_oarchive arch_out(f_out) ;
    arch_out << *this ;
    f_out.close() ;
}


void
ngsai::NormalizedKineticModel::load(const std::string& path)
{
    std::ifstream f_in(path) ;
    boost::archive::text_iarchive arch_in(f_in) ;
    arch_in >> *this ;
    f_in.close() ;
}




void 
ngsai::NormalizedKineticModel::addTo(
            ngsai::RawKineticModel&) const
{
   throw std::invalid_argument(
        "NormalizedKineticModel error : cannot add a "
        "NormalizedKineticModel to a RawKineticModel") ;
}


void 
ngsai::NormalizedKineticModel::addTo(
            ngsai::NormalizedKineticModel& model) const
{
    if(not m_is_init)
    {   throw std::invalid_argument(
            "NormalizedKineticModel error : this model has "
            "not been initialized") ;
    }
     if(not model.isInit())
    {   throw std::invalid_argument(
        "NormalizedKineticModel error : given model has "
        "not been initialized") ;
    }
    if(this->size() != model.size())
    {   char msg[4096] ; 
        sprintf(msg,
                "NormalizedKineticModel error : cannot "
                "add models having different sizes "
                "(%zu, %zu)",
                this->size(),
                model.size()) ;
        throw(std::invalid_argument(msg)) ;
    }

    try
    {   for(size_t i=0; i<this->size(); i++)
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
ngsai::NormalizedKineticModel::addTo(
            ngsai::PairWiseKineticModel&) const
{   throw std::invalid_argument(
        "NormalizedKineticModel error : cannot add a "
        "NormalizedKineticModel to a "
        "PairWiseKineticModel") ;
}


void
ngsai::NormalizedKineticModel::addTo(
            ngsai::PairWiseNormalizedKineticModel&) const
{   throw std::invalid_argument(
        "NormalizedKineticModel error : cannot add a "
        "NormalizedKineticModel to a "
        "PairWiseNormalizedKineticModel") ;
}


void
ngsai::NormalizedKineticModel::addTo(
            ngsai::DiPositionKineticModel&) const
{   throw std::invalid_argument(
        "NormalizedKineticModel error : cannot add a "
        "NormalizedKineticModel to a "
        "DiPositionKineticModel") ;
}


void
ngsai::NormalizedKineticModel::addTo(
            ngsai::DiPositionNormalizedKineticModel&) const
{   throw std::invalid_argument(
        "NormalizedKineticModel error : cannot add a "
        "NormalizedKineticModel to a "
        "DiPositionNormalizedKineticModel") ;
}