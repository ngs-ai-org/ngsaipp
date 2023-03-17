#include <ngsaipp/epigenetics/PairWiseNormalizedKineticModel.hpp>

#include <vector>
#include <string>
#include <fstream>
#include <stdexcept>
#include <boost/archive/text_iarchive.hpp>      // boost::archive::text_iarchive
#include <boost/archive/text_oarchive.hpp>      // boost::archive::text_oarchive
#include <boost/serialization/utility.hpp>      // stl data structures serialization

#include <ngsaipp/epigenetics/KmerMap.hpp>
#include <ngsaipp/epigenetics/model_utility.hpp>  // ngsai::normalize_kinetics()


ngsai::PairWiseNormalizedKineticModel::PairWiseNormalizedKineticModel()
    : PairWiseKineticModel(),
      m_has_kmermap(false),
      m_bckg_model()
{ ; }


ngsai::PairWiseNormalizedKineticModel::PairWiseNormalizedKineticModel(
                        const ngsai::KmerMap& bckg_model)
    : PairWiseKineticModel(),
      m_has_kmermap(true),
      m_bckg_model(bckg_model)
{ ; }


ngsai::PairWiseNormalizedKineticModel::PairWiseNormalizedKineticModel(
                        const ngsai::PairWiseNormalizedKineticModel& other)
    : PairWiseKineticModel(other),
      m_has_kmermap(other.m_has_kmermap),
      m_bckg_model(other.m_bckg_model)
{ ; }


ngsai::PairWiseNormalizedKineticModel::PairWiseNormalizedKineticModel(
                        ngsai::PairWiseNormalizedKineticModel&& other)
    : PairWiseKineticModel(other),
      m_has_kmermap(std::move(other.m_has_kmermap)),
      m_bckg_model(std::move(other.m_bckg_model))
{ ; }


ngsai::PairWiseNormalizedKineticModel::~PairWiseNormalizedKineticModel()
{ ; }


ngsai::PairWiseNormalizedKineticModel&
ngsai::PairWiseNormalizedKineticModel::operator = (
                        const ngsai::PairWiseNormalizedKineticModel& other)
{   ngsai::PairWiseKineticModel::operator=(other) ;
    m_bckg_model = other.m_bckg_model ;
    m_has_kmermap = other.m_has_kmermap ;
    return *this ;
}


ngsai::PairWiseNormalizedKineticModel&
ngsai::PairWiseNormalizedKineticModel::operator = (
                        ngsai::PairWiseNormalizedKineticModel&& other)
{   ngsai::PairWiseKineticModel::operator=(other) ;
    m_bckg_model = std::move(other.m_bckg_model) ;
    m_has_kmermap = other.m_has_kmermap ;
    return *this ;
}

bool
ngsai::PairWiseNormalizedKineticModel::operator == (
                    const ngsai::PairWiseNormalizedKineticModel& other) const
{   if(ngsai::PairWiseKineticModel::operator!=(other))
    {   return false ; }

   if(m_has_kmermap != other.m_has_kmermap)
   {    return false ; }

   if(m_bckg_model != other.m_bckg_model)
   {    return false ; }

    return true ;
}


bool
ngsai::PairWiseNormalizedKineticModel::operator != (
                    const ngsai::PairWiseNormalizedKineticModel& other) const
{   return not ((*this) == other) ; }


void 
ngsai::PairWiseNormalizedKineticModel::add(
                      const ngsai::KineticSignal& kinetics)
{   
    // because of the normalization, the kinetics vector will contain
    // 0's on the first size_kmer_half and last size_kmer_half positions
    // so kinetic signal must contain m_size + 2*size_kmer_half values
    // and only the central values will be added
    size_t size_extended = this->getKineticSignalRequiredSize() ;
    
    if(not m_has_kmermap)
    {   throw std::runtime_error("PairWiseNormalizedKineticModel error : an "
                                 "instance without KmerMap cannot be "
                                 "trained.") ;
    }
    if(not m_is_init)
    {   throw std::runtime_error("PairWiseNormalizedKineticModel error : "
                                 "cannot add data to a model that has not "
                                 "been initialized") ;
    }
    if(kinetics.size() != size_extended)
    {   char msg[4096]  ;
        sprintf(msg,
               "PairWiseNormalizedKineticModel error : number of IPD / PWD "
               "(%zu) must be %zu to match match model size after "
               "normalization (%zu)",
               kinetics.size(),
               size_extended,
               m_size) ;
        throw std::invalid_argument(msg) ;
    }
    if(not kinetics.isComplete())
    {   throw std::invalid_argument("PairWiseNormalizedKineticModel error : "
                                    "KineticSignal does not have complete "
                                    "data") ;
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

    PairWiseKineticModel::add(kinetics_norm) ;
}


void
ngsai::PairWiseNormalizedKineticModel::add(
                            const KineticModel& model)
{   model.addTo(*this) ; }


ngsai::PairWiseNormalizedKineticModel*
ngsai::PairWiseNormalizedKineticModel::copy() const
{
    return new ngsai::PairWiseNormalizedKineticModel(*this) ;
}


double
ngsai::PairWiseNormalizedKineticModel::logLikelihood(
                                const ngsai::KineticSignal& kinetics) const
{   
    // because of the normalization, the kinetics vector will contain
    // 0's on the first size_kmer_half and last size_kmer_half positions
    // so kinetic signal must contain m_size + 2*size_kmer_half values
    // and only the central values will be added
    size_t size_extended  = this->getKineticSignalRequiredSize() ;
    
    // model side errors
    if(not m_has_kmermap)
    {   throw std::runtime_error("PairWiseNormalizedKineticModel error : an "
                                 "instance without KmerMap cannot compute "
                                 "loglikelihood.") ;
    }
    if(not m_is_init)
    {   throw std::runtime_error("PairWiseNormalizedKineticModel error : "
                                 "cannot compute a loglikelihood from a model "
                                 "that has not been initialized") ;
    }
    if(not m_is_density)
    {   throw std::runtime_error("PairWiseNormalizedKineticModel error : "
                                 "cannot compute a loglikelihood from a model "
                                 "that does not contain signal densities") ;
    }
    if(not m_is_log)
    {   throw std::runtime_error("PairWiseNormalizedKineticModel error : "
                                 "cannot compute a loglikelihood from a model "
                                 "that does not contain signal log "
                                 "densities") ;

    }

    // cannot compute from something incomplete
    if(not kinetics.isComplete())
    {   return std::nan("") ; }

    // kinetic side errors
    if(kinetics.size() != size_extended)
    {   char msg[4096]  ;
        sprintf(msg,
               "PairWiseNormalizedKineticModel error : size of IPD / PWD / sequence " 
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

    return PairWiseKineticModel::logLikelihood(kinetics_norm) ;
}


size_t
ngsai::PairWiseNormalizedKineticModel::getKineticSignalRequiredSize() const
{   if(not this->isInit())
    {   return 0 ; }
    
    size_t size_kmer_half = m_bckg_model.getKmerSize() / 2 ;
    size_t size_extended  = m_size + 2*size_kmer_half ;
    return size_extended ;    
}

void
ngsai::PairWiseNormalizedKineticModel::save(const std::string& path) const
{   std::ofstream f_out(path) ;
    boost::archive::text_oarchive arch_out(f_out) ;
    arch_out << *this ;
    f_out.close() ;
}


void
ngsai::PairWiseNormalizedKineticModel::load(const std::string& path)
{
    std::ifstream f_in(path) ;
    boost::archive::text_iarchive arch_in(f_in) ;
    arch_in >> *this ;
    f_in.close() ;
}


void 
ngsai::PairWiseNormalizedKineticModel::addTo(
            ngsai::RawKineticModel&) const
{   throw std::invalid_argument(
        "PairWiseNormalizedKineticModel error : cannot add "
        "a PairWiseNormalizedKineticModel to a "
        "RawKineticModel") ;
}


void
ngsai::PairWiseNormalizedKineticModel::addTo(
            ngsai::NormalizedKineticModel&) const
{   throw std::invalid_argument(
        "PairWiseNormalizedKineticModel error : cannot add "
        "a PairWiseNormalizedKineticModel to a"
        "NormalizedKineticModel") ;
}


void
ngsai::PairWiseNormalizedKineticModel::addTo(
            ngsai::PairWiseKineticModel&) const
{   throw std::invalid_argument(
        "PairWiseNormalizedKineticModel error : cannot add a "
        "PairWiseNormalizedKineticModel to a "
        "PairWiseKineticModel") ;
}


void 
ngsai::PairWiseNormalizedKineticModel::addTo(
        ngsai::PairWiseNormalizedKineticModel& model) const
{
    if(not m_is_init)
    {   throw std::invalid_argument(
            "PairWiseNormalizedKineticModel error : this "
            "model has not been initialized") ;
    }
     if(not model.isInit())
    {   throw std::invalid_argument(
        "PairWiseNormalizedKineticModel error : given "
        "model has not been initialized") ;
    }
    if(this->size() != model.size())
    {   char msg[4096] ; 
        sprintf(msg,
                "PairWiseNormalizedKineticModel error : "
                "cannot add models having different sizes "
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
ngsai::PairWiseNormalizedKineticModel::addTo(
            ngsai::DiPositionKineticModel&) const
{   throw std::invalid_argument(
        "PairWiseNormalizedKineticModel error : cannot add "
        "a PairWiseNormalizedKineticModel to a "
        "DiPositionKineticModel") ;
}


void
ngsai::PairWiseNormalizedKineticModel::addTo(
            ngsai::DiPositionNormalizedKineticModel&) const
{   throw std::invalid_argument(
        "PairWiseNormalizedKineticModel error : cannot add a "
        "PairWiseNormalizedKineticModel to a "
        "DiPositionNormalizedKineticModel") ;
}
