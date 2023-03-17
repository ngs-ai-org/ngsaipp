#include <ngsaipp/epigenetics/KineticClassifier.hpp>

#include <iostream>
#include <string>
#include <vector>
#include <list>
#include <utility>            // std::move
#include <stdexcept>          // std::invalid_argument
#include <cmath>              // std::nan()
#include <pbbam/BamRecord.h>
#include <boost/histogram.hpp>


#include <ngsaipp/epigenetics/KineticModel.hpp>
#include <ngsaipp/epigenetics/KineticSignal.hpp>
#include <ngsaipp/epigenetics/CcsKineticExtractor.hpp>
#include <ngsaipp/genome/CpGRegion.hpp>  // ngsai::genome::CpGRegion
#include <ngsaipp/genome/constants.hpp>  // ngsai::genome::strand


namespace bh = boost::histogram ;


void
print_vector(std::ostream& stream, const std::vector<double>& v)
{   for(const auto& x : v)
    {   stream << x << " " ; }
}


ngsai::KineticClassifier::KineticClassifier()
    : m_window_size(0),
      m_model_meth(nullptr),
      m_model_unmeth(nullptr)
{ ; }


ngsai::KineticClassifier::KineticClassifier(
                                ngsai::KineticModel* model_meth,
                                ngsai::KineticModel* model_unmeth)
    : KineticClassifier()
{   this->setModels(model_meth, model_unmeth) ;
}


ngsai::KineticClassifier::KineticClassifier(
                                const ngsai::KineticClassifier& other)
    : m_window_size(other.m_window_size),
      m_model_meth(nullptr),
      m_model_unmeth(nullptr)
{
    m_model_meth   = other.m_model_meth->copy() ;
    m_model_unmeth = other.m_model_unmeth->copy() ;
}


ngsai::KineticClassifier::KineticClassifier(ngsai::KineticClassifier&& other)
    : m_window_size (std::move(other.m_window_size)),
      m_model_meth(std::move(other.m_model_meth)),
      m_model_unmeth(std::move(other.m_model_unmeth))
{ ; }


ngsai::KineticClassifier::~KineticClassifier()
{   this->freeModels() ; }


ngsai::KineticClassifier& 
ngsai::KineticClassifier::operator = (const KineticClassifier& other)
{   this->freeModels() ;
    
    m_window_size  = other.m_window_size ;
    m_model_meth   = other.m_model_meth->copy() ;
    m_model_unmeth = other.m_model_unmeth->copy() ;
    return *this ;
}


ngsai::KineticClassifier& 
ngsai::KineticClassifier::operator = (KineticClassifier&& other)
{   this->freeModels() ;

    m_window_size  = std::move(other.m_window_size) ;
    m_model_meth   = std::move(other.m_model_meth) ;
    m_model_unmeth = std::move(other.m_model_unmeth) ;
    return *this ;
}


void
ngsai::KineticClassifier::setModels(ngsai::KineticModel* model_meth,
                                    ngsai::KineticModel* model_unmeth)
{   
    this->freeModels() ;
    m_model_meth   = model_meth ;
    m_model_unmeth = model_unmeth ;

    // check models and turns them into log density
    this->checkModels() ;

    // m_window_size = model_meth->size() ;
    m_window_size = model_meth->getKineticSignalRequiredSize() ;
}



std::pair<double,double>
ngsai::KineticClassifier::classify(
                    const ngsai::genome::CpGRegion & cpg,
                    const std::list<PacBio::BAM::BamRecord>& ccss,
                    double prob_meth,
                    double prob_unmeth) const
{   
    // the models were not given
    if(m_window_size == 0)
    {   std::string msg("classification is not possible, no model was given") ;
        throw std::runtime_error(msg.c_str()) ;
    }
    
    if(cpg.size() != 2)
    {   char msg[256] ;
        sprintf(msg, 
                "CpG coordinates don't fit a CpG in size (%s)",
                cpg.toString().c_str()) ;
        throw std::invalid_argument(msg) ;
    }

    double prob_sum = prob_meth + prob_unmeth ; 
    prob_meth       = prob_meth   / prob_sum ;
    prob_unmeth     = prob_unmeth / prob_sum ;

    ngsai::CcsKineticExtractor kinetic_extractor ;

    size_t win_size_half = m_window_size / 2 ;

    // region from which the kinetics should be extracted
    // on + strand C corresponds to start pos of bed entry
    ngsai::genome::GenomeRegion cpg_p = cpg.forwardCoordinates() ;
    cpg_p.start -= win_size_half ;
    cpg_p.end   += win_size_half - 1 ;
    // on - strand C corresponds to end-1 pos of bed entry
    ngsai::genome::GenomeRegion cpg_m = cpg.reverseCoordinates() ;
    cpg_m.start -= win_size_half - 1 ;
    cpg_m.end   += win_size_half ;

    // compute average IPD and PWD vectors
    double n_fw = 0. ;
    double n_rv = 0. ;
    std::string seq_fw ;
    std::string seq_rv ;
    std::vector<double> ipd_fw_m(m_window_size, 0.) ;
    std::vector<double> pwd_fw_m(m_window_size, 0.) ;   
    std::vector<double> ipd_rv_m(m_window_size, 0.) ;
    std::vector<double> pwd_rv_m(m_window_size, 0.) ;  


    std::vector<uint16_t> ipd_16 ;
    std::vector<uint16_t> pwd_16 ;

    for(const auto& ccs : ccss)
    {   // forward strand
        if(kinetic_extractor.extract(ccs, cpg_p))
        {   if(n_fw == 0.)
            {   seq_fw = kinetic_extractor.getSequence() ; }
            ipd_16 = kinetic_extractor.getIPD() ;
            pwd_16 = kinetic_extractor.getPWD() ;
            for(size_t i=0; i<ipd_16.size(); i++)
            {   ipd_fw_m[i] += static_cast<double>(ipd_16[i]) ;
                pwd_fw_m[i] += static_cast<double>(pwd_16[i]) ;
            }
            n_fw += 1. ;
        }
        // reverse strand
        if(kinetic_extractor.extract(ccs, cpg_m))
        {   if(n_rv == 0.)
            {   seq_rv = kinetic_extractor.getSequence() ; }
            ipd_16 = kinetic_extractor.getIPD() ;
            pwd_16 = kinetic_extractor.getPWD() ;
            for(size_t i=0; i<ipd_16.size(); i++)
            {   ipd_rv_m[i] += static_cast<double>(ipd_16[i]) ;
                pwd_rv_m[i] += static_cast<double>(pwd_16[i]) ;
            }
            n_rv += 1 ;
        }
    }
    
    ngsai::KineticSignal kinetics ;
    if(n_fw != 0.)
    {   for(size_t i=0; i<m_window_size; i++)
        {   ipd_fw_m[i] /= n_fw ;
            pwd_fw_m[i] /= n_fw ;
        }
        kinetics.setIPDFw(std::move(ipd_fw_m)) ;
        kinetics.setPWDFw(std::move(pwd_fw_m)) ;
        kinetics.setSequenceFw(std::move(seq_fw)) ;
    }
    if(n_rv != 0.)
    {   for(size_t i=0; i<m_window_size; i++)
        {   ipd_rv_m[i] /= n_rv ;
            pwd_rv_m[i] /= n_rv ;
        }
        kinetics.setIPDRv(std::move(ipd_rv_m)) ;
        kinetics.setPWDRv(std::move(pwd_rv_m)) ;
        kinetics.setSequenceRv(std::move(seq_rv)) ;
    }
    

    // compute posterior probabilities
    // return nan if IPD and/or PWD are missing
    try
    {   return this->classify(kinetics, prob_meth, prob_unmeth) ; }
    catch(std::exception& e)
    {   std::cerr << "KineticClassifier error : error while classifying" 
                  << std::endl
                  << e.what()
                  << std::endl ;
        std::cerr << "region :" << cpg << std::endl ;
        std::cerr << "seq fw: " << kinetics.getSequenceFw() 
                  << std::endl ;
        std::cerr << "seq rv: " << kinetics.getSequenceRv()
                  << std::endl ;
        std::cerr << "IPD fw: " ; print_vector(std::cerr, 
                                               kinetics.getIPDFw()) ;
        std::cerr << std::endl ;
        std::cerr << "IPD rv: " ; print_vector(std::cerr,
                                               kinetics.getIPDRv()) ;
        std::cerr << std::endl ;
        std::cerr << "PWD fw: " ; print_vector(std::cerr,
                                               kinetics.getPWDFw()) ;
        std::cerr << std::endl ;
        std::cerr << "PWD rv: " ; print_vector(std::cerr,
                                               kinetics.getPWDRv()) ;
        std::cerr << std::endl ;
        throw e ;
    }
}


std::pair<double, double>
ngsai::KineticClassifier::classify(const ngsai::KineticSignal& kinetics,
                                   double prob_meth,
                                   double prob_unmeth) const
{   // the models were not given
    if(m_window_size == 0)
    {   std::string msg("KineticClassifier error : classification is not "
                        "possible, no model was given") ;
        throw std::runtime_error(msg) ;
    }
    
    // don't check kinetics here, let the models check what they need
    // if(kinetics.size() != m_window_size)
    // {   char msg[4096] ;
    //     sprintf(msg,
    //             "KineticClassifier error : number of IPDs / PWDs (%zu) do not "
    //             "fit the window size (%zu)",
    //             kinetics.size(), 
    //             m_window_size) ;
    //     throw std::invalid_argument(msg) ;
    // }
    
    // normalize prior prob
    double prob_sum = prob_meth + prob_unmeth ; 
    prob_meth   = prob_meth   / prob_sum ;
    prob_unmeth = prob_unmeth / prob_sum ;
    
    // compute loglikelihoods
    // will return nan if no IPD and/or PWD
    std::pair<double,double> ll ;
    // meth
    ll.first  = m_model_meth->logLikelihood(kinetics) ;
    // unmeth
    ll.second = m_model_unmeth->logLikelihood(kinetics) ;

    // std::cerr << ll.first << "\t" << ll.second << std::endl ;

    // rescale log likelihood
    double c = std::max(ll.first, ll.second) ;
    ll.first  -= c ;
    ll.second -= c ;

    // std::cerr << ll.first << "\t" << ll.second << std::endl ;
    // std::cerr << "----------------------------------------------------" << std::endl ;

    // compute posterior probabilities
    double post_prob_meth   = std::exp(ll.first)  * prob_meth ;
    double post_prob_unmeth = std::exp(ll.second) * prob_unmeth ;
    double post_prob_sum    = post_prob_meth + post_prob_unmeth ;
    post_prob_meth /= post_prob_sum ;
    post_prob_unmeth /= post_prob_sum ;

    return std::make_pair(post_prob_meth, post_prob_unmeth) ;
}


void
ngsai::KineticClassifier::checkModels() const
{   
    // check sizes are the same
    if(m_model_meth->size() != m_model_unmeth->size())
    {   char msg[4096] ;
        sprintf(msg, "the model sizes must be equal (%zu, %zu)",
                m_model_meth->size(),
                m_model_unmeth->size()) ;
        throw std::invalid_argument(msg) ; 
    }

    // check all models are init
    if((not m_model_meth->isInit()) or
       (not m_model_unmeth->isInit()))
    {
        std::string msg("at least one of the model has not been initialised") ;
        throw std::runtime_error(msg.c_str()) ;
    }
    
    // check all models have been normalized to densities
    if((not m_model_meth->isDensity()) or
       (not m_model_unmeth->isDensity()))
    {
        std::string msg("at least one of the model does not contain a "
                        "density") ;
        throw std::runtime_error(msg.c_str()) ;
    }

    // check that all models contains log of prob densities
    if((not m_model_meth->isLog()) or
       (not m_model_unmeth->isLog()))
    {
        std::string msg("at least one of the model does not contain log "
                         "probabilities") ;
        throw std::runtime_error(msg.c_str()) ;
    }
}


void
ngsai::KineticClassifier::freeModels()
{
    if(m_model_meth != nullptr)
    {   delete m_model_meth ; }    
    m_model_meth = nullptr ;
    if(m_model_unmeth != nullptr)
    {   delete m_model_unmeth ;  } 
    m_model_unmeth = nullptr ;
    m_window_size = 0 ;
}