#include <ngsaipp/epigenetics/DiPositionKineticModel.hpp>

#include <utility>
#include <stdexcept>                   // std::invalid_argument, std::runtime_error, std::out_of_range
#include <set>                         // std::set
#include <map>                         // std::map

#include <ngsaipp/io/bed_io.hpp>                        // ngsai::BedReader, ngsai::BedRecord
#include <ngsaipp/genome/constants.hpp>                 // ngsai::genome::strand
#include <ngsaipp/epigenetics/CcsKineticExtractor.hpp>  // ngsai::CcsKineticExtractor
#include <ngsaipp/epigenetics/KineticSignal.hpp>        // ngsai::KineticSignal


namespace bh = boost::histogram ;



ngsai::DiPositionKineticModel::DiPositionKineticModel()
    : PairWiseKineticModel()
{ ; }


ngsai::DiPositionKineticModel::DiPositionKineticModel(
                    const ngsai::DiPositionKineticModel& other)
    : PairWiseKineticModel(other)
{ ; }


ngsai::DiPositionKineticModel::DiPositionKineticModel(
                    ngsai::DiPositionKineticModel&& other)
    : PairWiseKineticModel(other)
{ ; }


ngsai::DiPositionKineticModel::~DiPositionKineticModel()
{ ; }


ngsai::DiPositionKineticModel&
ngsai::DiPositionKineticModel::operator = 
                    (const ngsai::DiPositionKineticModel& other)
{   PairWiseKineticModel::operator=(other) ;
    return *this ;
}


ngsai::DiPositionKineticModel&
ngsai::DiPositionKineticModel::operator = 
                    (ngsai::DiPositionKineticModel&& other)
{   PairWiseKineticModel::operator=(other) ;
    return *this ;
}


bool
ngsai::DiPositionKineticModel::operator == (
                    const ngsai::DiPositionKineticModel& other) const
{   return ngsai::PairWiseKineticModel::operator==(other) ; }


bool
ngsai::DiPositionKineticModel::operator != (
                    const ngsai::DiPositionKineticModel& other) const
{   return ngsai::PairWiseKineticModel::operator!=(other) ; }


void
ngsai::DiPositionKineticModel::add(
                            const KineticModel& model)
{   model.addTo(*this) ; }


ngsai::DiPositionKineticModel*
ngsai::DiPositionKineticModel::copy() const
{
    return new ngsai::DiPositionKineticModel(*this) ;
}


void
ngsai::DiPositionKineticModel::computeMeaningfulPairs()
{   m_pair_nb = m_size - 1 ;
    m_pair_index.clear() ;
    for(size_t i=0; i<m_size-1; i++)
    {   m_pair_index.emplace_back(std::make_pair(i,i+1)) ; }
}


void 
ngsai::DiPositionKineticModel::addTo(
            ngsai::RawKineticModel&) const
{
   throw std::invalid_argument(
        "DiPositionKineticModel error : cannot add "
        "a DiPositionKineticModel to a "
        "RawKineticModel") ;
}


void
ngsai::DiPositionKineticModel::addTo(
            ngsai::NormalizedKineticModel&) const
{   throw std::invalid_argument(
        "DiPositionKineticModel error : cannot add "
        "a DiPositionKineticModel to a"
        "NormalizedKineticModel") ;
}


void
ngsai::DiPositionKineticModel::addTo(
                ngsai::PairWiseKineticModel&) const
{   throw std::invalid_argument(
        "DiPositionKineticModel error : cannot add a "
        "DiPositionKineticModel to a "
        "PairWiseKineticModel") ;
}


void
ngsai::DiPositionKineticModel::addTo(
            ngsai::PairWiseNormalizedKineticModel&) const
{   throw std::invalid_argument(
        "DiPositionKineticModel error : cannot add "
        "a DiPositionKineticModel to a "
        "PairWiseNormalizedKineticModel") ;
}


void 
ngsai::DiPositionKineticModel::addTo(
            ngsai::DiPositionKineticModel& model) const
{
    if(not m_is_init)
    {   throw std::invalid_argument(
            "DiPositionKineticModel error : this "
            "model has not been initialized") ;
    }
     if(not model.isInit())
    {   throw std::invalid_argument(
        "DiPositionKineticModel error : given "
        "model has not been initialized") ;
    }
    if(this->size() != model.size())
    {   char msg[4096] ; 
        sprintf(msg,
                "DiPositionKineticModel error : "
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
ngsai::DiPositionKineticModel::addTo(
            ngsai::DiPositionNormalizedKineticModel&) const
{   throw std::invalid_argument(
        "DiPositionKineticModel error : cannot add a "
        "DiPositionKineticModel to a "
        "DiPositionNormalizedKineticModel") ;
}