#include <ngsaipp/epigenetics/DiPositionNormalizedKineticModel.hpp>

#include <vector>
#include <string>
#include <fstream>
#include <stdexcept>

#include <ngsaipp/epigenetics/KmerMap.hpp>


ngsai::DiPositionNormalizedKineticModel::DiPositionNormalizedKineticModel()
    : PairWiseNormalizedKineticModel()
{ ; }


ngsai::DiPositionNormalizedKineticModel::DiPositionNormalizedKineticModel(
                        const ngsai::KmerMap& bckg_model)
    : PairWiseNormalizedKineticModel(bckg_model)
{ ; }


ngsai::DiPositionNormalizedKineticModel::DiPositionNormalizedKineticModel(
                        const ngsai::DiPositionNormalizedKineticModel& other)
    : PairWiseNormalizedKineticModel(other)
{ ; }


ngsai::DiPositionNormalizedKineticModel::DiPositionNormalizedKineticModel(
                        ngsai::DiPositionNormalizedKineticModel&& other)
    : PairWiseNormalizedKineticModel(other)
{ ; }


ngsai::DiPositionNormalizedKineticModel::~DiPositionNormalizedKineticModel()
{ ; }


ngsai::DiPositionNormalizedKineticModel&
ngsai::DiPositionNormalizedKineticModel::operator = (
                        const ngsai::DiPositionNormalizedKineticModel& other)
{   ngsai::PairWiseNormalizedKineticModel::operator=(other) ;
    return *this ;
}


ngsai::DiPositionNormalizedKineticModel&
ngsai::DiPositionNormalizedKineticModel::operator = (
                        ngsai::DiPositionNormalizedKineticModel&& other)
{   ngsai::PairWiseNormalizedKineticModel::operator=(other) ;
    return *this ;
}

bool
ngsai::DiPositionNormalizedKineticModel::operator == (
                    const ngsai::DiPositionNormalizedKineticModel& other) const
{   return ngsai::PairWiseNormalizedKineticModel::operator==(other) ; }

bool
ngsai::DiPositionNormalizedKineticModel::operator != (
                    const ngsai::DiPositionNormalizedKineticModel& other) const
{   return ngsai::PairWiseNormalizedKineticModel::operator!=(other) ; }


void
ngsai::DiPositionNormalizedKineticModel::add(
                                const KineticModel& model)
{   model.addTo(*this) ; }


ngsai::DiPositionNormalizedKineticModel*
ngsai::DiPositionNormalizedKineticModel::copy() const
{
    return new ngsai::DiPositionNormalizedKineticModel(*this) ;
}


void
ngsai::DiPositionNormalizedKineticModel::computeMeaningfulPairs()
{   m_pair_nb = m_size - 1 ;
    m_pair_index.clear() ;
    for(size_t i=0; i<m_size-1; i++)
    {   m_pair_index.emplace_back(std::make_pair(i,i+1)) ; }
}


void 
ngsai::DiPositionNormalizedKineticModel::addTo(
            ngsai::RawKineticModel&) const
{
   throw std::invalid_argument(
        "DiPositionNormalizedKineticModel error : cannot "
        "add a DiPositionNormalizedKineticModel to a "
        "RawKineticModel") ;
}


void
ngsai::DiPositionNormalizedKineticModel::addTo(
            ngsai::NormalizedKineticModel&) const
{   throw std::invalid_argument(
        "DiPositionNormalizedKineticModel error : cannot "
        "add a DiPositionNormalizedKineticModel to a"
        "NormalizedKineticModel") ;
}


void
ngsai::DiPositionNormalizedKineticModel::addTo(
            ngsai::PairWiseKineticModel&) const
{   throw std::invalid_argument(
        "DiPositionNormalizedKineticModel error : cannot "
        "add a DiPositionNormalizedKineticModel to a "
        "PairWiseKineticModel") ;
}


void
ngsai::DiPositionNormalizedKineticModel::addTo(
            ngsai::PairWiseNormalizedKineticModel&) const
{   throw std::invalid_argument(
        "DiPositionNormalizedKineticModel error : cannot "
        "add a DiPositionNormalizedKineticModel to a "
        "PairWiseNormalizedKineticModel") ;
}


void
ngsai::DiPositionNormalizedKineticModel::addTo(
            ngsai::DiPositionKineticModel&) const
{   throw std::invalid_argument(
        "DiPositionNormalizedKineticModel error : cannot "
        "add a DiPositionNormalizedKineticModel to a "
        "DiPositionKineticModel") ;
}


void 
ngsai::DiPositionNormalizedKineticModel::addTo(
    ngsai::DiPositionNormalizedKineticModel& model) const
{
    if(not m_is_init)
    {   throw std::invalid_argument(
            "DiPositionNormalizedKineticModel error : this "
            "model has not been initialized") ;
    }
     if(not model.isInit())
    {   throw std::invalid_argument(
        "DiPositionNormalizedKineticModel error : given "
        "model has not been initialized") ;
    }
    if(this->size() != model.size())
    {   char msg[4096] ; 
        sprintf(msg,
                "DiPositionNormalizedKineticModel error : "
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