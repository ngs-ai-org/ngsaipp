#include <ngsaipp/epigenetics/KineticModel.hpp>
#include <utility>


ngsai::KineticModel::KineticModel()
    : m_size(0),
      m_is_init(false),
      m_is_density(false),
      m_is_log(false)
{ ; }


ngsai::KineticModel::KineticModel(size_t size,
                                  bool is_init,
                                  bool is_density,
                                  bool is_log)
    : m_size(size),
      m_is_init(is_init),
      m_is_density(is_density),
      m_is_log(is_log)
{ ; }


ngsai::KineticModel::KineticModel(const ngsai::KineticModel& other)
    : m_size(other.m_size),
      m_is_init(other.m_is_init),
      m_is_density(other.m_is_density),
      m_is_log(other.m_is_log)
{ ; }


ngsai::KineticModel::KineticModel(ngsai::KineticModel&& other)
    : m_size(std::move(other.m_size)),
      m_is_init(std::move(other.m_is_init)),
      m_is_density(std::move(other.m_is_density)),
      m_is_log(std::move(other.m_is_log))
{ ; }


ngsai::KineticModel::~KineticModel()
{ ; }


bool
ngsai::KineticModel::operator == (const ngsai::KineticModel& other) const
{
  if((m_size       == other.m_size ) and
     (m_is_init    == other.m_is_init) and
     (m_is_density == other.m_is_density) and
     (m_is_log     == other.m_is_log))
  { return true ; }
  return false ;
}

bool
ngsai::KineticModel::operator != (const ngsai::KineticModel& other) const
{ return not ((*this) == other) ; }

bool
ngsai::KineticModel::isInit() const
{   return m_is_init ; }


bool
ngsai::KineticModel::isDensity() const
{   return m_is_density ; }


bool
ngsai::KineticModel::isLog() const
{   return m_is_log ; }


size_t
ngsai::KineticModel::size() const
{   return m_size ; }


size_t
ngsai::KineticModel::getKineticSignalRequiredSize() const
{   if(not this->isInit())
    {   return 0 ; }
    return m_size ;
}