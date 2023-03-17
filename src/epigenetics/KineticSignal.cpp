#include <ngsaipp/epigenetics/KineticSignal.hpp>

#include <vector>
#include <string>
#include <stdexcept>
#include <utility>



ngsai::KineticSignal::KineticSignal()
    : m_size(0),
      m_sequence_fw(),
      m_sequence_rv(),
      m_ipd_fw(),
      m_ipd_rv(),
      m_pwd_fw(),
      m_pwd_rv(),
      m_seq_fw_set(false),
      m_seq_rv_set(false),
      m_ipd_fw_set(false),
      m_ipd_rv_set(false),
      m_pwd_fw_set(false),
      m_pwd_rv_set(false)
{ ; }


ngsai::KineticSignal::KineticSignal(const std::string& seq_fw,
                                    const std::string& seq_rv,
                                    const std::vector<double>& ipd_fw,
                                    const std::vector<double>& ipd_rv,
                                    const std::vector<double>& pwd_fw,
                                    const std::vector<double>& pwd_rv)
    : m_size(seq_fw.size()),
      m_sequence_fw(seq_fw),
      m_sequence_rv(seq_rv),
      m_ipd_fw(ipd_fw),
      m_ipd_rv(ipd_rv),
      m_pwd_fw(pwd_fw),
      m_pwd_rv(pwd_rv),
      m_seq_fw_set(true),
      m_seq_rv_set(true),
      m_ipd_fw_set(true),
      m_ipd_rv_set(true),
      m_pwd_fw_set(true),
      m_pwd_rv_set(true)
{   
    if(m_sequence_fw.size() != m_sequence_rv.size())
    {   char msg[4096] ;
        sprintf(msg,
               "KineticSignal error : fw and rv sequences have different "
               "lengths (%zu, %zu)",
                m_sequence_fw.size(), m_sequence_rv.size()) ;
        throw std::invalid_argument(msg) ;
    }
    if(m_ipd_fw.size() != m_ipd_rv.size())
    {   char msg[4096] ;
        sprintf(msg,
               "KineticSignal error : fw and rv IPDs have different "
               "lengths (%zu, %zu)",
                m_ipd_fw.size(), m_ipd_rv.size()) ;
        throw std::invalid_argument(msg) ;
    }
    if(m_pwd_fw.size() != m_pwd_rv.size())
    {   char msg[4096] ;
        sprintf(msg,
               "KineticSignal error : fw and rv PWDs have different "
               "lengths (%zu, %zu)",
                m_pwd_fw.size(), m_pwd_rv.size()) ;
        throw std::invalid_argument(msg) ;
    }
    if(m_sequence_fw.size() != m_ipd_fw.size())
    {   char msg[4096] ;
        sprintf(msg,
               "KineticSignal error : sequences and IPDs have different "
               "lengths (%zu, %zu)",
                m_sequence_fw.size(), m_ipd_fw.size()) ;
        throw std::invalid_argument(msg) ;
    }
    if(m_sequence_fw.size() != pwd_fw.size())
    {   char msg[4096] ;
        sprintf(msg,
               "KineticSignal error : sequences and PWDs have different "
               "lengths (%zu, %zu)",
                m_sequence_fw.size(), m_pwd_fw.size()) ;
        throw std::invalid_argument(msg) ;
    }
}


ngsai::KineticSignal::KineticSignal(std::string&& seq_fw,
                                    std::string&& seq_rv,
                                    std::vector<double>&& ipd_fw,
                                    std::vector<double>&& ipd_rv,
                                    std::vector<double>&& pwd_fw,
                                    std::vector<double>&& pwd_rv)
    : m_size(seq_fw.size()),
      m_sequence_fw(seq_fw),
      m_sequence_rv(seq_rv),
      m_ipd_fw(ipd_fw),
      m_ipd_rv(ipd_rv),
      m_pwd_fw(pwd_fw),
      m_pwd_rv(pwd_rv),
      m_seq_fw_set(true),
      m_seq_rv_set(true),
      m_ipd_fw_set(true),
      m_ipd_rv_set(true),
      m_pwd_fw_set(true),
      m_pwd_rv_set(true)
{   if(m_sequence_fw.size() != m_sequence_rv.size())
    {   char msg[4096] ;
        sprintf(msg,
               "KineticSignal error : fw and rv sequences have different "
               "lengths (%zu, %zu)",
                m_sequence_fw.size(), m_sequence_rv.size()) ;
        throw std::invalid_argument(msg) ;
    }
    if(m_ipd_fw.size() != m_ipd_rv.size())
    {   char msg[4096] ;
        sprintf(msg,
               "KineticSignal error : fw and rv IPDs have different "
               "lengths (%zu, %zu)",
                m_ipd_fw.size(), m_ipd_rv.size()) ;
        throw std::invalid_argument(msg) ;
    }
    if(m_pwd_fw.size() != m_pwd_rv.size())
    {   char msg[4096] ;
        sprintf(msg,
               "KineticSignal error : fw and rv PWDs have different "
               "lengths (%zu, %zu)",
                m_pwd_fw.size(), m_pwd_rv.size()) ;
        throw std::invalid_argument(msg) ;
    }
    if(m_sequence_fw.size() != m_ipd_fw.size())
    {   char msg[4096] ;
        sprintf(msg,
               "KineticSignal error : sequences and IPDs have different "
               "lengths (%zu, %zu)",
                m_sequence_fw.size(), m_ipd_fw.size()) ;
        throw std::invalid_argument(msg) ;
    }
    if(m_sequence_fw.size() != pwd_fw.size())
    {   char msg[4096] ;
        sprintf(msg,
               "KineticSignal error : sequences and PWDs have different "
               "lengths (%zu, %zu)",
                m_sequence_fw.size(), m_pwd_fw.size()) ;
        throw std::invalid_argument(msg) ;
    }
}


ngsai::KineticSignal::KineticSignal(const ngsai::KineticSignal& other)
    : m_size(other.m_size),
      m_sequence_fw(other.m_sequence_fw),
      m_sequence_rv(other.m_sequence_rv),
      m_ipd_fw(other.m_ipd_fw),
      m_ipd_rv(other.m_ipd_rv),
      m_pwd_fw(other.m_pwd_fw),
      m_pwd_rv(other.m_pwd_rv),
      m_seq_fw_set(other.m_seq_fw_set),
      m_seq_rv_set(other.m_seq_rv_set),
      m_ipd_fw_set(other.m_ipd_fw_set),
      m_ipd_rv_set(other.m_ipd_rv_set),
      m_pwd_fw_set(other.m_pwd_fw_set),
      m_pwd_rv_set(other.m_pwd_rv_set)
{   // no need to do checks, they have been done when constructing other
    ; 
}


ngsai::KineticSignal::KineticSignal(ngsai::KineticSignal&& other)
    : m_size(std::move(other.m_size)),
      m_sequence_fw(std::move(other.m_sequence_fw)),
      m_sequence_rv(std::move(other.m_sequence_rv)),
      m_ipd_fw(std::move(other.m_ipd_fw)),
      m_ipd_rv(std::move(other.m_ipd_rv)),
      m_pwd_fw(std::move(other.m_pwd_fw)),
      m_pwd_rv(std::move(other.m_pwd_rv)),
      m_seq_fw_set(std::move(other.m_seq_fw_set)),
      m_seq_rv_set(std::move(other.m_seq_rv_set)),
      m_ipd_fw_set(std::move(other.m_ipd_fw_set)),
      m_ipd_rv_set(std::move(other.m_ipd_rv_set)),
      m_pwd_fw_set(std::move(other.m_pwd_fw_set)),
      m_pwd_rv_set(std::move(other.m_pwd_rv_set))
{   // no need to do checks, they have been done when constructing other
    ; 
}


ngsai::KineticSignal::~KineticSignal()
{ ; }


ngsai::KineticSignal&
ngsai::KineticSignal::operator=(const ngsai::KineticSignal& other)
{   m_size        = other.m_size ;
    m_sequence_fw = other.m_sequence_fw ;
    m_sequence_rv = other.m_sequence_rv ;
    m_ipd_fw      = other.m_ipd_fw ;
    m_ipd_rv      = other.m_ipd_rv ;
    m_pwd_fw      = other.m_pwd_fw ;
    m_pwd_rv      = other.m_pwd_rv ;
    m_seq_fw_set  = other.m_seq_fw_set ;
    m_seq_rv_set  = other.m_seq_rv_set ;
    m_ipd_fw_set  = other.m_ipd_fw_set ;
    m_ipd_rv_set  = other.m_ipd_rv_set ;
    m_pwd_fw_set  = other.m_pwd_fw_set ;
    m_pwd_rv_set  = other.m_pwd_rv_set ;
    // no need to do checks, they have been done when constructing other
    return *this ;
}


ngsai::KineticSignal&
ngsai::KineticSignal::operator=(ngsai::KineticSignal&& other)
{   m_size        = std::move(other.m_size) ;
    m_sequence_fw = std::move(other.m_sequence_fw) ;
    m_sequence_rv = std::move(other.m_sequence_rv) ;
    m_ipd_fw      = std::move(other.m_ipd_fw) ;
    m_ipd_rv      = std::move(other.m_ipd_rv) ;
    m_pwd_fw      = std::move(other.m_pwd_fw) ;
    m_pwd_rv      = std::move(other.m_pwd_rv) ;
    m_seq_fw_set  = std::move(other.m_seq_fw_set) ;
    m_seq_rv_set  = std::move(other.m_seq_rv_set) ;
    m_ipd_fw_set  = std::move(other.m_ipd_fw_set) ;
    m_ipd_rv_set  = std::move(other.m_ipd_rv_set) ;
    m_pwd_fw_set  = std::move(other.m_pwd_fw_set) ;
    m_pwd_rv_set  = std::move(other.m_pwd_rv_set) ;
    // no need to do checks, they have been done when constructing other
    return *this ;
}


void
ngsai::KineticSignal::setSequenceFw(const std::string& seq)
{
    if((m_size != 0) and (m_size != seq.size()))
    {   char msg[4096] ;
        sprintf(msg,
                "KineticSignal error : forward sequence size must be %zu (%zu)",
                m_size,
                seq.size()) ;
        throw std::invalid_argument(msg) ;
    }
    m_sequence_fw = seq ;
    m_seq_fw_set = true ;
    // in case m_size was 0 set it
    m_size = m_sequence_fw.size() ;
}


void
ngsai::KineticSignal::setSequenceFw(std::string&& seq)
{
    if((m_size != 0) and (m_size != seq.size()))
    {   char msg[4096] ;
        sprintf(msg,
                "KineticSignal error : forward sequence size must be %zu (%zu)",
                m_size,
                seq.size()) ;
        throw std::invalid_argument(msg) ;
    }
    m_sequence_fw = seq ;
    m_seq_fw_set = true ;
    // in case m_size was 0 set it
    m_size = m_sequence_fw.size() ;
}


void
ngsai::KineticSignal::setSequenceRv(const std::string& seq)
{
    if((m_size != 0) and (m_size != seq.size()))
    {   char msg[4096] ;
        sprintf(msg,
                "KineticSignal error : reverse sequence size must be %zu (%zu)",
                m_size,
                seq.size()) ;
        throw std::invalid_argument(msg) ;
    }
    m_sequence_rv = seq ;
    m_seq_rv_set = true ;
    // in case m_size was 0 set it
    m_size = m_sequence_rv.size() ;
}


void
ngsai::KineticSignal::setSequenceRv(std::string&& seq)
{
    if((m_size != 0) and (m_size != seq.size()))
    {   char msg[4096] ;
        sprintf(msg,
                "KineticSignal error : reverse sequence size must be %zu (%zu)",
                m_size,
                seq.size()) ;
        throw std::invalid_argument(msg) ;
    }
    m_sequence_rv = seq ;
    m_seq_rv_set = true ;
    // in case m_size was 0 set it
    m_size = m_sequence_rv.size() ;
}


void
ngsai::KineticSignal::setIPDFw(const std::vector<double>& ipd)
{
    if((m_size != 0) and (m_size != ipd.size()))
    {   char msg[4096] ;
        sprintf(msg,
                "KineticSignal error : forward IPD size must be %zu (%zu)",
                m_size,
                ipd.size()) ;
        throw std::invalid_argument(msg) ;
    }
    m_ipd_fw     = ipd ;
    m_ipd_fw_set = true ;
    // in case m_size was 0 set it
    m_size = ipd.size() ;
}


void
ngsai::KineticSignal::setIPDFw(std::vector<double>&& ipd)
{
    if((m_size != 0) and (m_size != ipd.size()))
    {   char msg[4096] ;
        sprintf(msg,
                "KineticSignal error : forward IPD size must be %zu (%zu)",
                m_size,
                ipd.size()) ;
        throw std::invalid_argument(msg) ;
    }
    m_ipd_fw     = ipd ;
    m_ipd_fw_set = true ;
    // in case m_size was 0 set it
    m_size = ipd.size() ;
}


void
ngsai::KineticSignal::setIPDRv(const std::vector<double>& ipd)
{
    if((m_size != 0) and (m_size != ipd.size()))
    {   char msg[4096] ;
        sprintf(msg,
                "KineticSignal error : reverse IPD size must be %zu (%zu)",
                m_size,
                ipd.size()) ;
        throw std::invalid_argument(msg) ;
    }
    m_ipd_rv     = ipd ;
    m_ipd_rv_set = true ;
    // in case m_size was 0 set it
    m_size = ipd.size() ;
}


void
ngsai::KineticSignal::setIPDRv(std::vector<double>&& ipd)
{
    if((m_size != 0) and (m_size != ipd.size()))
    {   char msg[4096] ;
        sprintf(msg,
                "KineticSignal error : reverse IPD size must be %zu (%zu)",
                m_size,
                ipd.size()) ;
        throw std::invalid_argument(msg) ;
    }
    m_ipd_rv     = ipd ;
    m_ipd_rv_set = true ;
    // in case m_size was 0 set it
    m_size = ipd.size() ;
}


void
ngsai::KineticSignal::setPWDFw(const std::vector<double>& pwd)
{
    if((m_size != 0) and (m_size != pwd.size()))
    {   char msg[4096] ;
        sprintf(msg,
                "KineticSignal error : forward PWD size must be %zu (%zu)",
                m_size,
                pwd.size()) ;
        throw std::invalid_argument(msg) ;
    }
    m_pwd_fw     = pwd ;
    m_pwd_fw_set = true ;
    // in case m_size was 0 set it
    m_size = pwd.size() ;
}


void
ngsai::KineticSignal::setPWDFw(std::vector<double>&& pwd)
{
    if((m_size != 0) and (m_size != pwd.size()))
    {   char msg[4096] ;
        sprintf(msg,
                "KineticSignal error : forward PWD size must be %zu (%zu)",
                m_size,
                pwd.size()) ;
        throw std::invalid_argument(msg) ;
    }
    m_pwd_fw     = pwd ;
    m_pwd_fw_set = true ;
    // in case m_size was 0 set it
    m_size = pwd.size() ;
}


void
ngsai::KineticSignal::setPWDRv(const std::vector<double>& pwd)
{
    if((m_size != 0) and (m_size != pwd.size()))
    {   char msg[4096] ;
        sprintf(msg,
                "KineticSignal error : reverse PWD size must be %zu (%zu)",
                m_size,
                pwd.size()) ;
        throw std::invalid_argument(msg) ;
    }
    m_pwd_rv     = pwd ;
    m_pwd_rv_set = true ;
    // in case m_size was 0 set it
    m_size = pwd.size() ;
}


void
ngsai::KineticSignal::setPWDRv(std::vector<double>&& pwd)
{
    if((m_size != 0) and (m_size != pwd.size()))
    {   char msg[4096] ;
        sprintf(msg,
                "KineticSignal error : reverse PWD size must be %zu (%zu)",
                m_size,
                pwd.size()) ;
        throw std::invalid_argument(msg) ;
    }
    m_pwd_rv     = pwd ;
    m_pwd_rv_set = true ;
    // in case m_size was 0 set it
    m_size = pwd.size() ;
}


const std::string&
ngsai::KineticSignal::getSequenceFw() const
{   return m_sequence_fw ; }


const std::string&
ngsai::KineticSignal::getSequenceRv() const
{   return m_sequence_rv ; }


const std::vector<double>&
ngsai::KineticSignal::getIPDFw() const
{   return m_ipd_fw ; }


const std::vector<double>&
ngsai::KineticSignal::getIPDRv() const
{   return m_ipd_rv ; }


const std::vector<double>&
ngsai::KineticSignal::getPWDFw() const
{   return m_pwd_fw ; }


const std::vector<double>&
ngsai::KineticSignal::getPWDRv() const
{   return m_pwd_rv ; }


std::vector<double>
ngsai::KineticSignal::getIPDMean() const
{   
    // no IPD to average
    if((not m_ipd_fw_set) and (not m_ipd_rv_set))
    {   return std::vector<double>() ; }

    std::vector<double> ipds_m(m_size, 0.) ;
    double n = 0. ;
    if(m_ipd_fw_set)
    {   for(size_t i=0; i<m_size; i++)
        {   ipds_m[i] += m_ipd_fw[i] ; }
        n += 1. ;
    }
    if(m_ipd_rv_set)
    {   for(size_t i=0; i<m_size; i++)
        {   ipds_m[i] += m_ipd_rv[i] ; }
        n += 1. ;
    }

    for(size_t i=0; i<m_size; i++)
    {   ipds_m[i] /= n ; }

    return ipds_m ;
}


std::vector<double>
ngsai::KineticSignal::getPWDMean() const
{   
    // no PWD to average
    if((not m_pwd_fw_set) and (not m_pwd_rv_set))
    {   return std::vector<double>() ; }
    
    std::vector<double> pwds_m(m_size, 0.) ;
    double n = 0. ;
    if(m_pwd_fw_set)
    {   for(size_t i=0; i<m_size; i++)
        {   pwds_m[i] += m_pwd_fw[i] ; }
        n += 1. ;
    }
    if(m_pwd_rv_set)
    {   for(size_t i=0; i<m_size; i++)
        {   pwds_m[i] += m_pwd_rv[i] ; }
        n += 1. ;
    }

    for(size_t i=0; i<m_size; i++)
    {   pwds_m[i] /= n ; }

    return pwds_m ;
}


size_t
ngsai::KineticSignal::size() const
{   return m_size ; }


bool
ngsai::KineticSignal::isFwComplete() const
{   return m_seq_fw_set and m_ipd_fw_set and m_pwd_fw_set ; }


bool
ngsai::KineticSignal::isRvComplete() const
{   return m_seq_rv_set and m_ipd_rv_set and m_pwd_rv_set ; }


bool
ngsai::KineticSignal::isComplete() const
{   return this->isFwComplete() and this->isRvComplete() ; }