#include <ngsaipp/epigenetics/CcsKineticExtractor.hpp>

#include <vector>
#include <list>
#include <string>
#include <pbbam/BamRecord.h>
#include <limits>
#include <iostream>

#include <ngsaipp/genome/GenomeRegion.hpp>  // ngsai::genome::GenomeRegion
#include <ngsaipp/genome/constants.hpp>     // ngsai::genome::strand
#include <ngsaipp/dna/dna_utility.hpp>      // ngsai::dna::get_reverse_complement()

void print_cigar_expanded_(std::ostream& stream,
                           const std::vector<PacBio::BAM::CigarOperation>& cigar)
{   for(const auto& op : cigar)
    {   stream << op.Char() ; }
}

void print_cigar_expanded_(std::ostream& stream,
                           const std::list<PacBio::BAM::CigarOperation>& cigar)
{   for(const auto& op : cigar)
    {   stream << op.Length() << op.Char() ; }
}



ngsai::CcsKineticExtractor::CcsKineticExtractor()
    : m_window_seq(),
      m_window_ipd(),
      m_window_pwd(),
      m_cigar(),
      m_read_i(),
      m_ref_i(),
      m_region_start_ref(),
      m_region_end_ref(),
      m_region_start_ccs(),
      m_region_end_ccs()
{ ; }


ngsai::CcsKineticExtractor::CcsKineticExtractor(
                        const ngsai::CcsKineticExtractor& other)
    : m_window_seq(other.m_window_seq),
      m_window_ipd(other.m_window_pwd),
      m_window_pwd(other.m_window_pwd),
      m_cigar(other.m_cigar),
      m_read_i(other.m_read_i),
      m_ref_i(other.m_ref_i),
      m_region_start_ref(other.m_region_start_ref),
      m_region_end_ref(other.m_region_end_ref),
      m_region_start_ccs(other.m_region_start_ccs),
      m_region_end_ccs(other.m_region_end_ccs)
{ ; }


ngsai::CcsKineticExtractor::CcsKineticExtractor(
                        ngsai::CcsKineticExtractor&& other)
    : m_window_seq(std::move(other.m_window_seq)),
      m_window_ipd(std::move(other.m_window_pwd)),
      m_window_pwd(std::move(other.m_window_pwd)),
      m_cigar(std::move(other.m_cigar)),
      m_read_i(std::move(other.m_read_i)),
      m_ref_i(std::move(other.m_ref_i)),
      m_region_start_ref(std::move(other.m_region_start_ref)),
      m_region_end_ref(std::move(other.m_region_end_ref)),
      m_region_start_ccs(std::move(other.m_region_start_ccs)),
      m_region_end_ccs(std::move(other.m_region_end_ccs))
{ ; }


ngsai::CcsKineticExtractor::~CcsKineticExtractor()
{ ; }


ngsai::CcsKineticExtractor& 
ngsai::CcsKineticExtractor::operator = (
                        const ngsai::CcsKineticExtractor& other)
{
    m_window_seq       = other.m_window_seq ;
    m_window_ipd       = other.m_window_ipd ;
    m_window_pwd       = other.m_window_pwd ;
    m_cigar            = other.m_cigar ;
    m_read_i           = other.m_read_i ;
    m_ref_i            = other.m_ref_i ;
    m_region_start_ref = other.m_region_start_ref ;
    m_region_end_ref   = other.m_region_end_ref ;
    m_region_start_ccs = other.m_region_start_ccs ;
    m_region_end_ccs   = other.m_region_end_ccs ;
    return *this ;
}


ngsai::CcsKineticExtractor& 
ngsai::CcsKineticExtractor::operator = (
                        ngsai::CcsKineticExtractor&& other)
{
    m_window_seq        = std::move(other.m_window_seq) ;
    m_window_ipd        = std::move(other.m_window_ipd) ;
    m_window_pwd        = std::move(other.m_window_pwd) ;
    m_cigar             = std::move(other.m_cigar) ;
    m_read_i            = std::move(other.m_read_i) ;
    m_ref_i             = std::move(other.m_ref_i) ;
    m_region_start_ref  = std::move(other.m_region_start_ref) ;
    m_region_end_ref    = std::move(other.m_region_end_ref) ;
    m_region_start_ccs  = std::move(other.m_region_start_ccs) ;
    m_region_end_ccs    = std::move(other.m_region_end_ccs) ;
    return *this ;
}


bool 
ngsai::CcsKineticExtractor::extract(const PacBio::BAM::BamRecord& ccs,
                                    const ngsai::genome::GenomeRegion& window)
{   
    // reset fields of interest
    m_read_i = 0 ;  
    m_ref_i = 0 ;
    m_region_start_ref = 0 ;
    m_region_end_ref = 0 ;
    m_region_start_ccs = 0 ;
    m_region_end_ccs = 0 ;
    m_window_ipd.clear() ;
    m_window_pwd.clear() ;
    m_window_seq.clear() ;
    // m_window_cigar.clear() ;

    // obviously, the window cannot be found if not mapped
    if(not ccs.IsMapped())
    {   return false ; }

    // extract CIGAR string
    PacBio::BAM::Cigar cigar_v = ccs.CigarData() ;
    m_cigar = std::list<PacBio::BAM::CigarOperation>(cigar_v.begin(),
                                                     cigar_v.end()) ;

    // find region in read in REF FW coordinates
    if(not this->findWindow(ccs, window))
    {   return false ; }
    

    // strand on which window and CCS are
    bool is_region_fw = window.strand == ngsai::genome::FORWARD ;
    bool is_read_fw   = ccs.AlignedStrand() == 
                          PacBio::BAM::Strand::FORWARD ;

    // read sequence/IPD are given in NATIVE orientation, not 
    // in GENOMIC (aka REFERENCE FORWARD) orientation.
    // read_window_coord are given in REFERENCE FORWARD orientation
    // if read maps on REFERENCE REVERSE, read_window_coord need to 
    // be translated
    //    map in REFERENCE fw -> all ok
    //             from      to
    //    ----------|---------|---------------------> REF fw
    //        ------|---------|------>                CIGAR    
    //        ------|---------|------>                Sequence
    //        0 1 2 3 4 5 6 7 8 9 10
    //
    //    map in REFERENCE rv -> translation
    //             from      to
    //    ----------|---------|---------------------> REF fw
    //        0 1 2 3 4 5 6 7 8 9 10
    //        ------|---------|------>                CIGAR    
    //       <------|---------|------                 Sequence
    //              to'     from'
    //         10 9 8 7 6 5 4 3 2 1 0
    
    // extract sequence
    // GENOMIC orientation means as in BAM (always in REF fw orientation)
    std::string seq = 
        ccs.Sequence(PacBio::BAM::Orientation::GENOMIC) ;
    // slice
    m_window_seq = 
            std::string(seq.begin() + m_region_start_ccs,
                        seq.begin() + m_region_end_ccs) ;

    // length of read
    size_t read_len = seq.size() ;

    // extracts IPD and PW
    // NATIVE orientation means as it was for unaligned CCS (oriented 
    // as subread SEQ was)
    // normally alignment do NOT alter the direction of IPD
    // and PWD in BAM. Thus NATIVE returns IPD/PWD as they were for unaligned 
    // CCS.
    if(is_region_fw)
    {   // CCS maps in fw
        if(is_read_fw)
        {   // IPD are decoded
            std::vector<uint16_t> ipd = 
                ccs.ForwardIPD(PacBio::BAM::Orientation::NATIVE).Data() ;
            // CCS with too few passes have no IPD
            if(ipd.size() == 0)
            {   return false ; }
            m_window_ipd = std::vector<uint16_t>(
                            ipd.begin() + m_region_start_ccs,
                            ipd.begin() + m_region_end_ccs) ;
            // Pulse width
            std::vector<uint16_t> pwd = 
                ccs.ForwardPulseWidth(PacBio::BAM::Orientation::NATIVE).Data() ;
            m_window_pwd = std::vector<uint16_t>(
                            pwd.begin() + m_region_start_ccs,
                            pwd.begin() + m_region_end_ccs) ;
        }
        // CCS maps in rv
        else
        {   // IPD are decoded
            std::vector<uint16_t> ipd = 
                ccs.ReverseIPD(PacBio::BAM::Orientation::NATIVE).Data() ;
            // CCS with too few passes have no IPD
            if(ipd.size() == 0)
            {   return false ; }
            // slice
            m_window_ipd = std::vector<uint16_t>(
                            ipd.begin() + m_region_start_ccs,
                            ipd.begin() + m_region_end_ccs) ;
            
            // Pulse width
            std::vector<uint16_t> pwd = 
                ccs.ReversePulseWidth(PacBio::BAM::Orientation::NATIVE).Data() ;
            // slice
            m_window_pwd = std::vector<uint16_t>(
                            pwd.begin() + m_region_start_ccs,
                            pwd.begin() + m_region_end_ccs) ;
        }
    }
    // window is on rev strand
    else
    {   // CCS maps in fw
        if(is_read_fw)
        {   
            // IPD are decoded
            std::vector<uint16_t> ipd = 
                ccs.ReverseIPD(PacBio::BAM::Orientation::NATIVE).Data() ;
            // CCS with too few passes have no IPD
            if(ipd.size() == 0)
            {   return false ; }
             // convert [start,end) on fw strand coord to reverse strand coord
            size_t to_   = read_len - m_region_start_ccs ;
            size_t from_ = read_len - m_region_end_ccs ;
            m_region_start_ccs = from_ ;
            m_region_end_ccs = to_ ;
            // slice
            m_window_ipd = std::vector<uint16_t>(
                        ipd.begin() + m_region_start_ccs,
                        ipd.begin() + m_region_end_ccs) ;

            // Pulse width
            std::vector<uint16_t> pwd = 
                ccs.ReversePulseWidth(PacBio::BAM::Orientation::NATIVE).Data() ;
            // slice
            m_window_pwd = std::vector<uint16_t>(
                        pwd.begin() + m_region_start_ccs,
                        pwd.begin() + m_region_end_ccs) ;

            // reverse complement of sequence to fit window orientation (rv)
            m_window_seq = ngsai::dna::get_reverse_complement(m_window_seq) ;
        }
        // read maps in rv
        else
        {   // IPD are decoded
            std::vector<uint16_t> ipd = 
                ccs.ForwardIPD(PacBio::BAM::Orientation::NATIVE).Data() ;
            // CCS with too few passes have no IPD
            if(ipd.size() == 0)
            {   return false ; }
            // convert [start,end) on fw strand coord to reverse strand coord
            size_t to_   = read_len - m_region_start_ccs ;
            size_t from_ = read_len - m_region_end_ccs ;
            m_region_start_ccs = from_ ;
            m_region_end_ccs = to_ ;

            // slice
            m_window_ipd = std::vector<uint16_t>(
                        ipd.begin() + m_region_start_ccs,
                        ipd.begin() + m_region_end_ccs) ;
            
            // Pulse width
            std::vector<uint16_t> pwd = 
                ccs.ForwardPulseWidth(PacBio::BAM::Orientation::NATIVE).Data() ;
            // slice
            m_window_pwd = std::vector<uint16_t>(
                        pwd.begin() + m_region_start_ccs,
                        pwd.begin() + m_region_end_ccs) ;

            // reverse complement of sequence to fit window orientation (rv)
            m_window_seq = ngsai::dna::get_reverse_complement(m_window_seq) ;
        }
    }

    return true ;
}


std::string
ngsai::CcsKineticExtractor::getSequence() const
{   return m_window_seq ; }


std::vector<uint16_t> 
ngsai::CcsKineticExtractor::getIPD() const
{   return m_window_ipd ; }


std::vector<uint16_t> 
ngsai::CcsKineticExtractor::getPWD() const
{   return m_window_pwd ; }


bool
ngsai::CcsKineticExtractor::findWindow(
                                    const PacBio::BAM::BamRecord& ccs,
                                    const ngsai::genome::GenomeRegion& window)
{   
    // check region and read are on same chromosome
    if(ccs.ReferenceName() != window.chrom)
    {   return false ; }

    // pos where read map on ref
    size_t ref_start_align = ccs.ReferenceStart() ;
    size_t ref_end_align   = ccs.ReferenceEnd() ;

    // pos region on ref
    m_region_start_ref  = window.start ;
    m_region_end_ref    = window.end ;

    // check that read spans the entire region
    // -> start and end can be found in the read
    if((ref_start_align > m_region_start_ref) or
       (ref_end_align   < m_region_end_ref))
    {   return false ; }

    // pointer to ref sequence
    m_ref_i  = ref_start_align ;
    // pointer to read sequence
    m_read_i = 0 ;

    // parse alignment until 1st postion aligned to ref
    if(not this->findAlignmentStart())
    {   return false ; }
    
    // parse alignment until region starts
    if(not this->findWindowStart())
    {   return false ; }
    m_region_start_ccs = m_read_i ;

    // parse alignment until region ends
    if(not this->findWindowEnd())
    {   return false ; }
    m_region_end_ccs = m_read_i ;

    return true ;
}


bool
ngsai::CcsKineticExtractor::findAlignmentStart()
{   bool found = false ;
    while(m_cigar.size())
    {   auto cigar_element = m_cigar.front() ;
        auto op = cigar_element.Type() ;
        auto len = cigar_element.Length() ;
        switch(op)
        {   // alignment did not start yet
            case PacBio::Data::CigarOperationType::SOFT_CLIP:
                m_read_i += len ;
                m_cigar.pop_front() ;
                break ;
            // alignment starts here
            default:
                found = true ;
                break ;
        }
        if(found)
        {   break ; }
    }
    return found ;
}


bool
ngsai::CcsKineticExtractor::findWindowStart()
{   
    bool found = false ;  // end position has been found
    bool done  = false ;  // must terminate, does not mean found is true

    size_t read_i_new = m_read_i ;
    size_t ref_i_new  = m_ref_i ;

    while(m_cigar.size())
    {
        // the stretch of alignment positions to check
        PacBio::BAM::CigarOperation op(m_cigar.front()) ;
        m_cigar.pop_front() ;
        size_t op_len = op.Length() ;  // 6M -> 6
        
        // move pointers after entire operation
        switch (op.Type())
        {   // consums only query
            case PacBio::Data::CigarOperationType::SOFT_CLIP:
                // this->read_i++ ;
                // hit the unaligned end of the read, cannot find start
                // std::cerr << "findWindowStartNew S" << std::endl ;
                done = true ;
                break ;
            case PacBio::Data::CigarOperationType::INSERTION:
                // std::cerr << "findWindowStartNew I" << std::endl ;
                read_i_new += op_len ;
                break ;
            // consums only reference
            case PacBio::Data::CigarOperationType::DELETION:
                // std::cerr << "findWindowStartNew D" << std::endl ;
                ref_i_new += op_len ;
                break ;
            case PacBio::Data::CigarOperationType::REFERENCE_SKIP:
                // std::cerr << "findWindowStartNew N" << std::endl ;
                ref_i_new += op_len ;
                break ;
            // consums both
            case PacBio::Data::CigarOperationType::SEQUENCE_MATCH:
                // std::cerr << "findWindowStartNew =" << std::endl ;
                ref_i_new  += op_len ;
                read_i_new += op_len ;
                break ;
            case PacBio::Data::CigarOperationType::SEQUENCE_MISMATCH:
                // std::cerr << "findWindowStartNew X" << std::endl ;
                ref_i_new  += op_len ;
                read_i_new += op_len ;
                break ;
            case PacBio::Data::CigarOperationType::ALIGNMENT_MATCH:
                // std::cerr << "findWindowStartNew M" << std::endl ;
                ref_i_new += op_len ;
                read_i_new += op_len ;
                break ;
            default:
                char msg[4096] ;
                sprintf(msg, "unexpected cigar operation (%d%c)",
                        op.Length(), op.Char()) ;
                throw std::runtime_error(msg) ;
        }

        // hit soft clip
        if(done)
        {   break ; }
        // window start is right here
        if(ref_i_new == m_region_start_ref)
        {   m_read_i = read_i_new ;
            m_ref_i  = ref_i_new ;
            done =  true ;
            found = true ;
            break ;
        }
        // we went too far
        // if so, window start has been found, move pointers
        // back there and push back remaining of operation
        // in cigar
        if(ref_i_new > m_region_start_ref)
        {   
            size_t extra = ref_i_new - m_region_start_ref ;
            m_read_i = read_i_new - extra ;
            m_ref_i  = ref_i_new  - extra ;
            op.Length(extra) ;
            m_cigar.push_front(op) ;
            done =  true ;
            found = true ;
            break ;
        }
    }
    return found ;
}


bool
ngsai::CcsKineticExtractor::findWindowEnd()
{   
    bool found = false ;  // end position has been found
    bool done  = false ;  // must terminate, does not mean found is true

    size_t read_i_new = m_read_i ;
    size_t ref_i_new  = m_ref_i ;

    while(m_cigar.size())
    {   
        // the stretch of alignment positions to check
        PacBio::BAM::CigarOperation op(m_cigar.front()) ;
        m_cigar.pop_front() ;
        size_t op_len = op.Length() ;  // 6M -> 6

        // only matches are allowed in window
        if(op.Type() != PacBio::Data::CigarOperationType::SEQUENCE_MATCH and
           op.Type() != PacBio::Data::CigarOperationType::ALIGNMENT_MATCH)
        {   found = false ;
            break  ;
        }

        // move pointers after entire operation
        switch (op.Type())
        {   // consums only query
            case PacBio::Data::CigarOperationType::SOFT_CLIP:
                // this->read_i++ ;
                // hit the unaligned end of the read, cannot find start
                // std::cerr << "findWindowStartNew S" << std::endl ;
                done = true ;
                break ;
            case PacBio::Data::CigarOperationType::INSERTION:
                // std::cerr << "findWindowStartNew I" << std::endl ;
                read_i_new += op_len ;
                break ;
            // consums only reference
            case PacBio::Data::CigarOperationType::DELETION:
                // std::cerr << "findWindowStartNew D" << std::endl ;
                ref_i_new += op_len ;
                break ;
            case PacBio::Data::CigarOperationType::REFERENCE_SKIP:
                // std::cerr << "findWindowStartNew N" << std::endl ;
                ref_i_new += op_len ;
                break ;
            // consums both
            case PacBio::Data::CigarOperationType::SEQUENCE_MATCH:
                // std::cerr << "findWindowStartNew =" << std::endl ;
                ref_i_new += op_len ;
                read_i_new += op_len ;
                break ;
            case PacBio::Data::CigarOperationType::SEQUENCE_MISMATCH:
                // std::cerr << "findWindowStartNew X" << std::endl ;
                ref_i_new += op_len ;
                read_i_new += op_len ;
                break ;
            case PacBio::Data::CigarOperationType::ALIGNMENT_MATCH:
                // std::cerr << "findWindowStartNew M" << std::endl ;
                ref_i_new += op_len ;
                read_i_new += op_len ;
                break ;
            default:
                char msg[4096] ;
                sprintf(msg, "unexpected cigar operation (%d%c)",
                        op.Length(), op.Char()) ;
                throw std::runtime_error(msg) ;
        }

        // hit soft clip
        if(done)
        {   break ; }
        // window start is right here
        if(ref_i_new == m_region_end_ref)
        {   m_read_i = read_i_new ;
            m_ref_i  = ref_i_new ;
            done =  true ;
            found = true ;
            break ;
        }
        // we went too far
        // if so, window start has been found, move pointers
        // back there and push back remaining of operation
        // in cigar
        if(ref_i_new > m_region_end_ref)
        {   
            size_t extra = ref_i_new - m_region_end_ref ;
            m_read_i = read_i_new - extra ;
            m_ref_i  = ref_i_new  - extra ;
            op.Length(extra) ;
            m_cigar.push_front(op) ;
            done =  true ;
            found = true ;
            break ;
        }
    }

    return found ;
}
