#ifndef NGSAI_CCSKINETICEXTRACTOR_HPP
#define NGSAI_CCSKINETICEXTRACTOR_HPP

#include <vector>
#include <list>
#include <string>
#include <pbbam/BamRecord.h>

#include <ngsaipp/genome/GenomeRegion.hpp>



namespace ngsai {

    /*!
     * \brief The CcsKineticExtractor allows to extract inter pulse durations 
     * (IPDs) and pulse widths (PWDs) in given region from mapped PacBio CCSs. 
     * The region are always specified in reference coordinates (the reference 
     * onto which the CCS is mapped). The window specifies a region, on one 
     * of the reference strands, as a semi open interval [from, to). Then, if 
     * the CCS maps over this reference region, the CCS sequence, IPDs and 
     * PWDs values corresponding to the window are extracted. The sequence, 
     * IPDs and PWDs are always returned in 5' to 3' order with respect to the 
     * window strand (that is on the reference).
     */
    class CcsKineticExtractor
    {   
        public:
            
            /*!
             * \brief Constructs an object to extract kinetic values (IPDs 
             * or PWD) within windows from PacBio mapped CCSs.
             */
            CcsKineticExtractor() ;

            /*!
             * \brief Copy constructor.
             * \param window_size a reference to the object to copy.
             */
            CcsKineticExtractor(const ngsai::CcsKineticExtractor& other) ;

            /*!
             * \brief Move constructor.
             * \param window_size a rvalue to the object to move.
             */
            CcsKineticExtractor(ngsai::CcsKineticExtractor&& other) ;
            
            /*!
             * \brief Destructor
             */
            ~CcsKineticExtractor() ;

            /*!
             * \brief Assignment operator.
             * \param other A reference to the object to assign the value 
             * from.
             * \returns A reference to the current object.
             */
            ngsai::CcsKineticExtractor& 
            operator = (const ngsai::CcsKineticExtractor& other) ;

            /*!
             * \brief Move assignment operator.
             * \param other A rvalue to the object to move the value 
             * from.
             * \returns A reference to the current object.
             */
            ngsai::CcsKineticExtractor& 
            operator = (ngsai::CcsKineticExtractor&& other) ;

            /*!
             * \brief Extracts the kinetic data from the given mapped CCS 
             * over the window of interest. If the CCS is not mapped, does 
             * not contain kinetics, does not span entirely the window or if 
             * the alignment over the window is not perfect (only matches 
             * allowed), nothing is extracted.
             * \param ccs A mapped CCS to extract the kinetic data from.
             * \param window The coordinates [from,to) of the window of 
             * interest, expressed in reference forward coordinates.
             * \returns whether something could be extracted or not.
             */ 
            bool 
            extract(const PacBio::BAM::BamRecord& ccs,
                    const ngsai::genome::GenomeRegion& window) ;

            /*!
             * \brief Returns the sequence in 5' to 3' 
             * order. This is only meaningful after a successful extract() 
             * call.
             * \returns The DNA sequence corresponding to the window 
             * in which the extraction was performed. The sequence is 
             * always returned in 5' to 3' with respect to the region it was 
             * taken from. For instance, if the region was on the + strand 
             * of the reference (reference forward), the 0th base of the 
             * sequence represents the 5'most base for the + strand. 
             * Similarily, if the region was on the - strand of the reference 
             * (reference reverse), the 0th base of the sequence represents 
             * the 5'most base for the - strand.
             */
            std::string
            getSequence() const ;

            /*!
             * \brief Returns the decoded extracted IPD values in 5' to 3' 
             * order. This is only meaningful after a successful extract() 
             * call.
             * \returns A vector of IPD values corresponding to the window 
             * in which the extraction was performed. The data IPDs are 
             * always returned in 5' to 3' with respect to the region it was 
             * taken from. For instance, if the region was on the + strand 
             * of the reference (reference forward), the 0th element of the 
             * IPDs represents the 5'most value for the + strand. Similarily, 
             * if the region was on the - strand of the reference (reference 
             * reverse), the 0th element of the IPDs represents the 5'most 
             * value for the - strand.
             */
            std::vector<uint16_t> 
            getIPD() const ;

            /*!
             * \brief Returns the decoded extracted PWD values in 5' to 3' 
             * order. This is only meaningful after a successful extract() 
             * call.
             * \returns A vector of PWD values corresponding to the window 
             * in which the extraction was performed. The data PWDs are 
             * always returned in 5' to 3' with respect to the region it was 
             * taken from. For instance, if the region was on the + strand 
             * of the reference (reference forward), the 0th element of the 
             * PWDs represents the 5'most value for the + strand. Similarily, 
             * if the region was on the - strand of the reference (reference 
             * reverse), the 0th element of the PWDs represents the 5'most 
             * value for the - strand.
             */
            std::vector<uint16_t> 
            getPWD() const ;
        

        protected:
            
            /*!
             * \brief Parses the cigar of the CCS to find the given region 
             * - in reference coordinates - in the read. The range [from,to) 
             * corresponding to the region - in read coordinates - is stored 
             * in m_region_start_ref and m_region_end_ref.
             * A match with the given region will only be reported if the 
             * alignment of the read over the region is perfect (only matches
             * allowed). 
             * \param ccs the CCS of interest.
             * \param window the window of interest.
             * \returns whether the window coordinates could be found in the 
             * CCS.
             */
            bool findWindow(const PacBio::BAM::BamRecord& ccs,
                            const ngsai::genome::GenomeRegion& window) ;

            /*!
             * \brief Parses the CIGAR of the CCS and positions the CCS 
             * pointer on the position corresponding to the 1st position 
             * aligned on the reference (skips initial soft clips).
             * \returns whether the 1st position aligned on the reference 
             * could be found.
             */
            bool findAlignmentStart() ;

            /*!
             * \brief Parses the CIGAR of the CCS and positions the CCS 
             * pointer on the position of the CCS that is aligned with the 
             * 1st position of the reference window (it can be a match or a 
             * missmatch).
             * \returns Whether a position in the CCS aligned with the 1st 
             * position in the window could be found.
             */
            bool findWindowStart() ;

            /*!
             * \brief Parses the CIGAR of the CCS until finding the position 
             * of the CCS that is aligned to the last position in the 
             * reference window and position the CCS pointer on the position 
             * right after. All positions in the alignment must be matches, 
             * any indel or missmatch is forbidden. 
             * \returns Whether a perfect match along the entire reference 
             * window could be found.
             */
            bool findWindowEnd() ;
        private:
            /*!
             * \brief The window sequence, from 5' to 3', with respect to the 
             * window strand.
             */
            std::string m_window_seq ;
            /*!
             * \brief The window IPDs, from 5' to 3', with respect to the 
             * window strand.
             */
            std::vector<uint16_t>  m_window_ipd ;
            /*!
             * \brief The window PWDs, from 5' to 3', with respect to the 
             * window strand.
             */
            std::vector<uint16_t>  m_window_pwd ;



            /*!
             * \brief A list containing the CIGAR string operations of the 
             * CCS of interest.
             */ 
            std::list<PacBio::BAM::CigarOperation> m_cigar ;
            /*!
             * \brief The index of the current CIGAR operation to check in the
             * list. Also correspond to the position - in FORWARD orientation 
             * (as the CIGAR) - in the read.
             */ 
            size_t m_read_i ;
            /*!
             * \brief The position in the reference corresponding to the 
             * current CIGAR operation.
             */ 
            size_t m_ref_i ;
            /*!
             * \brief The coordinate of the region start on the reference.
             */
            size_t m_region_start_ref ;
            /*!
             * \brief The coordinate of the region end on the reference.
             */
            size_t m_region_end_ref ;
            /*!
             * \brief The coordinate of the region start on the CCS.
             */
            size_t m_region_start_ccs ;
            /*!
             * \brief The coordinate of the region end on the CCS.
             */
            size_t m_region_end_ccs ;

    } ;
}

#endif // NGSAI_CCSKINETICEXTRACTOR_HPP