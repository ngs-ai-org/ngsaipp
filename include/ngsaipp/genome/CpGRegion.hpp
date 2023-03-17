#ifndef NGSAI_GENOME_CPG_HPP
#define NGSAI_GENOME_CPG_HPP

#include <ngsaipp/genome/GenomeRegion.hpp>


namespace ngsai
{
    namespace genome
    {   
        /*!
         * \brief The CpGRegion class models the coordinates of a CpG on the 
         * genome. 
         * CpG are palyndromic sequences. Hence, if a CpG is located on a 
         * strand, another is located on the other strand (but in reverse 
         * orientation). For this reason, the strand information is always 
         * UNORIENTED.
         * A CpG is encoded as a two GenomicRegion containing [start,end) 
         * ranges (0-based) for each strand CpG. The starts are always the 
         * position of the leftmost C on both strands and the end positions 
         * the position just after the rightmost G on both strands. This 
         * corresponds to the following schema :
         *          start end
         *            |   |
         * 5' - N N N C G N N N - 3' fw strand
         * 3' - N N N G C N N N - 5' rv strand
         *            |   |
         *          start end
         */          
        class CpGRegion : public ngsai::genome::GenomeRegion
        {   
            public:
                /*!
                 * \brief Constructs and empy instance.
                 */
                CpGRegion() ;

                /*!
                 * \brief Copy constructor.
                 * \param cpg the region to copy from.
                 * \throws std::invalid_argument if something wrong is 
                 * detected.
                 */
                CpGRegion(const ngsai::genome::GenomeRegion& cpg) ;

                /*!
                 * \brief Move constructor.
                 * \param cpg the region to move from.
                 * \throws std::invalid_argument if something wrong is 
                 * detected.
                 */
                CpGRegion(ngsai::genome::GenomeRegion&& cpg) ;

                /*!
                 * \brief Constructs a CpG corresponding to the given genomic 
                 * coordinates.
                 * The strand is forced to UNORIENTED.
                 * \param chrom the name of the chromosome on which the CpG is 
                 * located.
                 * \param start the 0-based coordinate of the leftmost C.
                 * \param end the 0-based coordinate of the position just 
                 * after the rightmost G. 
                 * \throws std::invalid_argument if something wrong is 
                 * detected.
                 */
                CpGRegion(const std::string& chrom,
                          size_t start,
                          size_t end) ;
                /*!
                 * \brief Destructor.
                 */
                virtual
                ~CpGRegion() ;

                /*!
                 * \brief Returns the coordinates of the forward strand 
                 * CpG.
                 * \returns the forward strand CpG coordinates.
                 */
                ngsai::genome::GenomeRegion
                forwardCoordinates() const ;

                /*!
                 * \brief Returns the coordinates of the reverse strand 
                 * CpG.
                 * \returns the forward strand CpG coordinates.
                 */
                ngsai::genome::GenomeRegion
                reverseCoordinates() const ;

            protected:
                /*!
                 * \brief Checks that it has a length of 2bp and that the 
                 * strand is UNORIENTED.
                 * \throws std::invalid_argument if something wrong is 
                 * detected.
                 */
                void 
                sanity_check() const ; 
            
        } ;


    }  // namespace genome

}  // namespace ngsai


#endif // NGSAI_GENOME_CPG_HPP