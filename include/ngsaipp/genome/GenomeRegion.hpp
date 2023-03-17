#ifndef NGSAI_GENOME_GENOMICREGION_HPP
#define NGSAI_GENOME_GENOMICREGION_HPP

#include <string>
#include <ngsaipp/genome/constants.hpp>    // ngsai::genome::strand
#include <ngsaipp/utility/constants.hpp>   // ngsai::bool_ext


namespace ngsai
{   
    namespace genome
    {
        /*!
        * A class modeling a region on a genome. A region is defined as a 
        * segment of contigous bases, [start, end), defined by 0-based 
        * coordinates start and end.
        */
        class GenomeRegion
        {
            public:
                /*!
                 * \brief Constructs an empty object.
                 */
                GenomeRegion() ;

                /*!
                 * \brief Copy constructor.
                 * \param other the instance to copy from.
                 * \throws std::invalid_argument if something wrong is 
                 * detected.
                 */
                GenomeRegion(const ngsai::genome::GenomeRegion& other) ;

                /*!
                 * \brief Move constructor.
                 * \param other the instance to move from.
                 * \throws std::invalid_argument if something wrong is 
                 * detected.
                 */
                GenomeRegion(ngsai::genome::GenomeRegion&& other) ;

                /*!
                 * \brief Constructs an object with the given values.
                 * \param chrom the chromosome name.
                 * \param start the 0-based position of the 1st base inside 
                 * the range.
                 * \param end the 0-based  position of the past last base in 
                 * the range. The range is thus [start, end).
                 * \param strand the strand on which the region is located.
                 * \throws std::invalid_argument if something wrong is 
                 * detected.
                 */
                GenomeRegion(const std::string& chrom,
                             size_t start,
                             size_t end,
                             ngsai::genome::strand strand) ;

                /*!
                 * \brief Destructor.
                 */
                virtual
                ~GenomeRegion() ;

                /*!
                 * \brief Assigmnent operator.
                 */
                GenomeRegion&
                operator = (const GenomeRegion& other) ;

                /*!
                 * \brief Move assigmnent operator.
                 */
                GenomeRegion&
                operator = (GenomeRegion&& other) ;

                /*!
                 * \brief Compares if two GenomeRegion have equal values.
                 * \param other the BedRecord to compare the current instance 
                 * with.
                 * \returns True if all fields values are the same.
                 */
                bool
                operator == (const GenomeRegion& other) const ;

                /*!
                 * \brief Compares if two BedRecord have unequal values.
                 * \param other the BedRecord to compare the current instance 
                 * to.
                 * \returns True if at least one fields value is different.
                 */
                bool
                operator != (const GenomeRegion& other) const ;


                /*!
                 * \brief Overlap operator.
                 * Returns whether both regions overlap by at least one 
                 * position.
                 * \param other the other region.
                 * \return whether there is an overlap.
                 */
                bool
                operator | (const GenomeRegion& other) const ;

                /*!
                 * \brief Is upstream operator. Checks if the current object 
                 * islocated stricly upstream (no overlap) of the other region.
                 * The chromosome values are also tested. A smaller chromosome 
                 * value (natural order) is interpreted as upstream.
                 * \param other the other region.
                 * \return if the current region is upstream, without 
                 * overlapping, the other region.
                 */
                ngsai::bool_ext
                operator < (const GenomeRegion& other) const ;

                /*!
                 * \brief Is downstream operator. Checks if the current object 
                 * is located stricly downstream (no overlap) of the other 
                 * region. The chromosome values are also tested. A smaller 
                 * chromosome value (natural order) is interpreted as upstream.
                 * \param other the other region.
                 * \return if the current region is downstream, without 
                 * overlapping, the other region.
                 */
                ngsai::bool_ext
                operator > (const GenomeRegion& other) const ;

                /*!
                 * \brief Returns the size (length) of the region in bp.
                 * \returns the size of the region.
                 */
                virtual size_t
                size() const ;

                /*!
                 * \brief Returns the length of the overlap
                 * between both regions.
                 * \param other the second region.
                 * \return the length of the overlap
                 * between both regions.
                 */
                virtual int
                overlap_len(const GenomeRegion& other) const ;

                /*!
                 * \brief Creates a string representation of the instance.
                 * \returns the string representation of the record.
                 */
                virtual std::string
                toString() const ;

            
            private:
                /*!
                 * \brief Checks that start < end.
                 * \throws std::invalid_argument if something wrong is 
                 * detected.
                 */
                virtual void
                sanity_check() const ;

            public:
                /*
                 * \brief the chromosome name.
                 */
                std::string chrom ;
                /*!
                 * \brief The starting position of the feature in the 
                 * chromosome or scaffold. The first base in a chromosome is 
                 * numbered 0.
                 */
                size_t start ;
                /*!
                * \brief The ending position of the feature in the chromosome 
                * or scaffold. The end base is not included, the range is 
                * [start, end)
                */
                size_t end ;
                /*!
                * \brief The strand on which the bed interval is located.
                */
                ngsai::genome::strand strand ;
        } ;


        /*!
        * \brief Writes the given region to a stream.
        * \param stream the stream to write the record to.
        * \param record the record to write.
        * \returns a reference to the stream after the writing.
        */
        std::ostream& operator << (std::ostream& stream,
                                   const ngsai::genome::GenomeRegion& record) ;
    } // namespace genome

} // namespace ngsai


#endif // NGSAI_GENOME_GENOMICREGION_HPP