#ifndef NGSAI_IO_BEDRECORD_HPP
#define NGSAI_IO_BEDRECORD_HPP

#include <string>
#include <ngsaipp/genome/GenomeRegion.hpp>  // ngsai::genome::GenomeRegion
#include <ngsaipp/genome/constants.hpp>     // ngsai::genome::strand

namespace ngsai
{   
    /*!
    * A class to handle a record in a bed file. It can handle up to bed 6
    * format.
    */
    class BedRecord : public ngsai::genome::GenomeRegion
    {
        public:
            /*!
             * \brief Constructs an empty object.
             */
            BedRecord() ;

            /*!
             * \brief Copy constructor.
             * \param other the BedRecord to copy from.
             * \throws std::invalid_argument if something wrong is detected.
             */
            BedRecord(const BedRecord& other) ;

            /*!
             * \brief Move constructor.
             * \param other the BedRecord to move from.
             * \throws std::invalid_argument if something wrong is detected.
             */
            BedRecord(BedRecord&& other) ;

            /*!
             * \brief Constructs an object with the given values.
             * \param chrom the chromosome name.
             * \param start the 0-based position of the 1st base inside the 
             * range.
             * \param end the 0-based  position of the past last base in the 
             * range. The range is thus [start, end).
             * \param name the name of the region.
             * \param score the score of the region.
             * \param strand the strand on which the region is located.
             * \throws std::invalid_argument if something wrong is detected.
             */
            BedRecord(const std::string& chrom,
                      size_t start,
                      size_t end,
                      const std::string& name,
                      double score,
                      ngsai::genome::strand strand) ;

            /*!
             * \brief Destructor.
             */
            ~BedRecord() ;

            /*!
             * \brief Assigmnent operator.
             */
            BedRecord& operator = (const BedRecord& other) ;

            /*!
             * \brief Move assigmnent operator.
             */
            BedRecord& operator = (BedRecord&& other) ;

            /*!
             * \brief Compares if two BedRecord have equal values.
             * \param other the BedRecord to compare the current instance with.
             * \returns True if all fields values are the same.
             */
            bool operator == (const BedRecord& other) const ;

            /*!
             * \brief Compares if two BedRecord have unequal values.
             * \param other the BedRecord to compare the current instance to.
             * \returns True if at least one fields value is different.
             */
            bool operator != (const BedRecord& other) const ;

            /*!
             * \brief Creates a string representation of the instance.
             * \returns the string representation of the record.
             */
            std::string toString() const ;

        public:
            /*!
             * \brief Defines the name of the BED line.
             */
            std::string name ;
            /*!
             * \brief A score attributed to this feature.
             */
            double score ;
    } ;


    /*!
    * \brief Writes the given record to a stream.
    * \param stream the stream to write the record to.
    * \param record the record to write.
    * \returns a reference to the stream after the writing.
    */
    std::ostream& operator << (std::ostream& stream,
                               const ngsai::BedRecord& record) ;

} // namespace ngsai


#endif // NGSAI_IO_BEDRECORD_HPP