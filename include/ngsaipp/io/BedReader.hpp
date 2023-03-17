#ifndef NGSAI_IO_BEDREADER_HPP
#define NGSAI_IO_BEDREADER_HPP

#include <string>
#include <fstream>

#include <ngsaipp/io/BedRecord.hpp>


namespace ngsai
{
    /*!
     * \brief A class to perform sequencial reading from a BED 
     * file using the getNext() method. Currently BED 6 format 
     * is handled.
     */
    class BedReader
    {
        public:

            BedReader() = delete ;

            /*!
            * \brief Constructs an object to read from the given file.
            * \param path the path to the file to read.
            */
            BedReader(const std::string& path) ;

            /*!
            * \brief Destructor. Closes the file stream if necessary.
            */
            ~BedReader() ;

            /*!
            * \brief Reads the next record in the file.
            * \param record a reference to an object in which the read data 
            * will be stored.
            * \returns whether a record could be read. False indicates EOF.
            * \throws std::runtime_error if a format error is detected.
            */
            bool getNext(BedRecord& record) ;

            /*!
            * \brief Closes the file stream.
            */
            void close() ;
        
        protected:
            /*!
            * \brief Opens the stream to the file.
            * \throws std::runtime_error if the file could not be open.
            */
            void open() ;
        
        protected:
            /*!
            * \brief The path to the file to read.
            */
            std::string m_path ;
            /*!
            * \brief The stream to the file to read.
            */
            std::ifstream m_f_bed ;
            /*!
            * \brief whether the stream to the file is open.
            */
            bool m_is_open ;
    } ;

} // namespace ngsai

#endif // NGSAI_IO_BEDREADER_HPP