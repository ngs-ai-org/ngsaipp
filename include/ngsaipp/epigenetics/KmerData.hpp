#ifndef NGSAI_KMERDATA_HPP
#define NGSAI_KMERDATA_HPP

#include <string>
#include <vector>
#include <stdexcept>
#include <boost/serialization/access.hpp>
#include <boost/serialization/vector.hpp>


#include <ngsaipp/io/utility.hpp>


namespace ngsai
{
    /*!
     * \brief The KmerData class is a class allowing to model a kmer together 
     * with per sequence position IPD and PWD data. Each of the sequence 
     * position corresponds to the same data vector index.
     */
    class KmerData
    {
        public:
            /*!
             * \brief Constructs an empty object.
             */
            KmerData() ;

            /*!
             * \brief Copy constructor.
             * \param other the instance to copy from.
             * \throws std::invalid_argument if the lengths of the sequence 
             * and the IPD/PWD of the given instance to copy from don't match.
             */
            KmerData(const KmerData& other) ;

            /*!
             * \brief Move constructor.
             * \param other the instance to copy from.
             * \throws std::invalid_argument if the lengths of the sequence 
             * and the IPD/PWD of the given instance don't match.
             */
            KmerData(KmerData&& other) ;

            /*!
             * \brief Constructs an object with the given sequence but 
             * 0 value IPD and PWD data.
             * \param seq the sequence.
             */
            KmerData(const std::string& seq) ;

            /*!
             * \brief Constructs an object with the given sequence, IPDs,  
             * PWDs and number of occurences.
             * \param seq the sequence.
             * \param ipd a vector of per position IPDs.
             * \param ipd a vector of per position PWDs.
             * \throws std::invalid_argument if the lengths of the sequence 
             * and the IPD/PWD of the given instance to copy from don't match.
             */
            KmerData(const std::string& seq,
                     const std::vector<uint32_t>& ipd,
                     const std::vector<uint32_t>& pwd) ;
            
            /*!
             * \brief Destructor
             */
            ~KmerData() ;

            /*!
             * \brief Assignment operator.
             * \param other A reference to the object to assign the value 
             * from.
             * \throws std::invalid_argument if the lengths of the sequence 
             * and the IPD/PWD of the given instance to copy from don't match.
             * \returns A reference to the current object.
             */
            ngsai::KmerData& operator = (const ngsai::KmerData& other) ;

            /*!
             * \brief Assignment move operator.
             * \param other A reference to the object to assign the value 
             * from.
             * \throws std::invalid_argument if the lengths of the sequence 
             * and the IPD/PWD of the given instance to copy from don't match.
             * \returns A reference to the current object.
             */
            ngsai::KmerData& operator = (const ngsai::KmerData&& other) ;

            /*!
             * \brief Equality operator.
             * \param other The object to compare the current instance with.
             * \returns true if both instance fields have the same values.
             */
            bool operator == (const ngsai::KmerData& other) const ;

            /*!
             * \brief Unequality operator.
             * \param other The object to compare the current instance with.
             * \returns true if at least one field has a different value among 
             * both instances.
             */
            bool operator != (const ngsai::KmerData& other) const ;

            /*!
             * \brief Returns the size of the kmer sequence.
             * \returns the size of the kmer sequence.
             */
            size_t size() const ;

        public:
            /*!
             * \brief the sequence.
             */
            std::string sequence ;

            /*!
             * \brief the vector of IPD values
             */
            std::vector<uint32_t> ipd ;

             /*!
             * \brief the vector of PWD values
             */
            std::vector<uint32_t> pwd ;
        
        private:
            // for serialization
            friend class boost::serialization::access;

            /*!
             * \brief (De-)serializes the object using Boost serialization 
             * library.
             * \param Archive the archive to write/read the object to/from.
             * \param version Useless, for overriding reasons.
             */
            template<class Archive>
            void serialize(Archive & archive, const unsigned int version)
            {   if(version < 1)
                { ; }
                archive & this->sequence ;
                archive & this->ipd ;
                archive & this->pwd ;
            }       
    } ;

}  // namespace ngsai

/*!
 * \brief Writes the KmerData object to the given stream.
 * \param stream the stream of interest.
 * \param kmer the KmerData object of interest.
 * \returns a reference to the stream.
 */
std::ostream& operator << (std::ostream& stream,
                           const ngsai::KmerData& kmer) ;

#endif // NGSAI_KMERDATA_HPP