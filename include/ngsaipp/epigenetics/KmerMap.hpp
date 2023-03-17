#ifndef NGSAI_KMERMAP_HPP
#define NGSAI_KMERMAP_HPP

#include <string>
#include <utility>
#include <boost/serialization/access.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/serialization/string.hpp>
#include <boost/serialization/utility.hpp>

#include  <ngsaipp/epigenetics/KmerData.hpp>

namespace ngsai
{   
    // the KmerMap is a vector is these
    typedef std::pair<uint32_t,ngsai::KmerData> kmerBucket ;

    /*!
     * \brief The KmerMap is a template class allowing to create a map 
     * storing per kmer kinetic information. It is designed to allow O(k)
     * insertion and O(1) query operations using a Rabin-Karp hashing approach.
     * Each kmer is kinetic information is stored in a pair, together with the 
     * number of insertions in the map made for the corresponding kmer.
     * The memory footprint of the map is N**k where N is the number of 
     * character in the alphabet - the alphabet is defined in ngsai::base_code 
     * in dna/dna_utility.hpp) - and k is the kmer size. 
     */
    class KmerMap
    {   
        public:
            /*!
             * \brief The map insertion stategies. 
             * REPLACE : the kmer data already in the map are replaced with the 
             * inserted data. The number of times this kmer has been inserted 
             * is reset to 1.
             * SUM : the kmer data already in the map are summed with the 
             * inserted data. The number of times this kmer has been inserted 
             * is increased by 1.
             */
            enum insert_mode {REPLACE=0,
                              SUM,
                              N_MODES} ;
        public:

            /*!
             * \brief Default constructor. 
             * Creates an empty map.
             */
            KmerMap() ;

            /*!
             * \brief Construcs an empty map designed to contain N**kmer_size 
             * different kmers where N is the alphabet size. The Map memory 
             * allocation is performed here.
             * \param kmer_size the size of the kmers that will be 
             * stored.
             */
            KmerMap(size_t kmer_size) ;

            /*!
             * \brief Copy constructor.
             * \param other the instance to copy from.
             */
            KmerMap(const KmerMap& other) ;

            /*!
             * \brief Move constructor.
             * \param other the instance to move from.
             */
            KmerMap(KmerMap&& other) ;

            /*!
             * \brief Destructor.
             */
            ~KmerMap() ;
        
            /*!
             * \brief Assignment operator.
             * \param other A reference to the object to assign the value 
             * from.
             * \returns A reference to the current object.
             */
            ngsai::KmerMap& 
            operator = (const ngsai::KmerMap& other) ; 

            /*!
             * \brief Compares the content of the instances.
             * \param other an other instance to compare the current one 
             * with.
             * \returns whether both instances are equal.
             */
            bool 
            operator == (const ngsai::KmerMap& other) const ;   

            /*!
             * \brief Compares the content of the instances.
             * \param other an other instance to compare the current one 
             * with.
             * \returns whether both instances are different.
             */
            bool 
            operator != (const ngsai::KmerMap& other) const ;

            /*!
             * \brief Returns an iterator pointing to the first element 
             * of the map.
             * \returns An iterator to the first of the map.
             */
            std::vector<ngsai::kmerBucket>::iterator
            begin() ;

            /*!
             * \brief Returns an iterator pointing to the first element 
             * of the map.
             * \returns An iterator to the first of the map.
             */
            std::vector<ngsai::kmerBucket>::const_iterator
            begin() const ;

            /*!
             * \brief Returns an iterator pointing to past-the-last element 
             * of the map.
             * \returns An iterator to the end of the map.
             */
            std::vector<ngsai::kmerBucket>::iterator
            end() ;

            /*!
             * \brief Returns an iterator pointing to past-the-last element 
             * of the map.
             * \returns An iterator to the end of the map.
             */
            std::vector<ngsai::kmerBucket>::const_iterator
            end() const ;

            /*!
             * \brief Returns the number of kmers in the map.
             * \returns the number of kmers actually in the map.
             */
            size_t 
            size() const ;

            /*!
             * \brief Returns the total number of kmers that can be stored 
             * by the map.
             * \returns the total number of kmer that can be stored in the map.
             */
            size_t
            capacity() const ;

            /*!
             * \brief Resizes the map in order to store kmer of the given 
             * size. All previous data are lost.
             * \param kmer_size the new kmer size.
             */
            void
            resize(size_t kmer_size) ;

            /*!
             * \brief Searches the map for the kmer with the same sequence and 
             * returns an iterator pointing to it if such a kmer exists, 
             * otherwise returns an iterator pointing to the KmerMap::end. 
             * \returns An iterator to the element if it could be found, 
             * to KmerMap::end otherwise.
             */
            std::vector<ngsai::kmerBucket>::iterator 
            find(const KmerData& kmer) ;

            /*!
             * \brief Searches the map for the kmer with the same sequence and 
             * returns an iterator pointing to it if such a kmer exists, 
             * otherwise returns an iterator pointing to the KmerMap::end. 
             * \returns An iterator to the element if it could be found, 
             * to KmerMap::end otherwise.
             */
            std::vector<ngsai::kmerBucket>::const_iterator 
            find(const KmerData& kmer) const ;

            /*!
             * \brief Searches the map for each kmer in the sequence and
             * returns a vector of iterator to the corresponding kmers in the
             * map. If a kmer is not present in the map, the corresponding 
             * iterator points to KmerMap::end. The order of the iterators in 
             * the vector is the same as the order of the kmer occurence in 
             * the sequence.
             * \param sequence the sequence to search the kmers from.
             * \returns A vector of iterators pointing to the kmers in the map 
             * or to KmerMap::end for kmer not listed in the map.
             */
            std::vector<std::vector<ngsai::kmerBucket>::iterator> 
            find(const std::string& sequence) ;

             /*!
             * \brief Searches the map for each kmer in the sequence and
             * returns a vector of iterator to the corresponding kmers in the
             * map. If a kmer is not present in the map, the corresponding 
             * iterator points to KmerMap::end. The order of the iterators in 
             * the vector is the same as the order of the kmer occurence in 
             * the sequence.
             * \param sequence the sequence to search the kmers from.
             * \returns A vector of iterators pointing to the kmers in the map 
             * or to KmerMap::end for kmer not listed in the map.
             */
            std::vector<std::vector<ngsai::kmerBucket>::const_iterator> 
            find(const std::string& sequence) const ;


            /*!
             * \brief Inserts the given element in the map.
             * \param kmer a reference to the object to insert in the map.
             * \param mode the strategy to update the IPD/PWD if this kmer 
             * already exists in the map. If SUM is used, the kmer occurence 
             * counter is updated by 1.
             */
            void insert(const ngsai::KmerData& kmer,
                        insert_mode mode=insert_mode::REPLACE) ;

            /*!
             * \brief Inserts the given element in the map.
             * \param kmer a reference to the object to insert in the map.
             * \param mode the strategy to update the IPD/PWD if this kmer 
             * already exists in the map. If SUM is used, the kmer occurence 
             * counter is updated by 1.
             */
            void insert(const ngsai::KmerData&& kmer,
                        insert_mode mode=insert_mode::REPLACE) ;

            /*!
             * \brief Scans the sequence for kmers and insert each kmer 
             * inside the map. This is using Rabin-Karp's rolling hash for 
             * more efficient performances. The corresponding kmer IPDs/PDWs 
             * are extracted from the given IPDs/PWDs.
             * \param sequence the DNA sequence to get the kmers from.
             * \param ipd a vector of IPDs corresponding to each sequence 
             * position.
             * \param pwd a vector of PWDs corresponding to each sequence 
             * position.
             * \param mode the strategy to update the IPD/PWD if this kmer 
             * already exists in the map. If SUM is used, the kmer occurence 
             * counter is updated by 1.
             */
            void insert(const std::string& sequence,
                        const std::vector<uint32_t>& ipd,
                        const std::vector<uint32_t>& pwd,
                        insert_mode mode=insert_mode::REPLACE) ;

            /*!
             * \brief Returns a reference to the KmerData object with the
             * with the same sequence inside the map. If not such kmer exists 
             * in the map, this function throws an out_of_range exception.
             * \param a kmer with the sequence of interest.
             * \returns a reference to the KmerData object mapped 
             * with this sequence.
             */
            ngsai::kmerBucket& at(const ngsai::KmerData& kmer) ;

            /*!
             * \brief Returns a reference to the KmerData object with the
             * with the same sequence inside the map. If not such kmer exists 
             * in the map, this function throws an out_of_range exception.
             * \param a kmer with the sequence of interest.
             * \returns a reference to the KmerData object mapped 
             * with this sequence.
             */
            const ngsai::kmerBucket& at(const ngsai::KmerData& kmer) const ;

            /*!
             * \brief Returns the size of the kmers stored by this map.
             * \returns the size of the kmers in this map.
             */
            size_t getKmerSize() const ;
            
        private:
            // for serialization
            friend class boost::serialization::access;

            /*!
             * \brief (De-)serializes the object using Boost serialization 
             * library.
             * \param Archive the archive to write/read the object to/from.
             * \param version don't know.
             */
            template<class Archive>
            void serialize(Archive & archive, const unsigned int version)
            {   if(version < 1)
                { ; }
                archive & m_kmer_n ;
                archive & m_kmer_size ;
                archive & m_map ;
            }       
            
        protected:
            /*!
             * \brief The number of kmer actually in the map.
             */
            size_t m_kmer_n ;
            
            /*!
             * \brief The size of the kmer stored in the map.
             */
            size_t m_kmer_size ;

            /*!
             * \brief The vector containing the KmerData objects and the number 
             * insertions performed for this kmer. The vector has a length 
             * corresponding exactly to the number of possible kmer. Each 
             * object is inserted at the index given by a hash function 
             * working on the sequence.
             */
            std::vector<ngsai::kmerBucket> m_map ;
    } ;

} // namespace ngsai

/*!
 * \brief Writes the given KmerMap to the given stream.
 * \param stream a reference to the stream of interest.
 * \param m a KmerMap of interest.
 * \returns a reference to the stream after the writing.
 */
std::ostream& operator << (std::ostream& stream,
                               const ngsai::KmerMap& m) ;

#endif // NGSAI_KMERMAP_HPP