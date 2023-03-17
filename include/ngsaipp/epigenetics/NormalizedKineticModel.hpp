#ifndef NGSAI_NORMALIZEDKINETICMODEL_HPP
#define NGSAI_NORMALIZEDKINETICMODEL_HPP


// forward declarations
namespace ngsai
{
    class DiPositionKineticModel ;
    class DiPositionNormalizedKineticModel ;
    class PairWiseKineticModel ;
    class PairWiseNormalizedKineticModel ;
}


#include <ngsaipp/epigenetics/RawKineticModel.hpp>
#include <ngsaipp/epigenetics/KmerMap.hpp>
#include <ngsaipp/epigenetics/KineticSignal.hpp>
#include <vector>
#include <string>


namespace ngsai
{   
    /*!
     * \brief The NormalizedKineticModel class is made to model the IPD 
     * and PWD signal from PacBio CCS.
     * This model models a window of IPD/PWD signal of L consecutive bp. L is 
     * named the model size. 
     * Instead of modelling the raw kinetic signal like the RawKineticModel 
     * class, this class models kinetic ratios (normalized). Each stretch of 
     * kinetic values (IPD or PWD) x1,x2,...,xL is normalized by the expected 
     * kinetic signal given the DNA sequence s1,s2,...,sL associated to the 
     * kinetic values. In the end, as for the RawKineticModel class, 
     * the signal models contains 2*L different distributions of signal.
     * The distributions at each position are empiric and are made of 
     * 2*L histograms h_ipd(1), h_pwd(1), h_ipd(2), h_pwd(2) ..., h_ipd(L), 
     * h_pwd(L).
     * The model can be trained (feed with data) in order to fill the 
     * underlying histograms and estimate the signal distributions at each 
     * position in the window. If the model is feed with a vector of IPD 
     * values X x1, x2, ..., xL and a vector of PWD values Y y1, y2, ..., yL 
     * corresponding to a DNA sequence S s1, s2, ..., sL, each value xi and yi 
     * will be normalized into n(xi) n(yi) using a normalizing function 
     * accounting for the sequence context S. Then each histogram h_ipd(i) and 
     * h_pwd(i) will be added a count corresponding to the value n(xi) an 
     * n(yi) respectively. h_ipd(1) will be added a count corresponding to 
     * n(x1), h_ipd(2) will be added a count corresponding to n(x2), ..., 
     * h_ipd(L) will be added a count corresponding to n(xL). The same will 
     * apply for the PWD histograms with normalized Y values.
     * Once the histograms have been filled, the model can be normalized in 
     * order to transform the 2*L histograms into 2*L probability densities. 
     * At this moment, the model can compute the likelihood of observing a 
     * a strech of signal of length L with IPDs X x1, x2, ..., xL and PWDs Y  
     * y1, y2, ..., yL corresponding to a sequence S s1, s2, ..., sL. The 
     * likelihood corresponds to :
     * p(n(x1) | h_ipd(1)) * p(n(y1) | h_pwd(1)) *... * p(n(xL) | h_ipd(L)) * 
     * p(n(xL) | h_pwd(L)) where p(n(xi) | h(i)) means the 
     * probability of observing the normalized value n(x) at the i-th position 
     * of the window, according to the density modeled in the histogram h(i).
     */
    class NormalizedKineticModel : public ngsai::RawKineticModel
    {
        public:

            /*!
             * \brief Default constructor. For conveniency only.
             * Such instances will never be able to perform any training nor 
             * loglikelihood computations. Use the constructor with KmerMap 
             * instead.
             */ 
            NormalizedKineticModel() ;

            /*!
             * \brief Constructor.
             * \param bckg_model the model containing the expected mean IPD 
             * and PWD signal for each kmer.
             */
            NormalizedKineticModel(const ngsai::KmerMap& bckg_model) ;

            /*!
             * \brief Copy constructor.
             * \param other a reference to the object to copy from.
             */
            NormalizedKineticModel(const ngsai::NormalizedKineticModel& other) ;

            /*!
             * \brief Move constructor.
             * \param other a rvalue to the object to move.
             */
            NormalizedKineticModel(ngsai::NormalizedKineticModel&& other) ;

            /*!
             * \brief Destructor.
             */
            virtual
            ~NormalizedKineticModel() override ;

            /*!
             * \brief Assignment operator.
             * \param other a reference to the object to copy from.
             * \return a reference to the current instance.
             */
            ngsai::NormalizedKineticModel&
            operator = (const ngsai::NormalizedKineticModel& other) ;

            /*!
             * \brief Move assignment operator.
             * \param other a rvalue to the object to move.
             * \return a reference to the current instance.
             */
            ngsai::NormalizedKineticModel&
            operator = (ngsai::NormalizedKineticModel&& other) ;

            /*!
             * \brief Compares the content of the instances.
             * \param other an other instance to compare the current one 
             * with.
             * \returns whether both instances are equal.
             */
            bool 
            operator == (const ngsai::NormalizedKineticModel& other) const ;   

            /*!
             * \brief Compares the content of the instances.
             * \param other an other instance to compare the current one 
             * with.
             * \returns whether both instances are different.
             */
            bool 
            operator != (const ngsai::NormalizedKineticModel& other) const ;

            /*!
             * \brief Inserts the normalized averaged IPD and PWD values in 
             * the model. First the IPD and PWD are normalized (IPDRs / PWDRs) 
             * using the model kmer map and then averaged (fw and rv). 
             * Each element of the averaged IPDRs / PWDRs corresponds to a 
             * position in the model.
             * Note that the normalization process shortens the IPDRs / PWDRs 
             * signal stretches of the first K//2 and last K//2 values where 
             * K is the size of the kmer used for the normalization. 
             * Thus the kinetic values and sequence must have a length equal 
             * to S+2*K where S is the size of the model.  
             * \param kinetics an instance containing the kinetics (IPDs and 
             * PWDs) on the forward and reverse strand to the model. The 
             * forward and reverse strand kinetics will be average and the 
             * average values will be included in the model.
             * \throws std::invalid argument if the size of kinetics / 
             * sequence does not match the size of the model.
             * \throws std::runtime_error if the model was not initialized 
             * using set_parameters().
             * \throws std::invalid argument if the size of kinetics / 
             * sequences do not match the size of the model if the data are 
             * incomplete (sequence, IPD or PWD missing).
             */
            virtual
            void
            add(const ngsai::KineticSignal& kinetics) override ;

            /*!
             * \brief Adds the content of the given 
             * model to the current model.
             * \param model the model of interest.
             * \throws std::invalid_argument if the given 
             * model cannot be added.
             * \throws std::runtime_error if another 
             * error occurs.
             */
            virtual
            void
            add(const ngsai::KineticModel& model) override ;
                        
            /*!
             * \brief Creates a copy of the current instance and returns a
             * pointer to the copy.
             * \returns a pointer to the copy.
             */
            virtual
            ngsai::NormalizedKineticModel* 
            copy() const override ;

            /*!
             * \brief Computes the log likelihood of the normalized averaged 
             * kinetics given the model, p(ipdrs,pwdrs | model). First the IPDs 
             * and PWDs are normalized (IPDRs / PWDRs) using the model kmer 
             * map and then averaged (fw and rv). Each element of the averaged 
             * IPDRs / PWDRs corresponds to a position in the model.
             * Note that the normalization process shortens the IPDRs / PWDRs 
             * signal stretches of the first K//2 and last K//2 values where 
             * K is the size of the kmer used for the normalization. 
             * Thus the kinetic values and sequence must have a length equal 
             * to S+2*K where S is the size of the model.
             * This requires that the histograms have been normalized into log 
             * density histograms using density() and log().
             * \param kinetics an instance containing the kinetics (IPDs and 
             * PWDs) on the forward and reverse strand to the model. The 
             * forward and reverse strand kinetics will be normalized and 
             * averaged and the averaged values will used to compute the log 
             * likelihood.
             * \returns the likelihood of observing a strech of DNA with the 
             * given averaged IPDRs and PWDRs. Returns nan if the data are 
             * incomplete (sequence, IPD or PWD missing).
             * \throws std::runtime_error if the model was not initialized 
             * using set_parameters() or if the model has not been normalized 
             * to contain signal log densities using density() and log().
             * \throws std::invalid argument if the size of kinetics / 
             * sequences do not match the size.
             */
            virtual
            double
            logLikelihood(const ngsai::KineticSignal& kinetics) const override ;

            /*!
             * \brief Returns the required size of the KineticSignal instances 
             * needed to use with add() or logLikelihood(). The KineticSignal 
             * size is actually bigger than the model size. The normalization 
             * process trims exactly K//2 bp on each edge of the KineticSignal 
             * given, where K is the size of the kmer in the background model.
             * The required KineticSignal size is equal to L + 2*(K//2) where 
             * L is the size of the model.
             * \returns the size (in bp) of the KineticSignal that can be 
             * used to train the model of compute a loglikelihood from.
             */
            virtual
            size_t
            getKineticSignalRequiredSize() const override ;
            
            /*!
             * \brief Produces a string representation of the model. The 
             * background model is not included, only the histograms are.
             * \returns the string representation of the model.
             */
            virtual
            std::string
            toString() const override ;

            /*!
             * \brief Saves the instance using the serialization interface 
             * in the given file.
             * \param path the path to the file in which the instance will 
             * be saved.
             */
            virtual
            void
            save(const std::string& path) const override ;

            /*!
             * \brief Loads the instance using the serialization interface 
             * from the given file. All previous data are lost.
             * \param path the path to the file from which the instance will 
             * be loaded.
             */
            virtual
            void
            load(const std::string& path) override ;
        
        // double dispatch of add(const KineticModel&)
        public:
            /*!
             * \brief Adds the content of the current 
             * model to the given model.
             * \param model the model of interest.
             * \throws std::invalid_argument if the given 
             * model cannot be added.
             * \throws std::runtime_error if another 
             * error occurs.
             */
            virtual
            void
            addTo(
                ngsai::RawKineticModel& model) const override ;
            
            /*!
             * \brief Adds the content of the current 
             * model to the given model.
             * \param model the model of interest.
             * \throws std::invalid_argument if the given 
             * model cannot be added.
             * \throws std::runtime_error if another 
             * error occurs.
             */
            virtual
            void
            addTo(
                ngsai::NormalizedKineticModel& model) const override ;

            /*!
             * \brief Adds the content of the current 
             * model to the given model.
             * \param model the model of interest.
             * \throws std::invalid_argument if the given 
             * model cannot be added.
             * \throws std::runtime_error if another 
             * error occurs.
             */
            virtual
            void
            addTo(
                ngsai::PairWiseKineticModel& model) const override ;

            /*!
             * \brief Adds the content of the current 
             * model to the given model.
             * \param model the model of interest.
             * \throws std::invalid_argument if the given 
             * model cannot be added.
             * \throws std::runtime_error if another 
             * error occurs.
             */
            virtual
            void
            addTo(
                ngsai::PairWiseNormalizedKineticModel& model) const override ;

            /*!
             * \brief Adds the content of the current 
             * model to the given model.
             * \param model the model of interest.
             * \throws std::invalid_argument if the given 
             * model cannot be added.
             * \throws std::runtime_error if another 
             * error occurs.
             */
            virtual
            void
            addTo(
                ngsai::DiPositionKineticModel& model) const override ;

            /*!
             * \brief Adds the content of the current 
             * model to the given model.
             * \param model the model of interest.
             * \throws std::invalid_argument if the given 
             * model cannot be added.
             * \throws std::runtime_error if another 
             * error occurs.
             */
            virtual
            void
            addTo(
                ngsai::DiPositionNormalizedKineticModel& model) const override ;
        
        // for serialization
        protected:
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

                //https://theboostcpplibraries.com/boost.serialization-class-hierarchies
                archive & boost::serialization::base_object<ngsai::RawKineticModel>(*this) ; 
                archive & m_has_kmermap ;
                archive & m_bckg_model ;
            }

        protected:
            /*!
             * \brief whether the instance has a KmerMap
             */
            bool m_has_kmermap ;
            /*!
             * \brief the kmer background model to use to normalize the 
             * raw kinetics.
             */
            ngsai::KmerMap m_bckg_model ;
            
    } ;  // NormalizedKineticModel

} // namespace ngsai


#endif // NGSAI_NORMALIZEDKINETICMODEL_HPP