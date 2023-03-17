#ifndef NGSAI_PAIRWISENORMALIZEDKINETICMODEL_HPP
#define NGSAI_PAIRWISENORMALIZEDKINETICMODEL_HPP


// forward declarations
namespace ngsai
{
    class RawKineticModel ;
    class NormalizedKineticModel ;
    class DiPositionKineticModel ;
    class DiPositionNormalizedKineticModel ;
}


#include <ngsaipp/epigenetics/PairWiseKineticModel.hpp>
#include <ngsaipp/epigenetics/KmerMap.hpp>
#include <ngsaipp/epigenetics/KineticSignal.hpp>
#include <ngsaipp/epigenetics/RawKineticModel.hpp>
#include <ngsaipp/epigenetics/NormalizedKineticModel.hpp>
#include <ngsaipp/epigenetics/DiPositionKineticModel.hpp>
#include <string>


namespace ngsai
{   
    // forward declaration
    class DiPositionNormalizedKineticModel ;

    /*!
     * \brief The NormalizedKineticModel class is made to model the  IPD 
     * and PWD signal from PacBio CCS.
     * This model models a window of IPD/PWD signal of L consecutive bp. L is 
     * named the model size.
     * Instead of modelling the raw kinetic signal like the 
     * PairWiseKineticModel class, this class models kinetic ratios 
     * (normalized). Each stretch of kinetic values (IPD or PWD) x1,x2,...,xL 
     * is normalized by the expected kinetic signal given the DNA sequence 
     * s1,s2,...,sL associated to the kinetic values. In the end, as for the 
     * PairWiseKineticModel class, the signal models contains 2*N different 
     * distributions of signal modelling all possible  N = ((L*L) - L) / 2 
     * pair-wise position relation between any position p and q where p != q :
     * h_ipd_1_2, h_ipd_1_3, ..., h_ipd_1_L, h_ipd_2_3, ..., h_ipd_2_L, ..., 
     * h_ipd_L-1_L and h_pwd_1_2, h_pwd_1_3, ..., h_pwd_1_L, h_pwd_2_3, ..., 
     * h_pwd_2_L, ..., h_pwd_L-1_L.
     * The model can be trained (feed with data) in order to fill the 
     * underlying histograms and estimate the signal distributions at each 
     * position in the window. If the model is feed with a vector of IPD 
     * values X x1, x2, ..., xL and a vector of PWD values Y y1, y2, ..., yL 
     * corresponding to a DNA sequence S s1, s2, ..., sL, each value xi and yi 
     * will be normalized into n(xi) n(yi) using a normalizing function 
     * accounting for the sequence context S. Then each histogram h_ipd_p_q and 
     * h_pwd_p_q will be added a count corresponding to the pairs of values 
     * n(xp), n(xq) an n(yp), n(yq) respectively.
     * Once the histograms have been filled, the model can be normalized in 
     * order to transform the 2*N histograms into 2*N probability densities. 
     * At this moment, the model can compute the likelihood of observing a 
     * a strech of signal of length L with IPDs X x1, x2, ..., xL and PWDs Y  
     * y1, y2, ..., yL corresponding to a sequence S s1, s2, ..., sL. The 
     * likelihood corresponds to :
     * p(n(x1),n(x2) | h_ipd_1_2) * p(n(x1),n(x3) | h_ipd_1_3) *... * 
     * p(n(x1),n(xL) | h_ipd_1_L) * p(n(x2),n(x3) | h_ipd_2_3) * ... * 
     * p(n(x2),n(xL) | h_ipd_2_L) * ... * p(n(xL-1),n(xL) | h_ipd_L-1_L) * 
     * p(n(y1),n(y2) | h_pwd_1_2) * p(n(y1),n(y3) | h_pwd_1_3) * ... * 
     * p(n(y1),n(yL) | h_pwd_1_L) * p(n(y2),n(y3) | h_pwd_2_3) * ... * 
     * p(n(y2),n(yL) | h_pwd_2_L) * ... * p(n(yL-1),n(yL) | h_pwd_L-1_L) where 
     * p(n(xp), n(xq) | h_p_q) means the probability of observing the 
     * normalized value n(xp) at the position p and n(xq) at position q of the 
     * window, according to the density modeled in the histogram h_p_q.
     */
    class PairWiseNormalizedKineticModel : public ngsai::PairWiseKineticModel
    {   
        public:
            
            /*!
             * \brief Default constructor. For conveniency only.
             * Such instances will never be able to perform any training nor 
             * loglikelihood computations. Use the constructor with KmerMap 
             * instead.
             */ 
            PairWiseNormalizedKineticModel() ;

            /*!
             * \brief Constructor.
             * \param bckg_model the model containing the expected mean IPD 
             * and PWD signal for each kmer.
             */
            PairWiseNormalizedKineticModel(const ngsai::KmerMap& bckg_model) ;

            /*!
             * \brief Copy constructor.
             * \param other a reference to the object to copy from.
             */
            PairWiseNormalizedKineticModel(
                const ngsai::PairWiseNormalizedKineticModel& other) ;

            /*!
             * \brief Move constructor.
             * \param other a rvalue to the object to move.
             */
            PairWiseNormalizedKineticModel(
                ngsai::PairWiseNormalizedKineticModel&& other) ;

            /*!
             * \brief Destructor.
             */
            virtual ~PairWiseNormalizedKineticModel() override ;
            
            /*!
             * \brief Assignment operator.
             * \param other a reference to the object to copy from.
             * \return a reference to the current instance.
             */
            ngsai::PairWiseNormalizedKineticModel&
            operator = (const ngsai::PairWiseNormalizedKineticModel& other) ;

            /*!
             * \brief Move assignment operator.
             * \param other a rvalue to the object to move.
             * \return a reference to the current instance.
             */
            ngsai::PairWiseNormalizedKineticModel&
            operator = (ngsai::PairWiseNormalizedKineticModel&& other) ;

            /*!
             * \brief Compares the content of the instances.
             * \param other an other instance to compare the current one 
             * with.
             * \returns whether both instances are equal.
             */
            bool 
            operator == (
                const ngsai::PairWiseNormalizedKineticModel& other) const ;

            /*!
             * \brief Compares the content of the instances.
             * \param other an other instance to compare the current one 
             * with.
             * \returns whether both instances are different.
             */
            bool 
            operator != (
                const ngsai::PairWiseNormalizedKineticModel& other) const ;

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
             * \param kinetics an instance containing the kinetics (IPD and 
             * PWD) on the forward and reverse strand to the model. The 
             * forward and reverse strand kinetics will be average and the 
             * average values will be included in the model.
             * \throws std::invalid argument if the size of kinetics does not 
             * match the size of the model.
             * \throws std::runtime_error if the model was not initialized 
             * using set_parameters().
             * \throws std::invalid_argument if the size of kinetics is not 
             * equal to the model size or if the data are incomplete 
             * (sequence, IPD or PWD missing).
             */
            virtual void
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
            virtual ngsai::PairWiseNormalizedKineticModel*
            copy() const override ;

            /*!
             * \brief Computes the natural loglikelihood - of the average 
             * kinetics (fw and rv) - given the model, p(ipds,pwds | model). 
             * This requires that the histograms have been normalized into 
             * log density histograms using density() and log(). 
             * \param kinetics an instance containing the kinetics (IPD and 
             * PWD) on the forward and reverse strand to the model. The 
             * forward and reverse strand kinetics will be average and the 
             * average values will be included in the model. The DNA sequences 
             * are not necessary.
             * \returns the likelihood of observing a strech of DNA with the 
             * given average IPDs and PWDs. Returns nan if the data are 
             * incomplete (sequence, IPD or PWD missing).
             * \throws std::runtime_error if the model was not initialized 
             * using set_parameters() or if the model has not been normalized 
             * to contain signal log densities using density() and log().
             * \throws std::invalid_argument if the kinetics size does not 
             * match the model size.
             */
            virtual double
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
            virtual size_t
            getKineticSignalRequiredSize() const override ;

            /*!
             * \brief Saves the instance using the serialization interface 
             * in the given file.
             * \param path the path to the file in which the instance will 
             * be saved.
             */
            virtual void
            save(const std::string& path) const override ;

            /*!
             * \brief Loads the instance using the serialization interface 
             * from the given file. All previous data are lost.
             * \param path the path to the file from which the instance will 
             * be loaded.
             */
            virtual void
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
            void
            serialize(Archive & archive, const unsigned int version)
            {   if(version < 1)
                { ; }
                
                //https://theboostcpplibraries.com/boost.serialization-class-hierarchies
                archive & boost::serialization::base_object<ngsai::PairWiseKineticModel>(*this) ; 
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



    } ;  // PairWiseNormalizedKineticModel

}  // namespace ngsai

#endif // NGSAI_PAIRWISENORMALIZEDKINETICMODEL_HPP