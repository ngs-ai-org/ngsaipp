#ifndef NGSAI_KINETICMODEL_HPP
#define NGSAI_KINETICMODEL_HPP

// forward declarations
namespace ngsai
{
    class RawKineticModel ;
    class NormalizedKineticModel ;
    class DiPositionKineticModel ;
    class DiPositionNormalizedKineticModel ;
    class PairWiseKineticModel ;
    class PairWiseNormalizedKineticModel ;
}

#include <iostream>
#include <vector>
#include <boost/histogram/serialization.hpp>

#include <ngsaipp/epigenetics/KineticSignal.hpp>


namespace ngsai
{   
    /*!
     * \brief The KineticModel class is an abstract class 
     * providing an interface to derive for classes that 
     * implement a statistical model of the kinetic signal 
     * of PacBio CCSs.
     */
    class KineticModel
    {   
        public:
            /*!
             * \brief Constructor. Creates an empty model.
             */
            KineticModel() ;

            /*!
             * \brief Constructor with parameters.
             * \param size the size of the model.
             * \param is_init whether the model is 
             * initialised.
             * \param is_density whether the model is a 
             * density.
             * \param is_log whether the model contains log 
             * transformed values.
             */
            KineticModel(size_t size,
                         bool is_init,
                         bool is_density,
                         bool is_log) ;
            
            /*!
             * \brief Copy constructor.
             * \param other the instance to copy from.
             */
            KineticModel(const ngsai::KineticModel& other) ;

            /*!
             * \brief Move constructor.
             * \param other the instance to move from.
             */
            KineticModel(ngsai::KineticModel&& other) ;

            /*!
             * \brief Destructor.
             */
            virtual 
            ~KineticModel() ;

            /*!
             * \brief Compares the content of the instances.
             * \param other an other instance to compare 
             * the current one with.
             * \returns whether both instances are equal.
             */
            bool 
            operator == (
                const ngsai::KineticModel& other) const ;

            /*!
             * \brief Compares the content of the instances.
             * \param other an other instance to compare 
             * the current one with.
             * \returns whether both instances are different.
             */
            bool 
            operator != (
                const ngsai::KineticModel& other) const ;

            /*!
             * \brief Sets the model parameters and 
             * initialize it to accept data.
             * \param size the number of positions (bases) 
             * that will be modelled.
             * \param xmin the lower bound of 1st bin in 
             * each of the histogram axis.
             * \param xmax the upper bound of the last bin 
             * in each of the histogram axis.
             * \param nbins the number of bins in each of 
             * the histogram axis.
             * \param pseudo_counts a pseudo count that 
             * will be added to each bin of the histograms 
             * upon initialization.
             * \throw std::invalid_argument if size is not 
             * bigger than 0, if x_min is not smaller than 
             * x_max or if n_bins is not bigger than 
             * 0.
             */
            virtual
            void
            setParameters(size_t size,
                          double xmin,
                          double xmax,
                          size_t n_bins,
                          double pseudo_counts) = 0 ;
        
            /*!
             * Returns the histogram bin boundaries. All 
             * histograms have the same bin boundaries.
             * \returns the bin boundaries. Each pair 
             * represents the interval [start,end) of a bin.
             */
            virtual
            std::vector<std::pair<double,double>>
            getBinBoundaries() const = 0 ;

            /*!
             * \brief Returns if the model has been 
             * initialised.
             * \returns whether the models has already been 
             * initialised.
             */
            bool
            isInit() const ;

            /*!
             * \brief Returns whether the model has already 
             * been normalized to contain signal densities. 
             * If not, it contains signal counts.
             * \returns whether the model contains signal 
             * densities.
             */
            bool
            isDensity() const ;

            /*!
             * \brief Returns whether the model contains 
             * log probability densities.
             * \returns whether the model contains signal 
             * densities.
             */
            bool
            isLog() const ;

            /*!
             * \brief Returns the model size, that is the 
             * length in bp of the signal modeled.
             * \returns the model size
             */
            size_t
            size() const ;

            /*!
             * \brief Creates a copy of the current 
             * instance and returns a pointer to the copy.
             * \returns a pointer to the copy.
             */
            virtual
            ngsai::KineticModel*
            copy() const = 0 ;

            /*!
             * \brief Must implement the normalization 
             * model such that it contains densities of
             * probabilities and turn the m_is_density to 
             * true.
             * \throws std::invalid_argument if the model 
             * was not initialized using set_parameters().
             */
            virtual
            void 
            density() = 0 ;

            /*!
            * \brief Must apply a natural log transformation 
            * to the model content and turn the m_is_log 
            * flag to true.
            * \throws std::invalid_argument if the model 
            * was not initialized using set_parameters().
            */
            virtual
            void 
            log() = 0 ;

            /*!
             * \brief Must apply an exponential 
             * transformation to the model content and turn 
             * the m_is_log flag to false.
             * \throws std::invalid_argument if the model 
             * was not initialized using set_parameters().
             */
            virtual
            void 
            exp() = 0 ;

            /*!
             * \brief Adds the given kinetic signal in the 
             * model.
             * \param kinetic the kinetic data of interest.
             */
            virtual
            void
            add(const ngsai::KineticSignal& kinetic) = 0 ;

            /*!
             * \brief Should add the content of the given 
             * model to the current model.
             * \param model the model of interest.
             * \throws std::invalid_argument if the given 
             * model cannot be added.
             * \throws std::runtime_error if another 
             * error occurs.
             */
            virtual
            void
            add(const ngsai::KineticModel& model) = 0 ;

            /*!
             * \brief Computes the natural loglikelihood of 
             * the kinetics given the model, p(ipds,pwds | model). 
             * This requires that the histograms have been 
             * normalized into log density histograms using 
             * density() and log().
             * \param ipds a vector of IPD values 
             * corresponding to each position in the model.
             * \param pwds a vector of PWD values 
             * corresponding to each position in the model.
             * \param sequence the DNA sequence 
             * corresponding to the IPD and PWD kinetics.
             * \returns the likelihood of observing the 
             * given IPDs and PWDs given the model.
             * \throws std::runtime_error if the model was 
             * not initialized using set_parameters() or if 
             * the model has not been normalized to contain 
             * signal log densities using density() and log().
             * \throws std::invalid_argument if the 
             * kinetics size does not match the model size.
             */
            virtual
            double
            logLikelihood(
                const ngsai::KineticSignal& kinetic) const = 0 ;

            /*!
             * \brief Returns the required size of the 
             * KineticSignal instances needed to use with 
             * add() or logLikelihood(). The KineticSignal
             * size is equal to the model size.
             * \returns the size (in bp) of the 
             * KineticSignal that can be used to train the 
             * model of compute a loglikelihood from.
             */
            virtual
            size_t
            getKineticSignalRequiredSize() const ;

            /*!
             * \brief Produces a string representation of the model
             * \returns the string representation of the model.
             */
            virtual
            std::string
            toString() const = 0 ;

            /*!
             * \brief Saves the instance using the 
             * serialization interface in the given file.
             * \param path the path to the file in which 
             * the instance will be saved.
             */
            virtual
            void
            save(const std::string& path) const = 0 ;

            /*!
             * \brief Loads the instance using the 
             * serialization interface from the given file. 
             * All previous data are lost.
             * \param path the path to the file from which 
             * the instance will be loaded.
             */
            virtual
            void
            load(const std::string& path) = 0 ;

        // double dispatch of add(const KineticModel&)
        public:
            /*!
             * \brief Should add the content of the current 
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
                ngsai::RawKineticModel& model) const = 0 ;
            
            /*!
             * \brief Should add the content of the current 
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
                ngsai::NormalizedKineticModel& model) const = 0 ;

            /*!
             * \brief Should add the content of the current 
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
                ngsai::PairWiseKineticModel& model) const = 0 ;

            /*!
             * \brief Should add the content of the current 
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
                ngsai::PairWiseNormalizedKineticModel& model) const = 0 ;

            /*!
             * \brief Should add the content of the current 
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
                ngsai::DiPositionKineticModel& model) const = 0 ;

            /*!
             * \brief Should add the content of the current 
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
                ngsai::DiPositionNormalizedKineticModel& model) const = 0 ;

        protected:

            // for serialization
            friend class boost::serialization::access;

            /*!
             * \brief (De-)serializes the object using 
             * Boost serialization library.
             * \param archive the archive to write/read the 
             * object to/from.
             * \param version don't know.
             */
            template <class Archive>
            void
            serialize(Archive& archive, const unsigned int version)
            {   if(version < 1)
                { ; }

                archive & m_size ;
                archive & m_is_init ;
                archive & m_is_density ;
                archive & m_is_log ;
            }
        
        protected:
            /*!
             * \brief the model size in bp.
             */
            size_t m_size ;
            /*!
             * \brief whether the model has been 
             * initialized (x min/max, number of bins, 
             * pseudo counts, etc).
             */
            bool m_is_init ;
            /*!
             * \brief true if the histograms contain 
             * densities (as opposed to counts).
             */
            bool m_is_density ;
            /*!
             * \brief true if the histograms contain 
             * log-transformed values.
             */
            bool m_is_log ;

    } ; // KineticModel

}  // namespace ngsai


#endif // NGSAI_KINETICMODEL_HPP