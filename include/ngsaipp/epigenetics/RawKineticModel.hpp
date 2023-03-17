#ifndef NGSAI_RAWKINETICMODEL_HPP
#define NGSAI_RAWKINETICMODEL_HPP


// forward declarations
namespace ngsai
{
    class NormalizedKineticModel ;
    class DiPositionKineticModel ;
    class DiPositionNormalizedKineticModel ;
    class PairWiseKineticModel ;
    class PairWiseNormalizedKineticModel ;
}


#include <ngsaipp/epigenetics/KineticModel.hpp>
#include <ngsaipp/epigenetics/KmerMap.hpp>
#include <ngsaipp/epigenetics/KineticSignal.hpp>
#include <vector>
#include <string>
#include <boost/histogram.hpp>
#include <boost/histogram/serialization.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/serialization/string.hpp>
#include <boost/serialization/utility.hpp>


namespace bh = boost::histogram ;


/* a typedef for the 1D boost histograms to modelize 1D densities */
typedef bh::histogram<
            std::tuple<
                boost::histogram::axis::regular<
                    double, 
                    boost::use_default, 
                    boost::use_default, 
                    boost::use_default>
            >,
            boost::histogram::unlimited_storage<
                std::allocator<char>
            >
        > hist_1d_double ;


namespace ngsai
{   

    /*!
     * \brief The RawKineticModel class is made to model the raw IPD and PWD 
     * signal from PacBio CCS.
     * This model models a window of IPD/PWD signal of L consecutive bp. L is 
     * named the model size.
     * The signal in the Lbp window is modelled using 2*L different 
     * distributions of signal - L for IPD and L for PWD signal.
     * The distributions at each position are empiric. Overall the window 
     * signal is modelled by 2*L histograms h_ipd(1), h_ipd(2), ..., h_ipd(L) 
     * and h_pwd(1), h_pwd(2), ..., h_pwd(L).
     * The model can be trained (feed with data) in order to fill the 
     * underlying histograms and estimate the signal distributions at each 
     * position in the window. At this moment, the histograms are count 
     * histograms. If the model is feed with a vector of IPD values X x1, 
     * x2, ..., xL and a vector of PWD values Y y1, y2, ..., yL then each 
     * histogram h_ipd(i) and h_pwd(i) will be added a count corresponding to 
     * the value xi an yi respectively. h_ipd(1) will be added a count 
     * corresponding to x1, h_ipd(2) will be added a count corresponding to 
     * x2, ..., h_ipd(L) will be added a count corresponding to xL. The same 
     * will apply for the PWD histograms with Y values.
     * Once the histograms have been filled, the model can be normalized in 
     * order to transform the 2*L histograms into 2*L probability densities. 
     * At this moment, the model can compute the likelihood of observing a 
     * a strech of signal of length L with IPDs X x1, x2, ..., xL and PWDs Y  
     * y1, y2, ..., yL. The likelihood corresponds to :
     * p(x1 | h_ipd(1)) * p(y1 | h_pwd(1)) *... * p(xL | h_ipd(L)) * 
     * p(xL | h_pwd(L)) where p(xi | h(i)) means the 
     * probability of observing the value x at the i-th position of the 
     * window, according to the density modeled in the histogram h(i).
     */
    class RawKineticModel : public ngsai::KineticModel
    {  
        public:
            /*!
             * \brief Constructor. Creates an empty model.
             */
            RawKineticModel() ;

            /*!
             * \brief Copy constructor.
             * \param other a reference to the object to copy from.
             */
            RawKineticModel(const ngsai::RawKineticModel& other) ;

            /*!
             * \brief Move constructor.
             * \param other a rvalue to the object to move.
             */
            RawKineticModel(ngsai::RawKineticModel&& other) ;

            /*!
             * \brief Destructor.
             */
            virtual
            ~RawKineticModel() override ;
            
            /*!
             * \brief Assignment operator.
             * \param other a reference to the object to copy from.
             * \return a reference to the current instance.
             */
            ngsai::RawKineticModel&
            operator = (const ngsai::RawKineticModel& other) ;

            /*!
             * \brief Move assignment operator.
             * \param other a rvalue to the object to move.
             * \return a reference to the current instance.
             */
            ngsai::RawKineticModel&
            operator = (ngsai::RawKineticModel&& other) ;

            /*!
             * \brief Compares the content of the instances.
             * \param other an other instance to compare the current one 
             * with.
             * \returns whether both instances are equal.
             */
            bool 
            operator == (const ngsai::RawKineticModel& other) const ;   

            /*!
             * \brief Compares the content of the instances.
             * \param other an other instance to compare the current one 
             * with.
             * \returns whether both instances are different.
             */
            bool 
            operator != (const ngsai::RawKineticModel& other) const ;

            /*!
             * \brief Sets the model parameters and initialize it to accept 
             * data.
             * \param size the number of positions (bases) that will be 
             * modelled.
             * \param xmin the lower bound of 1st bin in the histograms.
             * \param xmax the upper bound of the last bin in the histograms.
             * \param nbins the number of bins in the histograms.
             * \param pseudo_counts a pseudo count that will be added to each 
             * bin of the histograms upon initialization.
             * \throw std::invalid_argument if size is not bigger than 0, if 
             * x_min is not smaller than x_max or if n_bins is not bigger than 
             * 0.
             */
            virtual
            void
            setParameters(size_t size,
                          double xmin,
                          double xmax,
                          size_t n_bins,
                          double pseudo_counts) override ;
            
            /*!
             * Returns the histogram bin boundaries. All histograms
             * have the same bin boundaries.
             * \returns the bin boundaries. Each pair represents the interval 
             * [start,end) of a bin.
             */
            virtual
            std::vector<std::pair<double,double>>
            getBinBoundaries() const override ;
           
            /*!
             * \brief Inserts the given kinetic values at the given position 
             * in the model.
             * \param ipd the IPD value to insert in the model.
             * \param pwd the PWD value to insert in the model.
             * \param position the position at which the given kinetic should 
             * be inserted in the model.
             * \throws std::invalid_argument if the position is not compatible 
             * with the size of the model.
             * \throws std::runtime_error if the model was not initialized 
             * using set_parameters().
             * \throws std::out_of_range if the position is out of range.
             */
            void
            add(double ipd, double pwd, size_t position) ;

            /*!
             * \brief Inserts the averaged IPD and PWD (fw and rv) values 
             * in the model. Each element of the averaged IPD / PWD 
             * corresponds to a position in the model.
             * \param kinetics an instance containing the kinetics (IPD and 
             * PWD) on the forward and reverse strand to the model. The 
             * forward and reverse strand kinetics will be average and the 
             * average values will be included in the model. The DNA sequences 
             * are not necessary.
             * \throws std::invalid argument if the size of kinetics does not 
             * match the size of the model.
             * \throws std::runtime_error if the model was not initialized 
             * using set_parameters().
             * \throws std::invalid_argument if the size of kinetics is not 
             * equal to the model size or if the data are incomplete 
             * (sequence, IPD or PWD missing).
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
            ngsai::RawKineticModel*
            copy() const override ;

            /*!
             * \brief Normalizes the count histograms to density histograms.
             * \throws std::invalid_argument if the model was not initialized 
             * using set_parameters().
             */
            virtual
            void
            density() override ;

            /*!
            * \brief Applies a natural log transformation to the histogram 
            * content.
            * \throws std::invalid_argument if the model was not initialized 
            * using set_parameters().
            */
            virtual
            void 
            log() override ;

            /*!
             * \brief Applies an exponential transformation to the histogram 
             * content.
             * \throws std::invalid_argument if the model was not initialized 
             * using set_parameters().
             */
            virtual
            void 
            exp() override ;

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
            virtual
            double
            logLikelihood(const ngsai::KineticSignal& kinetics) const override ;
            
            /*!
             * \brief Produces a string representation of the model
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
        
        public:
            friend std::ostream& 
            operator << (std::ostream& stream,
                         const ngsai::RawKineticModel& model) ; 

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
                archive & boost::serialization::base_object<ngsai::KineticModel>(*this) ; 
                archive & m_histograms_ipd ;
                archive & m_histograms_pwd ;
            }

        protected:
            /*!
             * \brief the IPD signal histograms at each position of the model.
             */
            std::vector<hist_1d_double> m_histograms_ipd ;

            /*!
             * \brief the PWD signal histograms at each position of the model.
             */
            std::vector<hist_1d_double> m_histograms_pwd ;
    } ;

   
    /*!
     * \brief Writes the given RawKineticModel to the output stream. Does 
     * nothing if the model is not initialized.
     * \param stream an output stream to write to.
     * \param model a model of interest.
     * \returns a reference to the output stream the model was written to.
     */
    std::ostream& 
    operator << (std::ostream& stream, const ngsai::RawKineticModel& model) ;
}

#endif // NGSAI_RAWKINETICMODEL_HPP