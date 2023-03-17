#ifndef NGSAI_PAIRWISEKINETICMODEL_HPP
#define NGSAI_PAIRWISEKINETICMODEL_HPP


// forward declarations
namespace ngsai
{
    class RawKineticModel ;
    class NormalizedKineticModel ;
    class DiPositionKineticModel ;
    class DiPositionNormalizedKineticModel ;
    class PairWiseNormalizedKineticModel ;
}


#include <ngsaipp/epigenetics/KineticModel.hpp>
#include <ngsaipp/epigenetics/KineticSignal.hpp>
#include <vector>
#include <boost/histogram.hpp>


namespace bh = boost::histogram ;

/* a typedef for the 2D boost histograms to modelize 2D densities */
typedef bh::histogram<
    std::tuple<
        boost::histogram::axis::regular<
            double,
            boost::use_default,
            boost::use_default,
            boost::use_default>,
        boost::histogram::axis::regular<
            double,
            boost::use_default,
            boost::use_default,
            boost::use_default>
        >,
        boost::histogram::unlimited_storage<
            std::allocator<char>
        >
    > hist_2d_double ;



namespace ngsai
{
    /*!
     * \brief The PairWiseKineticModel class is made to model the raw IPD 
     * and PWD signal from PacBio CCS.
     * This model models a window of IPD/PWD signal of L consecutive bp. L is 
     * named the model size.
     * This model accounts for pair-wise relaitions ships between any two 
     * positions in the window. There are exactly N = ((L*L) - L) / 2 
     * pair-wise relations.
     * The signal in the Lbp window is modelled using 2*N different 
     * distributions of signal - N for IPD and N for PWD signal. Each 
     * distribution is modelled by a 2D histogram representing the signal at 
     * one position p in the window with respect to the signal found at 
     * the paired position q in the window where p != q.  Overall the window 
     * signal is modelled by 2*N histograms h_ipd_p_q and h_pwd_p_q : 
     * h_ipd_1_2, h_ipd_1_3, ..., h_ipd_1_L, h_ipd_2_3, ..., h_ipd_2_L, ..., 
     * h_ipd_L-1_L and h_pwd_1_2, h_pwd_1_3, ..., h_pwd_1_L, h_pwd_2_3, ..., 
     * h_pwd_2_L, ..., h_pwd_L-1_L.
     * The model can be trained (feed with data) in order to fill the 
     * underlying histograms and estimate the signal distributions at each 
     * pair of positions in the window. At this moment, the histograms are 
     * count histograms.
     * Once the histograms have been filled, the model can be normalized in 
     * order to transform the 2*N histograms into 2*N probability densities. 
     * At this moment, the model can compute the likelihood of observing a 
     * a strech of signal of length L with IPDs X x1, x2, ..., xL and PWDs Y  
     * y1, y2, ..., yL. The likelihood corresponds to :
     * p(x1,x2 | h_ipd_1_2) * p(x1,x3 | h_ipd_1_3) *... * p(x1,xL | h_ipd_1_L) 
     * * p(x2,x3 | h_ipd_2_3) * ... * p(x2,xL | h_ipd_2_L) * ... * 
     * p(xL-1,xL | h_ipd_L-1_L) * p(y1,y2 | h_pwd_1_2) * p(y1,y3 | h_pwd_1_3) *
     * ... * p(y1,yL | h_pwd_1_L) * p(y2,y3 | h_pwd_2_3) * ... * 
     * p(y2,yL | h_pwd_2_L) * ... * p(yL-1,yL | h_pwd_L-1_L) where 
     * p(xp, xq | h_p_q) means the probability of observing the value xp at 
     * the position p and xq at position q of the window, according to the 
     * density modeled in the histogram h_p_q.
     */
    class PairWiseKineticModel : public ngsai::KineticModel
    {   
        public:
            /*!
             * \brief Constructor. Creates an empty model.
             */
            PairWiseKineticModel() ;

            /*!
             * \brief Copy constructor.
             * \param other a reference to the object to copy from.
             */
            PairWiseKineticModel(const ngsai::PairWiseKineticModel& other) ;

            /*!
             * \brief Move constructor.
             * \param other a rvalue to the object to move.
             */
            PairWiseKineticModel(ngsai::PairWiseKineticModel&& other) ;

            /*!
             * \brief Destructor.
             */
            virtual
            ~PairWiseKineticModel() override ;
            
            /*!
             * \brief Assignment operator.
             * \param other a reference to the object to copy from.
             * \return a reference to the current instance.
             */
            ngsai::PairWiseKineticModel&
            operator = (const ngsai::PairWiseKineticModel& other) ;

            /*!
             * \brief Move assignment operator.
             * \param other a rvalue to the object to move.
             * \return a reference to the current instance.
             */
            ngsai::PairWiseKineticModel&
            operator = (ngsai::PairWiseKineticModel&& other) ;

            /*!
             * \brief Compares the content of the instances.
             * \param other an other instance to compare the current one 
             * with.
             * \returns whether both instances are equal.
             */
            bool 
            operator == (const ngsai::PairWiseKineticModel& other) const ;   

            /*!
             * \brief Compares the content of the instances.
             * \param other an other instance to compare the current one 
             * with.
             * \returns whether both instances are different.
             */
            bool 
            operator != (const ngsai::PairWiseKineticModel& other) const ;

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
             * Returns the histogram bin boundaries. 
             * All axes of all histograms have the same bin boundaries, so a 
             * list corresponding to a single axis only is returned.
             * \returns the bin boundaries. Each pair represents the interval 
             * [start,end) of a bin.
             */
            virtual
            std::vector<std::pair<double,double>>
            getBinBoundaries() const override ;

            /*!
             * \brief Inserts the averaged IPD and PWD (fw and rv) values 
             * in the model. Each element of the averaged IPD / PWD 
             * corresponds to a position in the model.
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
            ngsai::PairWiseKineticModel*
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
             * average values will be included in the model.
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

        protected:
            /*!
             * \brief Computes and sets the 0-based pairs of position index
             * (m_pair_index and m_pair_nb) that are to be modelled given the 
             * current model size (m_size), e.g. 0,3 means the position 0 and 
             * 3 in the window.
             */
            virtual
            void
            computeMeaningfulPairs() ;

        // for serialization
        protected:
            friend class boost::serialization::access ;

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
                archive & m_pair_nb ;
                archive & m_pair_index ;
                archive & m_histograms_ipd ;
                archive & m_histograms_pwd ;
            }

        protected:
            /*!
             * The number of pairs to model within the window.
             */
            size_t m_pair_nb ;
            /*!
             * The index of all pairs of positions to model within the 
             * window. The index are 0-based and correspond to positions in 
             * the window, e.g. 0,3 means position 0 and position 3.
             */
            std::vector<std::pair<size_t,size_t>> m_pair_index ;
            /*!
             * \brief the IPD signal histograms for each pair of positions 
             * in the model. Memory is allocated for all L**2 pairs but only 
             * the meaningfull (pairs stored in m_pair_index) will contain 
             * something. 
             */
            std::vector<std::vector<hist_2d_double>> m_histograms_ipd ;

            /*!
             * \brief the PWD signal histograms for each pair of positions 
             * in the model. Memory is allocated for all L**2 pairs but only 
             * the meaningfull (pairs stored in m_pair_index) will contain 
             * something. 
             */
            std::vector<std::vector<hist_2d_double>> m_histograms_pwd ;

    } ;  // PairWiseKineticModel

    /*!
     * \brief Stores the content of a 2D histogram into a matrix. The 
     * first dimension of the matrix represents the 1st axis and the 2nd 
     * dimension the 2nd axis. The histogram bins going to -inf and +inf are 
     * also put in the matrix (first and last on both dimensions).
     * \param h the histogram of interest.
     * \returns a matrix representation of the histogram.
     */
    std::vector<std::vector<double>> hist2d_to_matrix(const hist_2d_double& h) ;

} // namespace ngsai

#endif // NGSAI_PAIRWISEKINETICMODEL_HPP