#ifndef NGSAI_DIPOSITIONNORMALIZEDKINETICMODEL_HPP
#define NGSAI_DIPOSITIONNORMALIZEDKINETICMODEL_HPP


// forward declarations
namespace ngsai
{
    class RawKineticModel ;
    class NormalizedKineticModel ;
    class DiPositionKineticModel ;
    class DiPositionNormalizedKineticModel ;
    class PairWiseKineticModel ;
}


#include <ngsaipp/epigenetics/PairWiseNormalizedKineticModel.hpp>
#include <ngsaipp/epigenetics/KineticSignal.hpp>
#include <string>


namespace ngsai
{   /*!
     * \brief The DiPositionNormalizedKineticModel class is made to model the 
     * normalized IPD  and PWD signal from PacBio CCS.
     * Instead of modelling the raw kinetic signal like the 
     * DiPositionKineticModel class, this class models kinetic ratios 
     * (normalized). Each stretch of kinetic values (IPD or PWD) x1,x2,...,xL 
     * is normalized by the expected kinetic signal given the DNA sequence 
     * s1,s2,...,sL associated to the kinetic values.
     * This model models a window of IPD/PWD signal of L consecutive bp. L is 
     * named the model size.
     * This model is a special case of the PairWiseNormalizedKineticModel. The 
     * PairWiseNormalizedKineticModel accounted for all N = ((L*L) - L) / 2  
     * pair-wise position interaction whereas this one only considers L-1 
     * interactions between any two directly neighbouring positions in the 
     * window. Except for this, everything is the same as in 
     * PairWiseNormalizedKineticModel. 
     */
    class DiPositionNormalizedKineticModel : 
                                public ngsai::PairWiseNormalizedKineticModel
    {   
        public:
            /*!
             * \brief Constructor. Creates an empty model.
             */
            DiPositionNormalizedKineticModel() ;

            /*!
             * \brief Constructor.
             * \param ipd_model the model containing the expected mean IPD 
             * signal for each kmer.
             * \param pwd_model the model containing the expected mean PWD 
             * signal for each kmer.
             */
            DiPositionNormalizedKineticModel(const ngsai::KmerMap& bckg_model) ;

            /*!
             * \brief Copy constructor.
             * \param other a reference to the object to copy from.
             */
            DiPositionNormalizedKineticModel(
                const ngsai::DiPositionNormalizedKineticModel& other) ;

            /*!
             * \brief Move constructor.
             * \param other a rvalue to the object to move.
             */
            DiPositionNormalizedKineticModel(
                ngsai::DiPositionNormalizedKineticModel&& other) ;

            /*!
             * \brief Destructor.
             */
            virtual
            ~DiPositionNormalizedKineticModel() override ;
            
            /*!
             * \brief Assignment operator.
             * \param other a reference to the object to copy from.
             * \return a reference to the current instance.
             */
            ngsai::DiPositionNormalizedKineticModel&
            operator = (const ngsai::DiPositionNormalizedKineticModel& other) ;

            /*!
             * \brief Move assignment operator.
             * \param other a rvalue to the object to move.
             * \return a reference to the current instance.
             */
            ngsai::DiPositionNormalizedKineticModel&
            operator = (ngsai::DiPositionNormalizedKineticModel&& other) ;

            /*!
             * \brief Compares the content of the instances.
             * \param other an other instance to compare the current one 
             * with.
             * \returns whether both instances are equal.
             */
            bool 
            operator == (
                const ngsai::DiPositionNormalizedKineticModel& other) const ;

            /*!
             * \brief Compares the content of the instances.
             * \param other an other instance to compare the current one 
             * with.
             * \returns whether both instances are different.
             */
            bool 
            operator != (
                const ngsai::DiPositionNormalizedKineticModel& other) const ;

            /*!
             * \brief Creates a copy of the current instance and returns a
             * pointer to the copy.
             * \returns a pointer to the copy.
             */
            virtual
            ngsai::DiPositionNormalizedKineticModel*
            copy() const override ;

            using PairWiseNormalizedKineticModel::add ;

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
             * current model size (m_size), e.g. 1,2 means the position 1 and 
             * 2 in the window.
             */
            virtual
            void
            computeMeaningfulPairs() override ;

    } ;  // PairWiseNormalizedKineticModel

}  // namespace ngsai

#endif // NGSAI_DIPOSITIONNORMALIZEDKINETICMODEL_HPP