#ifndef NGSAI_DIPOSITIONKINETICMODEL_HPP
#define NGSAI_DIPOSITIONKINETICMODEL_HPP

// forward declarations
namespace ngsai
{
    class RawKineticModel ;
    class NormalizedKineticModel ;
    class DiPositionNormalizedKineticModel ;
    class PairWiseNormalizedKineticModel ;
}

#include <ngsaipp/epigenetics/PairWiseKineticModel.hpp>


namespace ngsai
{
    /*!
     * \brief The DiPositionKineticModel class is made to model the raw IPD
     * and PWD signal from PacBio CCS.
     * This model models a window of IPD/PWD signal of L consecutive bp. L is 
     * named the model size.
     * This model is a special case of the PairWiseKineticModel. The 
     * PairWiseKineticModel accounted for all N = ((L*L) - L) / 2  pair-wise 
     * position interaction whereas this one only considers L-1 interactions 
     * between any two directly neighbouring positions in the window. Except 
     * for this, everything is the same as in PairWiseKineticModel. 
     */
    class DiPositionKineticModel : public ngsai::PairWiseKineticModel
    {   
        public:
            /*!
             * \brief Constructor. Creates an empty model.
             */
            DiPositionKineticModel() ;

            /*!
             * \brief Copy constructor.
             * \param other a reference to the object to copy from.
             */
            DiPositionKineticModel(const ngsai::DiPositionKineticModel& other) ;

            /*!
             * \brief Move constructor.
             * \param other a rvalue to the object to move.
             */
            DiPositionKineticModel(ngsai::DiPositionKineticModel&& other) ;

            /*!
             * \brief Destructor.
             */
            virtual
            ~DiPositionKineticModel() override ;
            
            /*!
             * \brief Assignment operator.
             * \param other a reference to the object to copy from.
             * \return a reference to the current instance.
             */
            ngsai::DiPositionKineticModel&
            operator = (const ngsai::DiPositionKineticModel& other) ;

            /*!
             * \brief Move assignment operator.
             * \param other a rvalue to the object to move.
             * \return a reference to the current instance.
             */
            ngsai::DiPositionKineticModel&
            operator = (ngsai::DiPositionKineticModel&& other) ;

            /*!
             * \brief Compares the content of the instances.
             * \param other an other instance to compare the current one 
             * with.
             * \returns whether both instances are equal.
             */
            bool 
            operator == (
                const ngsai::DiPositionKineticModel& other) const ;

            /*!
             * \brief Compares the content of the instances.
             * \param other an other instance to compare the current one 
             * with.
             * \returns whether both instances are different.
             */
            bool 
            operator != (
                const ngsai::DiPositionKineticModel& other) const ;

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

            using PairWiseKineticModel::add ;

            /*!
             * \brief Creates a copy of the current instance and returns a
             * pointer to the copy.
             * \returns a pointer to the copy.
             */
            virtual
            ngsai::DiPositionKineticModel*
            copy() const override ;

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
       
    } ;  // DiPositionKineticModel

} // namespace ngsai

#endif // NGSAI_DIPOSITIONKINETICMODEL_HPP