#ifndef NGSAI_KINETICCLASSIFIER_HPP
#define NGSAI_KINETICCLASSIFIER_HPP

#include <string>
#include <vector>
#include <list>
#include <boost/histogram.hpp>
#include <pbbam/BamRecord.h>

#include <ngsaipp/epigenetics/KineticModel.hpp>
#include <ngsaipp/epigenetics/KineticSignal.hpp>
#include <ngsaipp/genome/CpGRegion.hpp>



namespace ngsai
{
    /*!
     * The KineticClassifier class implements a baysian classifier that 
     * determine the methylation state of CpGs from PacBio CCS data. The 
     * classifier implements a 2 class model (methylated and unmethylated) 
     * and assigns a probability of belonging to each class, for each CpG, 
     * given the kinetics of the CCS mapping over the CpG region.
     */
    class KineticClassifier
    {
        public:
            
            /*!
             * \brief Constructor.
             */
            KineticClassifier() ;

            /*!
             * \brief Constructs a classifier using the given kinetic signal 
             * models for predicion. It takes ownership over the KineticModels.
             * \param model_meth the model for the methylated kinetic signal.
             * \param model_unmeth the model for the unmethylated kinetic 
             * signal.
             * \throws std::invalid_argument if the models don't have an 
             * identical size. 
             */
            KineticClassifier(ngsai::KineticModel* model_meth,
                              ngsai::KineticModel* model_unmeth) ;

            /*!
             * \brief Copy constructor.
             * \param other an instance to copy from.
             */
            KineticClassifier(const ngsai::KineticClassifier& other) ;

            /*!
             * \brief Move constructor.
             * \param other a rvalue to the object to move.
             */
            KineticClassifier(ngsai::KineticClassifier&& other) ;
            
            /*!
             * \brief Destructor. It will free all KineticModel memory!
             */
            ~KineticClassifier() ;

            /*!
             * \brief Assignment operator.
             * \param other an instance to assign from.
             * \returns A reference to the instance.
             */
            KineticClassifier&
            operator = (const KineticClassifier& other) ;

            /*!
             * \brief Move assignment operator.
             * \param other an instance to assign from.
             * \returns A reference to the instance.
             */
            KineticClassifier&
            operator = (KineticClassifier&& other) ;
    
            /*!
             * \brief Sets the classifier models.
             * \param model_meth the signal model for the methylated kinetics.
             * \param model_unmeth the signal model for the unmethylated 
             * kinetics.
             * \throws std::invalid_argument if all models are not initialised 
             * or if all models don't have the same size.
             */
            void
            setModels(ngsai::KineticModel* model_meth,
                      ngsai::KineticModel* model_unmeth) ;
            
            /*!
             * \brief Computes the posterior probabilities 
             * p(methylated | IPD,PWD) and p(unmethylated | IPD,PWD) for the 
             * given region, given the CCS read kinetics. The prior 
             * class probabilites p(methylated) and p(unmethylated) are 
             * normalized to sum up to 1.
             * \param cpg the coordinates of a CpG to classify, in genome 
             * coordinates. The CpG on both strand will be considered. The 
             * coordinates must be strandless. Start must be correspond to the 
             * coordinate of the C on the + strand and end to passed the G 
             * coordinate on the + strand.
             * \param ccss a list of CCS reads mapping on the CpG of interest. 
             * Their IPD and PWD kinetics will be used to compute the 
             * probability that the CpG is methylated or not methylated. 
             * CCS not mapping over the CpG or without a perfect alignment 
             * will be ignored.
             * \param prob_meth the a priori probability of the methylated 
             * class.
             * \param prob_unmeth the a priori probability of the unmethylated 
             * class.
             * \return a pair of posterior probability p(methylated | IPD,PWD) 
             * and p(unmethylated | IPD,PWD) respectively. The probabilities 
             * sum to 1.
             */
            std::pair<double, double>
            classify(const ngsai::genome::CpGRegion& cpg,
                     const std::list<PacBio::BAM::BamRecord>& ccss,
                     double prob_meth,
                     double prob_unmeth) const ;
            
            /*!
             * \brief Computes the posterior probabilities 
             * p(methylated | IPD,PWD) and p(unmethylated | IPD,PWD) for the 
             * given window signal.
             * \param kinetics an instance containing the kinetics (IPD and 
             * PWD) and DNA sequences on the forward and reverse strand of the 
             * region to classify.
             * \param prob_meth the a priori probability of the methylated 
             * class.
             * \param prob_unmeth the a priori probability of the unmethylated 
             * class.
             * \return a pair of posterior probability p(methylated | IPD,PWD) 
             * and p(unmethylated | IPD,PWD) respectively. The probabilities 
             * sum to 1.
             */
            std::pair<double, double>
            classify(const ngsai::KineticSignal& kinetics,
                     double prob_meth,
                     double prob_unmeth) const ;

        protected:
            /*!
             * \brief Checks that all model size are the same and that the 
             * models contain log probability densities, hence whether they 
             * can be used to compute likelihood.
             * \throws std::runtime_error if this is not the case.
             */ 
            void
            checkModels() const ;

            /*!
             * \brief Frees the model memory if the pointer point
             * to something.
             */
            void 
            freeModels() ;

        protected:
            /*!
             * \brief The window size in bases.
             */
            size_t m_window_size ;
            /*!
             * \brief The methylated kinetic signal model.
             */
            ngsai::KineticModel* m_model_meth ;
            /*!
             * \brief The unmethylated kinetic signal model.
             */
            ngsai::KineticModel* m_model_unmeth ;
    } ;

} // namespace ngsai


#endif // NGSAI_KINETICCLASSIFIER_HPP