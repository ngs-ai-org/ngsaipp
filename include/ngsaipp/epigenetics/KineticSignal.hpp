#ifndef NGSAI_KINETICSIGNAL_HPP
#define NGSAI_KINETICSIGNAL_HPP


#include <vector>
#include <string>


namespace ngsai
{   
    /*!
     * \brief the KineticSignal class is designed to be a data structure to 
     * store the DNA sequence and kinetic - IPD and PWD - information in a 
     * given strech of DNA. 
     * An instance has a size, that is the length of the DNA strech, in bp, 
     * that is modelled. The default contructor allow to create empy - 0bp -
     * streches and to add data later using the setter methods. In this case, 
     * the size is set at the 1st use of a setter. Except this 1st insertion -
     * on an empty instance - data inserted must have the same size as the 
     * instance.
     */
    class KineticSignal
    {   
        public:
            /*!
             * \brief Default constructor. Constructs an empty instance.
             */
            KineticSignal() ;

            /*!
             * \brief Constructs an instance containing the given information.
             * \param seq_fw the forward strand DNA sequence in 5' to 3' order.
             * \param seq_rv the reverse strand DNA sequence in 5' to 3' order.
             * \param ipd_fw the forward strand IPD signal in 5' to 3' order.
             * \param ipd_rv the reverse strand IPD signal in 5' to 3' order.
             * \param pwd_fw the forward strand PWD signal in 5' to 3' order.
             * \param pwd_rv the reverse strand PWD signal in 5' to 3' order.
             */
            KineticSignal(const std::string& seq_fw,
                       const std::string& seq_rv,
                       const std::vector<double>& ipd_fw,
                       const std::vector<double>& ipd_rv,
                       const std::vector<double>& pwd_fw,
                       const std::vector<double>& pwd_rv) ;

            /*!
             * \brief Constructs an instance containing the given information.
             * \param seq_fw the forward strand DNA sequence in 5' to 3' order.
             * \param seq_rv the reverse strand DNA sequence in 5' to 3' order.
             * \param ipd_fw the forward strand IPD signal in 5' to 3' order.
             * \param ipd_rv the reverse strand IPD signal in 5' to 3' order.
             * \param pwd_fw the forward strand PWD signal in 5' to 3' order.
             * \param pwd_rv the reverse strand PWD signal in 5' to 3' order.
             */
            KineticSignal(std::string&& seq_fw,
                       std::string&& seq_rv,
                       std::vector<double>&& ipd_fw,
                       std::vector<double>&& ipd_rv,
                       std::vector<double>&& pwd_fw,
                       std::vector<double>&& pwd_rv) ;
            
            /*!
             * \brief Copy constructor.
             * \param other the instance to copy from.
             */
            KineticSignal(const ngsai::KineticSignal& other) ;

             /*!
             * \brief Move constructor.
             * \param other the instance to move from.
             */
            KineticSignal(ngsai::KineticSignal&& other) ;

            /*!
             * \brief Destructor.
             */
            ~KineticSignal() ;

            /*!
             * \brief Assignment operator.
             * \param other a reference to the object to copy from.
             * \returns a reference to the current instance.
             */
            KineticSignal&
            operator = (const KineticSignal& other) ;

            /*!
             * \brief Move assignment operator.
             * \param other a rvalue to copy from.
             * \returns a reference to the current instance.
             */
            KineticSignal&
            operator = (KineticSignal&& other) ;

            /*!
             * \brief Sets the forward strand DNA sequence.
             * \param seq the new sequence in 5' to 3' order.
             */
            void
            setSequenceFw(const std::string& seq) ;

            /*!
             * \brief Sets the forward strand DNA sequence.
             * \param seq the new sequence in 5' to 3' order.
             */
            void
            setSequenceFw(std::string&& seq) ;

            /*!
             * \brief Sets the reverse strand DNA sequence.
             * \param seq the new sequence in 5' to 3' order.
             */
            void
            setSequenceRv(const std::string& seq) ;

            /*!
             * \brief Sets the reverse strand DNA sequence.
             * \param seq the new sequence in 5' to 3' order.
             */
            void
            setSequenceRv(std::string&& seq) ;

            /*!
             * \brief Sets the forward strand IPDs.
             * \param seq the new IPDs in 5' to 3' order.
             */
            void
            setIPDFw(const std::vector<double>& ipd) ;

            /*!
             * \brief Sets the forward strand IPDs.
             * \param seq the new IPDs in 5' to 3' order.
             */
            void
            setIPDFw(std::vector<double>&& ipd) ;

            /*!
             * \brief Sets the reverse strand IPDs.
             * \param seq the new IPDs in 5' to 3' order.
             */
            void
            setIPDRv(const std::vector<double>& ipd) ;

            /*!
             * \brief Sets the reverse strand IPDs.
             * \param seq the new IPDs in 5' to 3' order.
             */
            void
            setIPDRv(std::vector<double>&& ipd) ;
            
            /*!
             * \brief Sets the forward strand IPDs.
             * \param seq the new IPDs in 5' to 3' order.
             */
            void
            setPWDFw(const std::vector<double>& pwd) ;

            /*!
             * \brief Sets the forward strand IPDs.
             * \param seq the new IPDs in 5' to 3' order.
             */
            void 
            setPWDFw(std::vector<double>&& pwd) ;

            /*!
             * \brief Sets the reverse strand IPDs.
             * \param seq the new IPDs in 5' to 3' order.
             */
            void
            setPWDRv(const std::vector<double>& pwd) ;

            /*!
             * \brief Sets the reverse strand IPDs.
             * \param seq the new IPDs in 5' to 3' order.
             */
            void
            setPWDRv(std::vector<double>&& pwd) ;

            /*!
             * \brief Accesses the forward strand DNA sequence.
             * \returns the forward strand DNA sequence in 5' to 3' order.
             */
            const std::string&
            getSequenceFw() const ;

            /*!
             * \brief Accesses the reverse strand DNA sequence.
             * \returns the reverse strand DNA sequence in 5' to 3' order.
             */
            const std::string&
            getSequenceRv() const ;

            /*!
             * \brief Accesses the forward strand IDPs.
             * \returns the forward strand IPDs in 5' to 3' order.
             */
            const std::vector<double>&
            getIPDFw() const ;

            /*!
             * \brief Accesses the reverse strand IDPs.
             * \returns the reverse strand IPDs in 5' to 3' order.
             */
            const std::vector<double>&
            getIPDRv() const ;

            /*!
             * \brief Accesses the forward strand PWDs.
             * \returns the forward strand PWDs in 5' to 3' order.
             */
            const std::vector<double>&
            getPWDFw() const ;

            /*!
             * \brief Accesses the reverse strand PWDs.
             * \returns the reverse strand PWDs in 5' to 3' order.
             */
            const std::vector<double>&
            getPWDRv() const ;

            /*!
             * \brief Computes the mean IPD values across both strands. If 
             * IPDs are set for only one strand, the average are equal to 
             * this strand values. If no IPD are set on any strands, an empty 
             * vector is returned.
             * \returns the mean IPD values in 5' to 3' order. The vector 
             * is empty if no IPD have been set on any strand.
             */
            std::vector<double>
            getIPDMean() const ;

            /*!
             * \brief Computes the mean PWD values across both strands. If 
             * PWDs are set for only one strand, the average are equal to 
             * this strand values. If no PWD are set on any strands, an empty 
             * vector is returned.
             * \returns the mean PWD values in 5' to 3' order. The vector 
             * is empty if no IPD have been set on any strand.
             */
            std::vector<double>
            getPWDMean() const ;

            /*!
             * \brief Indicates if the forward strand has the sequence, IPD and 
             * PWD values set.
             * \return Whether the forward strand has all data set.
             */
            bool
            isFwComplete() const ;

            /*!
             * \brief Indicates if the reverse strand has the sequence, IPD and 
             * PWD values set.
             * \return Whether the reverse strand has all data set.
             */
            bool
            isRvComplete() const ;

            /*!
             * \brief Indicates if the both strands have the sequence, IPD and 
             * PWD values set.
             * \return Whether both strands have all data set.
             */
            bool
            isComplete() const ; 

            /*!
             * \brief Returns the size of the sequence and signal stored on 
             * a single strand, in bp.
             * \returns the size of the sequence and signal stored. 
             */
            size_t
            size() const ;

        
        protected:
            /*!
             * \brief The size of the sequence and signal stored 
             * on a single strand.
             */
            size_t m_size ;
            /*!
             * \brief The forward strand DNA sequence.
             */
            std::string m_sequence_fw ;
            /*!
             * \brief The reverse strand DNA sequence.
             */
            std::string m_sequence_rv ;
            /*!
             * \brief The forward strand IPD signal.
             */
            std::vector<double> m_ipd_fw ;
            /*!
             * \brief The reverse strand IPD signal.
             */
            std::vector<double> m_ipd_rv ;
            /*!
             * \brief The forward strand PWD signal.
             */
            std::vector<double> m_pwd_fw ;
            /*!
             * \brief The reverse strand PWD signal.
             */
            std::vector<double> m_pwd_rv ;
            /*!
             * \brief Whether a value was set for m_seq_fw.
             */
            bool m_seq_fw_set ;
            /*!
             * \brief Whether a value was set for m_seq_rv.
             */
            bool m_seq_rv_set ;
            /*!
             * \brief Whether a value was set for m_ipd_fw.
             */
            bool m_ipd_fw_set ;
            /*!
             * \brief Whether a value was set for m_ipd_rv.
             */
            bool m_ipd_rv_set ;
            /*!
             * \brief Whether a value was set for m_pwd_fw.
             */
            bool m_pwd_fw_set ;
            /*!
             * \brief Whether a value was set for m_pwd_rv.
             */
            bool m_pwd_rv_set ;

    } ; // KineticSignal


} // namespace ngsai


#endif // NGSAI_KINETICSIGNAL_HPP
