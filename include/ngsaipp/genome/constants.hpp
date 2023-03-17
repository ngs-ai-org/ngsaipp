#ifndef NGSAI_GENOME_CONSTANT_HPP
#define NGSAI_GENOME_CONSTANT_HPP

namespace ngsai {

    namespace genome {

        /*!
         * \brief Defines a strand on a DNA molecule.
         */
        enum strand {FORWARD=0,
                     REVERSE=1,
                     UNORIENTED=2} ;

        /*!
         * \brief Converts a char encoded strand into a ngsai::genome::strand 
         * value.
         * \param strand the strand encoded as a char, either '+', '-' or '.'
         * \return the corresponding ngsai::genome::strand value.
         * \throws std::invalid_argument if the char code corresponds to 
         * nothing.
         */
        ngsai::genome::strand char_to_strand(char strand) ;


        /*!
         * \brief Converts a ngsai::genome::strand value into the corresponding 
         * char value.
         * \param strand the strand encoded as ngsai::genome::strand.
         * \return the corresponding char strand value.
         * \throws std::invalid_argument if the char code corresponds to 
         * nothing.
         */
        char strand_to_char(ngsai::genome::strand strand) ;
    }  // namespace genome
}  // namespace ngsai


#endif // NGSAI_GENOME_CONSTANT_HPP