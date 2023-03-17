#ifndef NGSAI_DNA_CONSTANT_HPP
#define NGSAI_DNA_CONSTANT_HPP

#include <unordered_map>


namespace ngsai {

    namespace dna {
        

        /*!
        * \brief A map containing the reverse complement 
        * of each DNA base.
        */
        const std::unordered_map<char,char> base_pairing({{'A', 'T'},
                                                          {'C', 'G'},
                                                          {'G', 'C'},
                                                          {'T', 'A'},
                                                          {'N', 'N'}}) ;

        /*!
        * \brief A map containing an int code corresponding to each
        * DNA base.
        */
        const std::unordered_map<char,int> base_code({{'A', 0},
                                                      {'C', 1},
                                                      {'G', 2},
                                                      {'T', 3}}) ;

        /*!
        * \brief A map containing the DNA base corresponding to each
        * int code.
        */
        const std::unordered_map<int,char> code_base({{0, 'A'},
                                                      {1, 'C'},
                                                      {2, 'G'},
                                                      {3, 'T'}}) ;

        /*!
        * \brief The size of the DNA character set.
        */
        const size_t dna_alphabet_size = base_code.size() ;

    } // namespace dna
} // namespace ngsai

#endif // NGSAI_DNA_CONSTANT_HPP