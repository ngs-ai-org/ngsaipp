#ifndef UNITTESTS_EPIGNETICS_UTILITY_HPP
#define UNITTESTS_EPIGNETICS_UTILITY_HPP

#include <ngsaipp/epigenetics/KmerMap.hpp>


/*!
 * \brief Returns true if both maps have strictly identical field values.
 * \param map1 the 1st map to compare.
 * \param map2 the 2nd map to compare.
 * \return true if both maps have the same field values, false otherwise.
 */
bool
kmerMapEqual(const ngsai::KmerMap& map1,
             const ngsai::KmerMap& map2) ;


/*!
 * \brief Computes the scalar product of a matrix and a scalar.
 * \param m a matrix.
 * \param n a scalar.
 * \return the resulting matrix.
 */
std::vector<std::vector<double>>
operator *(const std::vector<std::vector<double>>& m, double n) ;


/*!
 * \brief Computes the scalar division of a matrix and a scalar.
 * \param m a matrix.
 * \param n a scalar.
 * \return the resulting matrix.
 */
std::vector<std::vector<double>>
operator /(const std::vector<std::vector<double>>& m, double n) ;


/*!
 * \brief Performs a log transformation on a matrix.
 * \param m a matrix.
 * \return the resulting matrix.
 */
std::vector<std::vector<double>>
log(const std::vector<std::vector<double>>& m) ;

#endif  // UNITTESTS_EPIGNETICS_UTILITY_HPP