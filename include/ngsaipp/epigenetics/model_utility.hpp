#ifndef NGSAI_MODEL_UTILITY_HPP
#define NGSAI_MODEL_UTILITY_HPP

#include <string>
#include <vector>
#include <utility>


#include <ngsaipp/epigenetics/KineticModel.hpp>
#include <ngsaipp/epigenetics/KmerMap.hpp>
#include <ngsaipp/io/BedRecord.hpp>


namespace ngsai
{
    /*!
     * \brief Normalizes the raw kinetic signal by the mean signal 
     * expected  from its sequence context. Note that the normalization 
     * process can only produce a value for the central values of 
     * the vectors. The first K/2 and last K/2 values - where K is the size 
     * of the kmer inside the model - will be dropped.
     * \param sequence the sequence from which the kinetics were taken,
     * \param ipd the raw IPDs corresponding to each base of the sequence.
     * \param pwd the raw PWDs corresponding to each base of the sequence.
     * \param model a map containing the mean expected signal for each 
     * possible kmer.
     * \throws std::invalid_argument if the kinetics and sequence size 
     * don't match.
     * \returns a pair containing the normalized IPD (first) and PWD 
     * (second) values.
     */
    std::pair<std::vector<double>, std::vector<double>>
    normalize_kinetics(const std::string& seq,
                       const std::vector<double>& ipd,
                       const std::vector<double>& pwd,
                       const ngsai::KmerMap& model) ;

    /*!
     * \brief Normalizes the raw kinetic signal by the mean signal 
     * expected  from its sequence context. Note that the normalization 
     * process can only produce a value for the central values of 
     * the vectors. The first K/2 and last K/2 values - where K is the size 
     * of the kmer inside the model - will be dropped.
     * \param sequence the sequence from which the kinetics were taken,
     * \param ipd the raw IPDs corresponding to each base of the sequence.
     * \param pwd the raw PWDs corresponding to each base of the sequence.
     * \param model a map containing the mean expected signal for each 
     * possible kmer.
     * \throws std::invalid_argument if the kinetics and sequence size 
     * don't match.
     * \returns a pair containing the normalized IPD (first) and PWD 
     * (second) values.
     */
    std::pair<std::vector<double>, std::vector<double>>
    normalize_kinetics(const std::string& seq,
                       const std::vector<uint16_t>& ipd,
                       const std::vector<uint16_t>& pwd,
                       const ngsai::KmerMap& model) ;

    /*!
     * \brief Trains a KineticModel, from the CCS overlapping the genomic 
     * windows listed as BedRecord. The model is trained on the subset of 
     * regions [regions_from,region_to).
     * For each CpG listed, the mean IPD and PWD values at each position within 
     * the window (defined in model) are computed from all CCS aligning at 
     * this position. The mean IPD and PWD values are then introduced in the 
     * model.
     * \param model a pointer to the kinetic model to train. This model must 
     * have been already initialised using setParameters().
     * \param region a list of genomic regions containing the CpG 
     * coordinates from which the signal must be learnt.
     * \param regions_from the index of the 1st region to use for training in 
     * the vector of regions.
     * \param regions the index of the past the last region to use for 
     * training in the vector of regions.
     * \param paths_bam the paths to the bam files containing the  
     * mapped CCS from which the signal must be extracted to learn the 
     * models.
     * \throws std::invalid_argument if the model was not yet initialised.
     */
    void
    train_KineticModel(ngsai::KineticModel* model,
                       const std::vector<ngsai::BedRecord>& regions,
                       size_t regions_from,
                       size_t regions_to,
                       const std::vector<std::string>& paths_bam) ;
    
    /*!
     * \brief Trains a KineticModel, from the CCS overlapping the genomic 
     * windows listed as BedRecord.
     * For each CpG listed, the mean IPD and PWD values at each position within 
     * the window (defined in model) are computed from all CCS aligning at 
     * this position. The mean IPD and PWD values are then introduced in the 
     * model.
     * \param model a pointer to the kinetic model to train. This model must 
     * have been already initialised using setParameters().
     * \param region a list of genomic regions containing the CpG 
     * coordinates from which the signal must be learnt.
     * \param paths_bam the paths to the bam files containing the  
     * mapped CCS from which the signal must be extracted to learn the 
     * models.
     * \throws std::invalid_argument if the model was not yet initialised.
     */
    void
    train_KineticModel(ngsai::KineticModel* model,
                       const std::vector<ngsai::BedRecord>& regions,
                       const std::vector<std::string>& paths_bam) ;
    

    /*!
     * \brief Trains a KineticModel, from the CCS overlapping the genomic 
     * windows listed in the BED file.
     * For each CpG listed, the mean IPD and PWD values at each position within 
     * the window (defined in model) are computed from all CCS aligning at 
     * this position. The mean IPD and PWD values are then introduced in the 
     * model.
     * \param model a pointer to the kinetic model to train. This model must 
     * have been already initialised using setParameters().
     * \param path_bed the path to the bed file containing the CpG 
     * coordinates from which the signal must be learnt.
     * \param paths_bam the paths to the bam files containing the  
     * mapped CCS from which the signal must be extracted to learn the 
     * models.
     * \throws std::invalid_argument if the model was not yet initialised.
     */
    void
    train_KineticModel(ngsai::KineticModel* model,
                       const std::string& path_bed,
                       const std::vector<std::string>& paths_bam) ;

}  // namespace ngsai

#endif // NGSAI_MODEL_UTILITY_HPP