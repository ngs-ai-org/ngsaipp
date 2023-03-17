#include <ngsaipp/epigenetics/model_utility.hpp>
#include <ngsaipp/epigenetics/KineticModel.hpp>
#include <ngsaipp/epigenetics/KmerMap.hpp>
#include <ngsaipp/epigenetics/KineticSignal.hpp>
#include <ngsaipp/epigenetics/CcsKineticExtractor.hpp>
#include <ngsaipp/io/bed_io.hpp>

#include <vector>
#include <string>
#include <stdexcept>
#include <pbbam/CompositeBamReader.h>  // GenomicIntervalCompositeBamReader
#include <pbbam/BamRecord.h>      



std::pair<std::vector<double>, std::vector<double>>
ngsai::normalize_kinetics(const std::string& seq,
                          const std::vector<double>& ipd,
                          const std::vector<double>& pwd,
                          const ngsai::KmerMap& model)
{   
    size_t l      = ipd.size() ;
    size_t k      = model.getKmerSize() ;
    size_t center = k / 2 ;
    size_t l_norm = l - 2*center ;

    // iterator to all models
    auto iters = model.find(seq) ;

    // normalized kinetic values
    // std::vector<double>  ipdr(l, 0.) ;
    // std::vector<double>  pwdr(l, 0.) ;
    std::vector<double>  ipdr(l_norm, 0.) ;
    std::vector<double>  pwdr(l_norm, 0.) ;
    
    // sum signal
    for(size_t start=0; start<l-k+1; start++)
    {   // model corresponding to subsequence starting at start
        auto iter = iters[start] ;
        // the model does not contain values for this sequence
        if(iter == model.end())
        {   continue ; }
        
        // ipdr[start+center] += static_cast<double>(ipd[start+center]) / 
        //                       static_cast<double>(iter->second.ipd[center]) ;
        // pwdr[start+center] += static_cast<double>(pwd[start+center]) / 
        //                       static_cast<double>(iter->second.pwd[center]) ;
        ipdr[start] += static_cast<double>(ipd[start+center]) / 
                       static_cast<double>(iter->second.ipd[center]) ;
        pwdr[start] += static_cast<double>(pwd[start+center]) / 
                       static_cast<double>(iter->second.pwd[center]) ;
    }
    
    return std::make_pair(ipdr, pwdr) ;
}


std::pair<std::vector<double>, std::vector<double>>
ngsai::normalize_kinetics(const std::string& seq,
                          const std::vector<uint16_t>& ipd,
                          const std::vector<uint16_t>& pwd,
                          const ngsai::KmerMap& model)
{
    std::vector<double> ipd_d(ipd.begin(), ipd.end()) ;
    std::vector<double> pwd_d(pwd.begin(), pwd.end()) ;
    return ngsai::normalize_kinetics(seq, ipd_d, pwd_d, model) ;
}


void
ngsai::train_KineticModel(ngsai::KineticModel* model,
                          const std::vector<ngsai::BedRecord>& regions,
                          size_t regions_from,
                          size_t regions_to,
                          const std::vector<std::string>& paths_bam)
{   
    if(not model->isInit())
    {   throw std::invalid_argument("training error : KineticModel was not " 
                                    "initialized") ;
    }
    if(regions_from >= regions_to)
    {   char msg[4096] ;
        sprintf(msg,
                "training error : from (%zu) must be < to (%zu)",
                regions_from,
                regions_to) ;
        throw std::invalid_argument(msg) ;
    }
    if(regions_to > regions.size())
    {   char msg[4096] ;
        sprintf(msg,
                "training error : from (%zu) must be < to (%zu)",
                regions_from,
                regions_to) ;
        throw std::invalid_argument(msg) ;
    }

    // size_t window_size = model->size() ;
    size_t window_size = model->getKineticSignalRequiredSize() ;

    // half the window size
    size_t win_size_half = window_size / 2 ;

    // bam readers
    PacBio::BAM::GenomicIntervalCompositeBamReader reader_bam(paths_bam) ;

    // extract kinetics at all regions and fill histograms
    ngsai::BedRecord cpg ;
    ngsai::CcsKineticExtractor extractor ;
    PacBio::BAM::BamRecord ccs ;
    for(size_t i= regions_from; i<regions_to; i++)
    {   cpg = regions[i] ;
        // only keep CpG on + and compute coordinates of + and - strand CpGs
        if(cpg.strand == ngsai::genome::REVERSE)
        {   continue ; }

        // CpG window on + strand
        ngsai::BedRecord window_p(cpg) ;
        window_p.start -= win_size_half ;
        window_p.end   += win_size_half - 1 ;
        // CpG window on - strand
        ngsai::BedRecord window_m(cpg) ;
        window_m.strand = ngsai::genome::REVERSE ;
        window_m.start -= win_size_half - 1 ;
        window_m.end   += win_size_half ;

        PacBio::BAM::GenomicInterval interval(cpg.chrom, 
                                              cpg.start,
                                              cpg.end) ;

        // compute mean IPD and PWD, on each strand, from CCS for this region
        double n_fw = 0. ;
        double n_rv = 0. ;
        std::string seq_fw ;
        std::string seq_rv ;
        std::vector<double> ipds_fw_m(window_size, 0.) ;
        std::vector<double> ipds_rv_m(window_size, 0.) ;
        std::vector<double> pwds_fw_m(window_size, 0.) ;
        std::vector<double> pwds_rv_m(window_size, 0.) ;
        reader_bam.Interval(interval) ;
        while(reader_bam.GetNext(ccs))
        {   // extract kinetics on + strand
            if(extractor.extract(ccs, window_p))
            {   if(n_fw == 0.)
                {   seq_fw = extractor.getSequence() ; }
                std::vector<uint16_t> ipds = extractor.getIPD() ;
                std::vector<uint16_t> pwds = extractor.getPWD() ;
                for(size_t i=0; i<ipds_fw_m.size(); i++)
                {   ipds_fw_m[i] += static_cast<double>(ipds[i]) ;
                    pwds_fw_m[i] += static_cast<double>(pwds[i]) ;
                }
                n_fw += 1. ;
            }
            // extract kinetics on - strand
            if(extractor.extract(ccs, window_m))
            {   if(n_rv == 0.)
                {   seq_rv = extractor.getSequence() ; }
                std::vector<uint16_t> ipds = extractor.getIPD() ;
                std::vector<uint16_t> pwds = extractor.getPWD() ;
                for(size_t i=0; i<ipds_rv_m.size(); i++)
                {   ipds_rv_m[i] += static_cast<double>(ipds[i]) ;
                    pwds_rv_m[i] += static_cast<double>(pwds[i]) ;
                }
                n_rv += 1. ;
            }
        }
        
        ngsai::KineticSignal kinetics ;
        if(n_fw != 0.)
        {   // compute mean
            for(size_t i=0; i<window_size; i++)
            {   ipds_fw_m[i] /= n_fw ;
                pwds_fw_m[i] /= n_fw ;
            }
            // insert
            kinetics.setSequenceFw(std::move(seq_fw)) ;
            kinetics.setIPDFw(std::move(ipds_fw_m)) ;
            kinetics.setPWDFw(std::move(pwds_fw_m)) ;
        }
        if(n_rv != 0.)
        {   // compute mean
            for(size_t i=0; i<window_size; i++)
            {   ipds_rv_m[i] /= n_rv ;
                pwds_rv_m[i] /= n_rv ;
            }
            // insert
            kinetics.setSequenceRv(std::move(seq_rv)) ;
            kinetics.setIPDRv(std::move(ipds_rv_m)) ;
            kinetics.setPWDRv(std::move(pwds_rv_m)) ;
        }
        if(kinetics.isComplete())
        {   model->add(kinetics) ; }
    }

}


void
ngsai::train_KineticModel(ngsai::KineticModel* model,
                          const std::vector<ngsai::BedRecord>& regions,
                          const std::vector<std::string>& paths_bam)
{   ngsai::train_KineticModel(model, regions, 0, regions.size(), paths_bam) ; }


void
ngsai::train_KineticModel(ngsai::KineticModel* model,
                          const std::string& path_bed,
                          const std::vector<std::string>& paths_bam)
{   
    if(not model->isInit())
    {   throw std::invalid_argument("training error : KineticModel was not " 
                                    "initialized") ;
    }

    // bam readers
    PacBio::BAM::GenomicIntervalCompositeBamReader reader_bam(paths_bam) ;

    // list all CpGs and fill histograms
    std::vector<ngsai::BedRecord> cpgs ;
    ngsai::BedRecord cpg ;
    ngsai::BedReader bed_reader(path_bed) ;
    while(bed_reader.getNext(cpg))
    {   
        // only keep CpG on + and compute coordinates of + and - strand CpGs
        if(cpg.strand == ngsai::genome::REVERSE)
        {   continue ; }

        cpgs.push_back(cpg) ;
    }
    bed_reader.close() ;

    // train model
    ngsai::train_KineticModel(model, cpgs, paths_bam) ;

}
