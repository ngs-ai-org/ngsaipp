#include <ngsaipp/io/BedReader.hpp>
#include <ngsaipp/io/BedRecord.hpp>

#include <string>
#include <vector>
#include <iostream>
#include <stdexcept>

#include <ngsaipp/genome/constants.hpp>         // ngsai::genome::char_to_strand()
#include <ngsaipp/utility/string_utility.hpp>   // ngsai::split()


ngsai::BedReader::BedReader(const std::string& path)
    : m_path(path),
      m_f_bed(),
      m_is_open(false)
{   this->open() ; }


ngsai::BedReader::~BedReader()
{   if(m_is_open)
    {   this->close() ; }
}


bool 
ngsai::BedReader::getNext(ngsai::BedRecord& record)
{   if(not m_is_open)
    {   char msg[4096] ;
        sprintf(msg, 
                "BedReader error : cannot read from %s, file is closed", 
                m_path.c_str()) ;
        throw std::runtime_error(msg) ;
    }
    std::string line ;
    std::getline(m_f_bed, line) ;
    // EOF?
    if(m_f_bed.eof())
    {   return false ; }

    // tokenize
    std::vector<std::string> fields = ngsai::split(line, '\t') ;
    size_t n_fields = fields.size() ;
    if(n_fields != 6)
    {   char msg[4096] ; 
        sprintf(msg,
                "BedReader error : bed entry must have 6 fields (%zu)",
                n_fields) ;
        throw std::runtime_error(msg) ;
    }

    // convert
    try
    {
        record = 
            std::move(
                ngsai::BedRecord(fields[0],
                                 std::stoi(fields[1]),
                                 std::stoi(fields[2]),
                                 fields[3],
                                 std::stod(fields[4]),
                                 ngsai::genome::char_to_strand(fields[5][0]))) ;
    }
    catch(const std::invalid_argument& e)
    {   throw std::runtime_error(e.what()) ; }

    // record.chrom  = fields[0] ;
    // record.start  = std::stoi(fields[1]) ;
    // record.end    = std::stoi(fields[2]) ;
    // record.name   = fields[3] ;
    // record.score  = std::stoi(fields[4]) ;
    // record.strand = ngsai::genome::char_to_strand(fields[5][0]) ;

    return true ;
}


void
ngsai::BedReader::close()
{   if(m_is_open)
    {   m_f_bed.close() ;
        m_is_open = false ;
    }
}


void
ngsai::BedReader::open()
{   m_f_bed.open(m_path, std::ios::in) ;
    if(m_f_bed.fail())
    {   char msg[4096] ;
        sprintf(msg,
                "BedReader error : could not open %s", 
                m_path.c_str()) ;
        throw std::runtime_error(msg) ;
    }
    m_is_open = true ;
}