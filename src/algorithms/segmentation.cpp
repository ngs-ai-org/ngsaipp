#include <ngsaipp/algorithms/segmentation.hpp>

#include <iostream>
#include <sstream>    // std::ostringstream
#include <list>
#include <vector>
#include <limits>



ngsai::algorithms::Segment::Segment()
    : start(),
      end(),
      l_score(),
      r_score()
{ ; }


ngsai::algorithms::Segment::Segment(size_t start,
                 size_t end,
                 double l_score,
                 double r_score)
    : start(start),
      end(end),
      l_score(l_score),
      r_score(r_score)
{ ; }


ngsai::algorithms::Segment::Segment(const Segment& other)
    : start(other.start),
      end(other.end),
      l_score(other.l_score),
      r_score(other.r_score)
{ ; }


ngsai::algorithms::Segment::~Segment()
{ ; }


ngsai::algorithms::Segment&
ngsai::algorithms::Segment::operator = (
                    const ngsai::algorithms::Segment& other)
{   start   = other.start ;
    end     = other.end ;
    l_score = other.l_score ;
    r_score = other.r_score ;
    return (*this) ;
}


bool
ngsai::algorithms::Segment::operator == (
                    const ngsai::algorithms::Segment& other)
{   if((start   == other.start) and
       (end     == other.end) and
       (l_score == other.l_score) and
       (r_score == other.r_score))
    {   return true ; }
    else
    {   return false ; }
}


bool
ngsai::algorithms::Segment::operator!=(
                    const ngsai::algorithms::Segment& other)
{   return not ((*this) == other) ; }


void
ngsai::algorithms::Segment::extend(
                    const ngsai::algorithms::Segment& other)
{   end = other.end ;
    r_score = other.r_score ;
}


std::string 
ngsai::algorithms::Segment::toString() const
{
    std::ostringstream oss ;
    oss << "<"
        << start << ", "
        << end   << ", "
        << l_score << ", "
        << r_score 
        << ">";
    return oss.str() ;
}


std::ostream& 
ngsai::algorithms::operator << (std::ostream& stream,
                                const ngsai::algorithms::Segment& s)
{   stream << s.toString() ;
    return stream ;
}


std::list<ngsai::algorithms::Segment>::iterator
ngsai::algorithms::find_lscore(
                        std::list<ngsai::algorithms::Segment>& segments, 
                        ngsai::algorithms::Segment& segment)
{
    std::list<ngsai::algorithms::Segment>::reverse_iterator iter = 
                                                            segments.rbegin() ;
    size_t j = segments.size() -1 ;
    while(iter != segments.rend())
    {   if(iter->l_score < segment.l_score)
        {   // std::cout << "lscore found " << j << std::endl ;
            break ;
        }
        j-- ;
        iter++ ;
    }
    return --(iter.base()) ;
}


bool
ngsai::algorithms::update_list(
                            std::list<ngsai::algorithms::Segment>& segments,
                            ngsai::algorithms::Segment& segment)
{   
    bool extended_last = false ;
    // std::cout << segments << std::endl ;

    // search element j with Lj < Lk
    auto iter = ngsai::algorithms::find_lscore(segments, segment) ;
    if(iter == segments.end())
    {   segments.push_back(segment) ; 
        extended_last = false ;
        // std::cout << "appending" << std::endl ;
    }
    else
    {   if(iter->r_score >= segment.r_score)
        {   segments.push_back(segment) ;
            extended_last = false ;
            // std::cout << "appending" << std::endl ;
        }
        else
        {   iter->extend(segment) ;
            while(segments.back() != (*iter))
            {   segments.pop_back() ; }
            extended_last = true ;
            // std::cout << "extension" << std::endl ;
        }
    }
    return extended_last ;
}


std::list<std::pair<size_t,size_t>> 
ngsai::algorithms::find_max_scoring_segments(const std::vector<double>& scores)
{   
    std::list<ngsai::algorithms::Segment> segments ;

    size_t segment_start = 0 ;
    size_t segment_end   = 0 ;
    double segment_lscore = 0. ;
    double segment_rscore = 0. ;
    
    for(size_t i=0; i<scores.size(); i++)
    {   segment_start = i ;
        segment_end   = i + 1 ;
        segment_rscore += scores[segment_end-1] ;

        // positive score segment starts here
        if(scores[i] > 0.)
        {   

            // search end of positive score segment
            while((i<scores.size()-1) and (scores[i+1] > 0.))
            {   segment_end++ ;
                i++ ;
                segment_rscore += scores[segment_end-1] ;
            }

            // candidate segment
            Segment candidate(segment_start,
                              segment_end,
                              segment_lscore,
                              segment_rscore) ;
            // std::cout << "candidate " << candidate.toString() << std::endl ;
            // bool extended_last = update_list(segments, candidate) ;
            // std::cout << extended_last << std::endl ;
            // std::cout << segments << std::endl ;
            // std::cout << std::string(40, '-') << std::endl << std::endl ;
            
            while(ngsai::algorithms::update_list(segments, candidate))
            {   // std::cout << segments << std::endl ;
                // std::cout << std::string(40, '-') << std::endl << std::endl ;
                candidate = segments.back() ; 
                // std::cout << "new candidate " << candidate.toString() << std::endl ;
                segments.pop_back() ;
            }
            // std::cout << segments << std::endl ;
            // std::cout << std::string(40, '-') << std::endl << std::endl ;
            
            
            // next position L-score is previous segment R-score
            segment_lscore = segment_rscore ;
        }
        else
        {   // this is next position L-score 
            segment_lscore += scores[i] ; 
        }
    }

    // convert to pairs
    std::list<std::pair<size_t,size_t>> l_segments ;
    while(segments.size())
    {   l_segments.push_back(std::make_pair(segments.front().start, 
                                            segments.front().end)) ;
        segments.pop_front() ;
    }

    return l_segments ;
}


std::list<std::pair<size_t,size_t>> 
ngsai::algorithms::find_min_length_cover(const std::vector<double>& scores,
                                         double alpha,
                                         size_t min_length_low,
                                         size_t min_length_high)
{   
    size_t n = scores.size() ;

    typedef enum state {low=0,           // low scoring segment state
                        high=1} state ;  // high scoring segment state 
    
    std::vector<state>   z(n) ;   // score label indicating low/high segment
    std::vector<size_t> t0(n) ; 
    std::vector<size_t> t1(n) ; 
    std::vector<double> u0(n) ;
    std::vector<double> u1(n) ;
    std::vector<double> v0(n) ;
    std::vector<double> v1(n) ;

    double m_inf = std::numeric_limits<double>::min() ;

    u0[0] = 0. ; u1[0] = scores[0] - alpha ;

    if(min_length_low > 1)
    {   v0[0] = m_inf ; }
    else
    {   v0[0] = u0[0] ; }

    if(min_length_high > 1)
    {   v1[0] = m_inf ; }
    else
    {   v1[0] = u1[0] ; }

    double s = 0. ;
    for(size_t i=1; i<n; i++)
    {   s += scores[i] ;

        if(i >= (min_length_high - 1))
        {   s -= scores[i-min_length_high+1] ; }

        if(u0[i-1] > v1[i-1])
        {   u0[i] = u0[i-1] ; t0[i] = 0 ; }
        else
        {   u0[i] = v1[i-1] ; t0[i] = 1 ; }

        if((v0[i-1]-alpha) > u1[i-1])
        {   u1[i] = scores[i] + v0[i-1] - alpha ; t1[i] = 0 ; }
        else
        {   u1[i] = scores[i] + u1[i-1] ; t1[i] = 1 ; }

        if(i >= (min_length_low - 1))
        {   v0[i] = u0[i-min_length_low +1] ; }
        else
        {   v0[i] = m_inf ; }

        if(i >= (min_length_high - 1))
        {   v1[i] = u1[i-min_length_high+1] + s ; }
        else
        {   v1[i] = m_inf ; }
    }

    size_t m ;
    if(v0[n-1] > v1[n-1])
    {   z[n-1] = state::low ; m = min_length_low ; }
    else
    {   z[n-1] = state::high ; m = min_length_high ; }

    // i is n-2, n-3, ..., 0
    for(size_t i=n-1; i-->0;)
    {   if(m > 1)
        {   z[i] = z[i+1] ; m -= 1 ; }
        else if(m == 1 and z[i+1] == state::low)
        {   if(t0[i+1] == 0)
            {   z[i] = state::low ; }
            else
            {   z[i] = state::high ; m = min_length_high ; }
        }
        else if((m == 1) and (z[i+1] == state::high))
        {   if(t1[i+1] == 1)
            {   z[i] = high ; }
            else
            {   z[i] = low ; m = min_length_low ; }
        }
    }

    // transform z into segments
    std::list<std::pair<size_t,size_t>> segments ;
    size_t segment_start = 0 ;
    size_t segment_end   = 0 ;
    state  segment_state = state::low ;
    for(size_t i=0; i<z.size(); i++)
    {   // start of high segment
        if((segment_state == state::low) and 
           (z[i] == state::high))
        {   segment_start = i ;
            segment_state = state::high ;
        }
        // end of high segment
        if((segment_state == state::high) and 
           (z[i] == state::low))
        {   segment_end = i ;
            segment_state = state::low ;
            segments.push_back(std::make_pair(segment_start,
                                              segment_end)) ;
        }
    }
    // end of last segment
    if(segment_state == state::high)
    {   segment_end = z.size() ;
        segments.push_back(std::make_pair(segment_start,
                                          segment_end)) ;
    }

    return segments ;
}