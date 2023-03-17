#ifndef NGSAI_ALGORITHMS_SEGMENTATION_HPP
#define NGSAI_ALGORITHMS_SEGMENTATION_HPP

#include <iostream>
#include <string>
#include <vector>
#include <list>


namespace ngsai
{
    namespace algorithms
    {
        /*!
         * The segment class is a data structure for the 
         * find_max_scoring_segments() function implementing Ruzzo and Tompa 
         * 1999's algorithm. 
         * The segment class allows to model a score sub-sequence inside a 
         * score super-sequence, together with the data required by Ruzzo and 
         * Tompa's algorithm. A sub-sequence encompasses the elements 
         * contained in the range [from,to) of the super-sequence.
         */
         class Segment
         {   public:
                 /*!
                 * \brief Constructs and empty instance.
                 */
                 Segment() ;
 
                 /*!
                 * \brief Value constructor.
                 * \param start the index of the 1st element of the 
                 * super-sequence that is part of the sub-sequence.
                 * \param end the index of the past-last element of the 
                 * super-sequence that is part of the sub-sequence.
                 * \param l_score the L-score of this sub-sequence.
                 * \param r_score the R-score of this sub-sequence.
                 */
                 Segment(size_t start,
                         size_t end,
                         double l_score,
                         double r_score) ;
 
                 /*!
                 * \brief Copy constructor.
                 * \param other an instance to copy from.
                 */
                 Segment(const Segment& other) ;
 
                 /*!
                 * \brief Destructor.
                 */
                 ~Segment() ;
 
                 /*!
                 * \brief Assignment operator.
                 * \param other an instance to assign from.
                 * \returns a reference to the current instance.
                 */
                 Segment&
                 operator = (const Segment& other) ;
 
                 /*!
                 * \brief Equality operator.
                 * \param other an instance to compare the current instance 
                 * with.
                 * \returns whether both instances have equal field values.
                 */
                 bool
                 operator == (const Segment& other) ;
 
                 /*!
                 * \brief Inequality operator.
                 * \param other an instance to compare the current instance 
                 * with.
                 * \returns whether both instances have at least one different 
                 * field value.
                 */
                 bool 
                 operator != (const Segment& other) ;
 
                 /*!
                 * \brief Given an other instance modelling a sub-sequence 
                 * located *after* the current instance sub-sequence, this 
                 * sub-sequence is modified to include up to the last element 
                 * in the other sub-sequence.
                 * \param other an other instance modellign a sub-sequence 
                 * located *after* the current instance sub-sequence. If this 
                 * is not the case, the behavior of this method has no sense.
                 */
                 void
                 extend(const Segment& other) ;
 
                 /*!
                 * \brief Creates a string representation of the current 
                 * instance.
                 * \returns the string representation of the current instance.
                 */
                 std::string
                 toString() const ;
 
             public:
                 /*!
                 * \brief the index of the 1st element of the super-sequence 
                 * contained in this sub-sequence.
                 */
                 size_t start ;
                 /*!
                 * \brief the index of the past-last element of the 
                 * super-sequence contained in this sub-sequence.
                 */
                 size_t end ;
                 /*!
                 * \brief the cumulative sum of all scores up to *but not 
                 * including* the leftmost score of this sub-sequence (also 
                 * called L-score).
                 */
                 double l_score ;
                 /*!
                 * \brief the cumulative sum of all scores up to and 
                 * *including* the rightmost score of this sub-sequence 
                 * (also called R-score). 
                 */
                 double r_score ;

        } ;  // class Segment

        /*! \brief Sends a Segment's reprensentation to stream.
         * \param stream a stream of interest. 
         * \param s an instance of Segment of interest. 
         */
        std::ostream& operator << (
                    std::ostream& stream,
                    const ngsai::algorithms::Segment& s) ;

        /*!
         * \brief Routine of find_max_scoring_segments().
         * Searches the segment list from right to left for the first 
         * segment having a L-score smaller than the candidate L-score. This 
         * corresponds to algorithm step 1.
         * \param segments the current list of maximum scoring segments.
         * \param segment a candidate maximum scoring segment.
         * \returns an iterator to the element in segments with a L-score 
         * smaller than the candidate segment L-score. Returns segments.end() 
         * if no such element could be found.
         */
        std::list<Segment>::iterator
        find_lscore(std::list<Segment>& segments, Segment& segment) ;

        /*!
         * \brief Routine of find_max_scoring_segments().
         * Process a candidate maximum scoring segment candidate and 
         * appends it into the list of already known segments or merge it with 
         * a or more already existing ones. This corresponds to steps 1 to 4 
         * in the algorithm desription.
         * \param segments the current list of maximum scoring segments.
         * \param segment a candidate maximum scoring segment.
         * \returns true if the candidate segment was merged with one or more 
         * already existing segments in segments, false if the candidate 
         * segment was appended at the end of segments.
         */
        bool
        update_list(std::list<Segment>& segments, Segment& segment) ;

        /*!
         * \brief Implementation Ruzzo and Tompa 1999's algorithm (A Linear 
         * Time ALgorithm for Finding all Maximal Scoring Subsequences).
         * Finds all maximum positive scoring segments within the given 
         * sequence of scores.
         * \param scores the sequence of positive and negative scores. 
         * Positive scores means good, negative scores means bad.
         * \returns a list of maximum scoring segments. Each pair in the 
         * list indicates a range [from,to) in the score sequence that 
         * delimitates a maximum scoring segment.  Each range can be retrieved 
         * using pair.first and pair.second respectively.
         */
        std::list<std::pair<size_t,size_t>> 
        find_max_scoring_segments(const std::vector<double>& scores) ;

        /*!
         * \brief Implementation of Csũros 2004's MinLength-Cover's algorithm 
         * (Maximum Scoring Segment Sets).
         * Finds a set of maximum positive scoring segments within the given 
         * sequence of scores given the constrains. This set of segments 
         * (k-cover) has the optimal score given the constrains.
         * \param the sequence of positive and negative scores. 
         * Positive scores means good, negative scores means bad.
         * \param alpha the segment set size penaly. The penalized segment set 
         * score is ŵ(S) = w(S) - a*|S| where w(S) is the raw segment set 
         * score and |S| is the number of segments in the set. 
         * \param min_length_low the minimum size of a low scoring segment.
         * \param min_length_high the minimum size of a high scoring segment.
         * \returns a list of maximum scoring segments. Each pair in the 
         * list indicates a range [from,to) in the score sequence that 
         * delimitates a maximum scoring segment.  Each range can be retrieved 
         * using pair.first and pair.second respectively.
         */
        std::list<std::pair<size_t,size_t>> 
        find_min_length_cover(const std::vector<double>& scores,
                              double alpha,
                              size_t min_length_low,
                              size_t min_length_high) ;
        
    } // namespace algorithms
    
} // namespace ngsai



#endif // NGSAI_ALGORITHMS_SEGMENTATION_HPP