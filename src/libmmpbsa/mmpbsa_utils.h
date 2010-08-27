/* 
 * Utility functions and type/constant definitions for the MMPBSA library.
 *
 * Created by David Coss <David.Coss@stjude.org> 2010
 */


#ifndef MMPBSA_UTILS_H
#define	MMPBSA_UTILS_H

#include <fstream>
#include <iomanip>
#include <cmath>
#include <cctype>
#include <vector>
#include <valarray>
#include <map>
#include <algorithm>

#include "mmpbsa_exceptions.h"
#include "StringTokenizer.h"

//MEAD
#include "MEAD/Coord.h"

#define MMPBSA_XML_TITLE "grid_queue"
#define MMPBSA_PI 3.14159265358979323846
#define MMPBSA_FORMAT std::setprecision(5)
#define DEFAULT_SCNB 2.0
#define DEFAULT_SCEE 1.2
#define DEFAULT_DIELC 1.0




typedef double mmpbsa_t;


namespace mmpbsa_utils {

    /**
     * When a parameter has a list of variables (separated by commas), this method
     * tokenizes that list and loades it into the array, in left to right order.
     *
     * @param values
     * @param array
     */
    int loadListArg(const std::string& values,std::vector<size_t>& array, const size_t& offset = 0);

    /**
     * Return center of interaction beetween two atom coordinate sets.
     * It is defined as the center of the union of the atoms of A that
     * are within CUTOFF of any atoms of B with the atoms of B within
     * CUTOFF of A's atoms.
     * 
     * @param acrds
     * @param bcrds
     * @param cutoff
     * @return
     */
    Coord * interaction_minmax(const std::valarray<mmpbsa_t>& acrds,
            const std::valarray<mmpbsa_t>& bcrds, const mmpbsa_t& cutoff = 4.0);

    /**
     * Performs a lookup of a radius of the given atom. The lookup is not as
     * simple as using the atom name as the key for the provided map, because
     * atom names may contain more than just the symbol for the atom. The method
     * uses a set of naming rules to find the radius.
     * 
     * @param atomName
     * @param radiusMap
     * @return
     */
    mmpbsa_t lookup_radius(const std::string& atomName,
            std::map<std::string,mmpbsa_t>& radiusMap)
            throw (mmpbsa::MMPBSAException);

    std::string toUpperCase(const std::string& bean);

    /**
     * Removes white space at the beginning and end of a string.
     * @param str
     */
    std::string trimString(const std::string& bean);

    /**
     * Determines if the given test value is contained in the array.
     * 
     * @param array
     * @param test
     * @return
     */
    template <class T>
    bool contains(const std::vector<T>& array, const T& test);

    
    /**
     * Returns a new vector containing elements of oldVector which are
     * greater than or equal to the given conditional. 
     * Note:Template class must define operator>=
     *
     * @param oldVector
     * @param conditional
     * @return
     */
    template <class T> std::vector<T> compress_ge(const std::vector<T>& oldVector,
        const T& conditional);
    template <class T> std::vector<T> compress_ge(const std::valarray<T>& oldVector,
            const std::slice& currSlice, const T& conditional);
    template <class T> std::vector<T> compress_ge(const std::valarray<T>& oldVector,
            const T& conditional);

    /**
     * Returns a valarray with elements that are 0 if the previous array's
     * element was non-zero and 1 if the previous array's element was 0
     *
     * @param array
     * @return std::valarray<T>
     */
    template <class T> std::valarray<bool> logical_not(const std::valarray<T>& array)
    {
        return array != T(0);
    }

    /**
     * Returns a vector whose elements are the contents of largerArray,
     * correspoding to the elements whose indices equal the values of subsetIndices.
     * 
     *
     * For example, if subsetIndices = {1,3,5}
     * and largerArray = {323,8,13,42,15,9,23,11}
     * Then take(largerArray,subsetIndices) would return
     * {8,42,9}
     * 
     * @param largerArray
     * @param subsetIndices
     * @return
     */
    template <class T> std::vector<T> take(const std::valarray<T>& largerArray,
        const std::valarray<size_t>& subsetIndices,const std::valarray<bool>& mask);

    /**
     * Returns an array whose elements correspond to the sum of all preceeding
     * elements of the orig array.
     *
     * For example, if orig = {1,2,3,4,5}
     * cumsum(orig) returns {
     * @param orig
     * @return
     */
    template <class T> std::valarray<T> cumsum(const std::valarray<T>& orig);

    template <class T> std::valarray<size_t> cumBoolSum(const std::valarray<bool>& orig);

    /**
     * Returns an array which alternates elements of the two provided arrays
     * 
     * @param left
     * @param right
     * @return std::valarray<T>
     */
    template <class T> std::valarray<T> zip(const std::valarray<T>& left,
            const std::valarray<T>& right);
    template <class T> std::valarray<T> zip(const std::valarray<T>& left, const std::slice& leftSlice,
         const std::valarray<T>& right, const std::slice& rightSlice);
    template <class T> std::valarray<T> zip(const std::vector<T>& left,
            const std::valarray<T>& right);
    template <class T> std::valarray<T> zip(const std::valarray<std::valarray<T> >& toBeZipped);

    /**
     * Returns a valarray copy of the vector with elements clically shifted by n.
     * n > 0 indicates a shift to the left.
     * n < 0 indicates a shift to the right.
     * 
     * @param orig
     * @param n
     * @return std::valarray<T>
     */
    template <class T> std::valarray<T> cshift(const std::vector<T>& orig,
            const int& n);

    /**
     * Returns the index of the first element that equals the reference value.
     * If the value is not found, the size of the array is returned.
     *
     * @param array
     * @param reference
     * @return int
     */
    template <class T> size_t find_first(const std::valarray<T>& array,
            const T& reference);

    template<class T> T* cross_product(const T* A, const T* B,
            const size_t& ndim = 3);

    template<class T> T dot_product(const T* A, const T* B,
            const size_t& ndim = 3);
    
    
}//end mmpbsa_utils namespace

#endif	/* MMPBSA_UTILS_H */

