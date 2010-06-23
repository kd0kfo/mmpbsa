/* 
 * File:   mmpbsa_utils.h
 * Author: dcoss
 *
 * Created on June 23, 2010, 9:08 AM
 */


#ifndef MMPBSA_UTILS_H
#define	MMPBSA_UTILS_H

typedef double mmpbsa_t;

#include <vector>
#include <valarray>

namespace mmpbsa_utils {

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
        return array != 0;
    }

    /**
     * Returns a valarray whos elements are the contents of largerArray,
     * correspoding to the elements whose indices equal the values of subsetIndices.
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
    template <class T> std::valarray<T> take(const std::valarray<T>& largerArray,
    const std::slice_array<T>& subsetIndices); //perhaps use indirect_array


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

    /**
     * Returns an array which alternates elements of the two provided arrays
     * 
     * @param left
     * @param right
     * @return std::valarray<T>
     */
    template <class T> std::valarray<T> zip(const std::valarray<T>& left,
            const std::valarray<T>& right);

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
    template <class T> int find_first(const std::valarray<T>& array,
            const T& reference);

}//end mmpbsa_utils namespace

#endif	/* MMPBSA_UTILS_H */

