#include "mmpbsa_utils.h"

template <class T> std::vector<T> mmpbsa_utils::compress_ge(const std::vector<T>& oldVector,
    const T& conditional)
{
    std::vector<T> returnMe;

    for(size_t i = 0;i<oldVector.size();i++)
        if(oldVector[i] >= conditional)
            returnMe.push_back(oldVector[i]);

    return returnMe;
}

template <class T> std::vector<T> mmpbsa_utils::compress_ge(const std::valarray<T>& oldVector,
    const std::slice& currSlice, const T& conditional)
{
    std::vector<T> returnMe;

    for(size_t i = currSlice.start();i<oldVector.size();i += currSlice.stride())
        if(oldVector[i] >= conditional)
            returnMe.push_back(oldVector[i]);

    return returnMe;/**/
}

template <class T> std::valarray<T> mmpbsa_utils::take(const std::valarray<T>& largerArray,
    const std::slice_array<T>& subsetIndices)
{
    std::valarray<T> returnMe = largerArray[subsetIndices];
    return returnMe;
}

template <class T> std::valarray<T> mmpbsa_utils::cumsum(const std::valarray<T>& orig)
{
    using std::slice;

    int size = orig.size();
    std::valarray<T> returnMe(size);

    returnMe[0] = orig[0];

    for(int i = 1;i<size;i++)
    {
        std::valarray<T> curr = orig[slice(0,i+1,1)];
        returnMe[i] = curr.sum();
    }

    return returnMe;
}

template <class T> std::valarray<T> mmpbsa_utils::zip(const std::valarray<T>& left,
            const std::valarray<T>& right)
{
    if(left.size() != right.size())
        throw MMPBSAException("When using zip, the two arrays must have the same size.",DATA_FORMAT_ERROR);

    int oldsize = left.size();
    std::valarray<T> returnMe(2*oldsize);
    size_t returnMeIndex = 0;
    for(size_t i = 0;i<oldsize;i++)
    {
        returnMe[returnMeIndex++] = left[i];
        returnMe[returnMeIndex++] = right[i];
    }
    return returnMe;
}

template <class T> std::valarray<T> mmpbsa_utils::zip(const std::valarray<T>& left, const std::slice& leftSlice,
         const std::valarray<T>& right, const std::slice& rightSlice)
{
    int oldsize = leftSlice.size();
    if(oldsize != rightSlice.size())
        throw MMPBSAException("Cannot zip two slices of difference sizes.",DATA_FORMAT_ERROR);
    
    std::valarray<T> returnMe(2*oldsize);
    size_t returnMeIndex = 0;
    for(size_t i = 0;i<oldsize;i++)
    {
        returnMe[returnMeIndex++] = left[leftSlice.start()+i*leftSlice.stride()];
        returnMe[returnMeIndex++] = right[rightSlice.start()+i*rightSlice.stride()];
    }
    return returnMe;

}

template <class T> std::valarray<T> mmpbsa_utils::zip(const std::vector<T>& left,
            const std::valarray<T>& right)
{
    if(left.size() != right.size())
        throw MMPBSAException("When using zip, the two arrays must have the same size.",DATA_FORMAT_ERROR);

    int oldsize = left.size();
    std::valarray<T> returnMe(2*oldsize);
    size_t returnMeIndex = 0;
    for(size_t i = 0;i<oldsize;i++)
    {
        returnMe[returnMeIndex++] = left[i];
        returnMe[returnMeIndex++] = right[i];
    }
    return returnMe;
}

template <class T> std::valarray<T> mmpbsa_utils::cshift(const std::vector<T>& orig,
            const int& n)
{
    int size = orig.size();
    if(abs(n) > size)
        throw MMPBSAException("cshift amount must be between 0 and the array "
                "size",DATA_FORMAT_ERROR);
    
    std::valarray<T> returnMe(size);
    int returnMeIndex = 0;
    int shift = n;
    if(n < 0)
        shift = size+n;
    for (int i = shift; i < size; i++)
        returnMe[returnMeIndex++] = orig[i];
    for (int i = 0; i < shift; i++)
        returnMe[returnMeIndex++] = orig[i];

    return returnMe;


}

template <class T> size_t mmpbsa_utils::find_first(const std::valarray<T>& array,
            const T& reference)
{
    size_t half_size = size_t(array.size()/2);
    for(size_t i = 0;i<half_size;i++)
    {
        if(array[i] == reference)
            return i;
        if(array[i+half_size] == reference)
            return i+half_size;
    }
    if(array[array.size()-1] == reference)//in case of an odd sized array.
        return array.size()-1;
    else
        return array.size();//reference not found.

}
