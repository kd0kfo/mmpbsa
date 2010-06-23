#include "mmpbsa_utils.h"

template <class T> std::vector<T> mmpbsa_utils::compress_ge(const std::vector<T>& oldVector,
    const T& conditional)
{
    std::vector<T> returnMe;

    for(int i = 0;i<oldVector.size();i++)
        if(oldVector[i] >= conditional)
            returnMe.push_back(oldVector[i]);

    return returnMe;
}

template <class T> std::vector<T> mmpbsa_utils::compress_ge(const std::valarray<T>& oldVector,
    const T& conditional)
{
    std::vector<T> returnMe;

    for(int i = 0;i<oldVector.size();i++)
        if(oldVector[i] >= conditional)
            returnMe.push_back(oldVector[i]);

    return returnMe;
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

    for(int i = 1;i<size;i++)
        returnMe[i] = orig[slice(0,i,1)].sum();

    return returnMe;
}

template <class T> std::valarray<T> mmpbsa_utils::zip(const std::valarray<T>& left,
            const std::valarray<T>& right)
{
    using std::slice;
    using std::slice_array;
    if(left.size() != right.size())
        throw MMPBSAException("When using zip, the two arrays must have the same size.",DATA_FORMAT_ERROR);

    int oldsize = left.size();
    std::valarray<T> returnMe(2*oldsize);
    for(int i = 0;i<oldsize;i++)
    {
        slice_array<T> curr = returnMe[slice(i,2,1)];
        curr[0] = left[i];
        curr[1] = right[i];
    }
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

template <class T> int find_first(const std::valarray<T>& array,
            const T& reference)
{
    int half_size = int(array.size()/2);
    for(int i = 0;i<half_size;i++)
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