#include "mmpbsa_utils.h"


template <class T> bool mmpbsa_utils::contains(const std::vector<T>& array, const T& test)
{
    for(size_t i = 0;i<array.size();i++)
        if(array[i] == test)
            return true;

    return false;
}

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

    return returnMe;
}

template <class T> std::vector<T> mmpbsa_utils::compress_ge(const std::valarray<T>& oldVector,
    const T& conditional)
{
    std::vector<T> returnMe;

    for(size_t i = 0;i<oldVector.size();i++)
        if(oldVector[i] >= conditional)
            returnMe.push_back(oldVector[i]);

    return returnMe;
}

template <class T> std::vector<T> mmpbsa_utils::take(const std::valarray<T>& largerArray,
    const std::valarray<size_t>& subsetIndices,const std::valarray<bool>& mask)
{
    std::vector<T> returnMe;
    for(size_t i = 0;i<mask.size();i++)
        if(mask[i])
        {
            if(subsetIndices[i] > largerArray.size())
            {
                std::ostringstream error;
                error << "An index (" << subsetIndices[i] << ") to be used in "
                        "taking from the larger array exceeds "
                        "the size (" << largerArray.size() << ") of the larger array.";
                throw mmpbsa::MMPBSAException(error,mmpbsa::INVALID_ARRAY_SIZE);
            }
            returnMe.push_back(largerArray[subsetIndices[i]]);
        }

    return returnMe;
}

template <class T> std::valarray<size_t> mmpbsa_utils::cumBoolSum(const std::valarray<bool>& orig)
{
    using std::slice;
    using std::valarray;

    int size = orig.size();
    valarray<size_t> returnMe(1,size);

    returnMe[0] = orig[0];

    for(int i = 1;i<size;i++)
    {
        returnMe[i] = returnMe[i-1]+orig[i];
    }

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
        throw mmpbsa::MMPBSAException("When using zip, the two arrays must have "
                "the same size.",mmpbsa::DATA_FORMAT_ERROR);

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
        throw mmpbsa::MMPBSAException("Cannot zip two slices of difference "
                "sizes.",mmpbsa::DATA_FORMAT_ERROR);

    using std::slice;

    std::valarray<T> returnMe(2*oldsize);
    for(size_t i = 0;i<oldsize;i++)
    {
        returnMe[2*i] = left[leftSlice.start()+i*leftSlice.stride()];
        returnMe[2*i+1] = right[rightSlice.start()+i*rightSlice.stride()];
    }
    return returnMe;

}

template <class T> std::valarray<T> mmpbsa_utils::zip(const std::vector<T>& left,
            const std::valarray<T>& right)
{
    if(left.size() != right.size())
        throw mmpbsa::MMPBSAException("When using zip, the two arrays must have "
                "the same size.",mmpbsa::DATA_FORMAT_ERROR);

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

template <class T> std::valarray<T> mmpbsa_utils::zip(const std::valarray<std::valarray<T> >& toBeZipped)
{
    size_t newSize = toBeZipped[0].size();//elements in each array
    for(size_t i = 1;i<toBeZipped.size();i++)
        if(toBeZipped[i].size() != newSize)
        {
            std::ostringstream error;
            error << "When using mmpbsa_utils::zip, all arrays "
                    "must have the same length. Here the expected length was "
                    << newSize << ", but one array had a length of " << toBeZipped[i].size();
            throw mmpbsa::MMPBSAException(error,mmpbsa::INVALID_ARRAY_SIZE);
        }

    std::valarray<T> returnMe(newSize*toBeZipped.size());
    size_t returnMeIndex = 0;
    for(size_t i = 0;i<newSize;i++)
        for(size_t j = 0;j<toBeZipped.size();j++)
            returnMe[returnMeIndex++] = toBeZipped[j][i];

    return returnMe;
}

template <class T> std::valarray<T> mmpbsa_utils::cshift(const std::vector<T>& orig,
            const int& n)
{
    int size = orig.size();
    if(abs(n) > size)
        throw mmpbsa::MMPBSAException("cshift amount must be between 0 and the array "
                "size",mmpbsa::DATA_FORMAT_ERROR);

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

template<class T> T* mmpbsa_utils::cross_product(const T* A, const T* B,const size_t& dim)
{

    if(dim != 3)
    {
        std::ostringstream error;
        error << "mmpbsa_utils::cross_product currently only works on 3 dimensional vectors."
                "Here, a " << dim << "-dimensional vector was given.";
        throw mmpbsa::MMPBSAException(error,mmpbsa::DATA_FORMAT_ERROR);
    }

    T * cross = new T[dim];

    cross[0] = A[1]*B[2]-A[2]*B[1];
    cross[1] = -A[0]*B[2]+A[2]*B[0];
    cross[2] = A[0]*B[1]-A[1]*B[0];

    return cross;//
}

template<class T> T mmpbsa_utils::dot_product(const T* A, const T* B, const size_t& ndim)
{
    T returnMe;
    for(size_t i = 0;i<ndim;i++)
        returnMe += A[0]*B[0];
}

//explicit instantiations
template std::valarray<size_t> mmpbsa_utils::zip<size_t>(std::valarray<size_t> const&, std::slice const&, std::valarray<size_t> const&, std::slice const&);
template bool mmpbsa_utils::contains<size_t>(std::vector<size_t, std::allocator<size_t> > const&, size_t const&);
template mmpbsa_t* mmpbsa_utils::cross_product<mmpbsa_t>(mmpbsa_t const*, mmpbsa_t const*, size_t const&);
template std::valarray<size_t> mmpbsa_utils::cumBoolSum<bool>(std::valarray<bool> const&);
template std::valarray<size_t> mmpbsa_utils::cumsum<size_t>(std::valarray<size_t> const&);
template std::vector<size_t, std::allocator<size_t> > mmpbsa_utils::take<size_t>(std::valarray<size_t> const&, std::valarray<size_t> const&, std::valarray<bool> const&);
template size_t mmpbsa_utils::find_first<int>(std::valarray<int> const&, int const&);
template std::valarray<size_t> mmpbsa_utils::cshift<size_t>(std::vector<size_t, std::allocator<size_t> > const&, int const&);
template std::valarray<size_t> mmpbsa_utils::zip<size_t>(std::vector<size_t, std::allocator<size_t> > const&, std::valarray<size_t> const&);

