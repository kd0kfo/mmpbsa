#include "mmpbsa_utils.h"

std::string mmpbsa_utils::toUpperCase(const std::string& bean)
{
    std::string returnMe = bean;
    for(size_t i = 0;i<bean.size();i++)
        returnMe[i] = std::toupper(bean[i]);
    return returnMe;
}

std::string mmpbsa_utils::trimString(const std::string& theString)
{
    std::string bean = theString;
    size_t lastpos = bean.find_last_not_of(" \t");
    size_t firstpos = bean.find_first_not_of(" \t");
    if(firstpos == bean.npos || lastpos == bean.npos)
    {
        bean = "";
    }
    else
        bean = bean.substr(firstpos,lastpos - firstpos + 1);

    return bean;
}

Coord * mmpbsa_utils::interaction_minmax(const std::valarray<mmpbsa_t>& acrds,
        const std::valarray<mmpbsa_t>& bcrds, const mmpbsa_t& cutoff)
{
    using std::valarray;
    using std::slice;
    using std::max;
    using std::min;
    
    if(acrds.size() % 3 != 0 || bcrds.size() % 3 != 0)
        throw MMPBSAException("interaction_minmax was supply coordinates "
                "that were not 3-D.",DATA_FORMAT_ERROR);

    mmpbsa_t cutsqrd = cutoff*cutoff;
    valarray<bool> aflags(false,size_t(acrds.size()/3));
    valarray<bool> bflags(false,size_t(bcrds.size()/3));

    mmpbsa_t rsqrd = 0;//valarray<mmpbsa_t> rsqrd(0.0,numCoords);
    for(size_t i = 0;i<bflags.size();i++)
    {
        for(size_t j = 0;j<aflags.size();j++)
        {
            rsqrd = pow(bcrds[3*i]-acrds[3*j],2) + pow(bcrds[3*i+1]-acrds[3*j+1],2)
                    + pow(bcrds[3*i+2]-acrds[3*j+2],2);
            if(rsqrd < cutsqrd)
            {
                bflags[i] = true;
                aflags[j] = true;
            }
        }

    }
    Coord max_corner,min_corner;
    bool seenFirstInteractor = false;
    for(size_t i = 0;i<aflags.size();i++)
    {
        if(aflags[i])
        {
            if(!seenFirstInteractor)
            {
                max_corner.x = acrds[3*i];
                max_corner.y = acrds[3*i+1];
                max_corner.z = acrds[3*i+2];
                min_corner = max_corner;
                seenFirstInteractor = true;
            }
            else
            {
                max_corner.x = max(acrds[3*i],mmpbsa_t(max_corner.x));
                max_corner.y = max(acrds[3*i+1],mmpbsa_t(max_corner.y));
                max_corner.z = max(acrds[3*i+2],mmpbsa_t(max_corner.z));
                min_corner.x = min(acrds[3*i],mmpbsa_t(min_corner.x));
                min_corner.y = min(acrds[3*i+1],mmpbsa_t(min_corner.y));
                min_corner.z = min(acrds[3*i+2],mmpbsa_t(min_corner.z));
            }
        }
    }

    for(size_t i = 0;i<bflags.size();i+=3)
    {
        if(bflags[i])
        {
            max_corner.x = max(bcrds[3*i],mmpbsa_t(max_corner.x));
            max_corner.y = max(bcrds[3*i+1],mmpbsa_t(max_corner.y));
            max_corner.z = max(bcrds[3*i+2],mmpbsa_t(max_corner.z));
            min_corner.x = min(bcrds[3*i],mmpbsa_t(min_corner.x));
            min_corner.y = min(bcrds[3*i+1],mmpbsa_t(min_corner.y));
            min_corner.z = min(bcrds[3*i+2],mmpbsa_t(min_corner.z));
        }
    }

    Coord * returnMe = new Coord[2];
    returnMe[0] = min_corner;
    returnMe[1] = max_corner;
    return returnMe;
}

mmpbsa_t mmpbsa_utils::lookup_radius(const std::string& atomName,
        const std::map<std::string,mmpbsa_t>& radiusMap)
            throw (MMPBSAException)
{
    using std::string;
    using std::map;
    using mmpbsa_utils::trimString;

    //If atomic radii are the same for element variations with common first letters, such
    //as Carbon-alpha and Carbon-beta, they are stored in the .siz by that
    //shared abbreviation, i.e. "C". Therefore if a radius is not found by the
    //whole name, remove one character and search again. A null string indicates
    //this failed as well.
    if(atomName.size() == 0)
        return -1;

    //A direct name match is preferred. Otherwise test for untrimmed keys and/or ambiguities.
    if(radiusMap.find(atomName) != radiusMap.end())
        return radiusMap.at(atomName);
    
    //Not found by atomName. Check atom name only entries
    std::vector<mmpbsa_t> possibleMatches;
    string theAtom = trimString(atomName);
    for(map<std::string,mmpbsa_t>::const_iterator it = radiusMap.begin();it != radiusMap.end();it++)
    {
        if(it->first == theAtom)//direct match
        {
            possibleMatches.push_back(it->second);
        }
    }

    //see what we found.
    if(possibleMatches.size() == 0)//see above note on shared radii
    {
        mmpbsa_t deeperSearch = lookup_radius(theAtom.erase(theAtom.size()-1),radiusMap);
        if(deeperSearch == -1)
        {
            char error[256];
            sprintf(error,"No radius found for '%s' in Radii Map",atomName.c_str());
            throw MMPBSAException(error,DATA_FORMAT_ERROR);
        }
        return deeperSearch;
    }
    if(possibleMatches.size() > 1)
    {
        char error[256];
        sprintf(error,"%d matches for '%s' in Radii Map",possibleMatches.size(), atomName.c_str());
        throw MMPBSAException(error,DATA_FORMAT_ERROR);
    }
    
    //iff there is one match.
    return possibleMatches[0];
}

