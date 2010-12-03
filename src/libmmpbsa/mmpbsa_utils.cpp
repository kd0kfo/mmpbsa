#include "mmpbsa_utils.h"

int mmpbsa_utils::loadListArg(const std::string& values,std::vector<size_t>& array, const size_t& offset)
{
    using mmpbsa_utils::StringTokenizer;
    StringTokenizer valTokens(values,",");
    int currValue = 0;
    while(valTokens.hasMoreTokens())
    {
        std::istringstream curr(valTokens.nextToken());
        curr >> currValue;
        array.push_back(size_t(currValue) - offset);
    }
    return 0;
}

std::string mmpbsa_utils::toUpperCase(const std::string& bean)
{
    std::string returnMe = bean;
    for(size_t i = 0;i<bean.size();i++)
        returnMe[i] = std::toupper(bean[i]);
    return returnMe;
}

std::string mmpbsa_utils::toLowerCase(const std::string& bean)
{
    std::string returnMe = bean;
    for(size_t i = 0;i<bean.size();i++)
        returnMe[i] = std::tolower(bean[i]);
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

float * mmpbsa_utils::interaction_minmax(const std::valarray<mmpbsa_t>& acrds,
        const std::valarray<mmpbsa_t>& bcrds, const mmpbsa_t& cutoff)
{
    using std::valarray;
    using std::slice;
    using std::max;
    using std::min;
    
    if(acrds.size() % 3 != 0 || bcrds.size() % 3 != 0)
        throw mmpbsa::MMPBSAException("interaction_minmax was supply coordinates "
                "that were not 3-D.",mmpbsa::DATA_FORMAT_ERROR);

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
    float max_corner[3],min_corner[3];
    bool seenFirstInteractor = false;
    for(size_t i = 0;i<aflags.size();i++)
    {
        if(aflags[i])
        {
            if(!seenFirstInteractor)
            {
                for(size_t j = 0;j<3;j++)
                    min_corner[j] = max_corner[j] = acrds[3*i+j];
                seenFirstInteractor = true;
            }
            else
            {
                for(size_t j = 0;j<3;j++)
                {
                    max_corner[j] = max(acrds[3*i+j],mmpbsa_t(max_corner[j]));
                    min_corner[j] = min(acrds[3*i+j],mmpbsa_t(min_corner[j]));
                }
            }
        }
    }

    for(size_t i = 0;i<bflags.size();i+=3)
    {
        if(bflags[i])
        {
            for(size_t j = 0;j<3;j++)
            {
                max_corner[j] = max(bcrds[3*i+j],mmpbsa_t(max_corner[j]));
                min_corner[j] = min(bcrds[3*i+j],mmpbsa_t(min_corner[j]));
            }
        }
    }

    float * returnMe = new float[6];//x,y,z for min and max respectively
    for(size_t i = 0;i<3;i++)
    {
        returnMe[i] = min_corner[i];
        returnMe[i+3] = max_corner[i];
    }
    return returnMe;
}
float mmpbsa_utils::lookup_radius(const std::string& atomName,
       const std::map<std::string,float>& radiusMap)
            throw (mmpbsa::MMPBSAException)
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
    {
    	return radiusMap.find(atomName)->second;
    }

    //Not found by atomName. Check atom name only entries
    std::vector<mmpbsa_t> possibleMatches;
    string theAtom = trimString(atomName);
    for(map<std::string,float>::const_iterator it = radiusMap.begin();it != radiusMap.end();it++)
    {
        if(it->first == theAtom)//direct match
        {
            possibleMatches.push_back(it->second);
        }
    }
    //see what we found.
    if(possibleMatches.size() == 0)//see above note on shared radii
    {
        float deeperSearch = lookup_radius(theAtom.erase(theAtom.size()-1),radiusMap);
        if(deeperSearch == -1)
        {
            std::string error = "No radius found for '%s' in Radii Map" + atomName;
            throw mmpbsa::MMPBSAException(error,mmpbsa::DATA_FORMAT_ERROR);
        }
        return deeperSearch;
    }
    if(possibleMatches.size() > 1)
    {
        std::ostringstream error;
        error << possibleMatches.size() << " matches for '" << atomName << "' in Radii Map";
        throw mmpbsa::MMPBSAException(error,mmpbsa::DATA_FORMAT_ERROR);
    }
    
    //iff there is one match.
    return float(possibleMatches[0]);
}

std::string mmpbsa_utils::get_human_time()
{
	time_t rawtime;
	struct tm * timeinfo;
	time ( &rawtime );
	timeinfo = localtime ( &rawtime );
	return asctime (timeinfo);
}

