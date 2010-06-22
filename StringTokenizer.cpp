#include "StringTokenizer.h"

utils::StringTokenizer::StringTokenizer(const std::string& string)
	{
		StringTokenizer temp(string," ",false);
		*this = temp;
	}

utils::StringTokenizer::StringTokenizer(const std::string& newString, const char * newDelim)
	{
		StringTokenizer temp(newString,newDelim,false);
		*this = temp;
	}

utils::StringTokenizer::StringTokenizer(const std::string& newString, const char * newDelim, bool keepDelim)
	{
		
	  delim = newDelim;
	  string = newString;
	  this->keepDelim = keepDelim;
	  index =  0;
	  tokenStream << newString;
	}

utils::StringTokenizer::StringTokenizer(const StringTokenizer& rhs)
	{
	  delim = rhs.getDelim();
	  index = rhs.index;
	  keepDelim = rhs.keepDelim;
	  string = rhs.string;
	  tokenStream << string;
	}

utils::StringTokenizer utils::StringTokenizer::operator=(const StringTokenizer& rhs)
	{
	  if(this == &rhs)
	    return *this;
	  delim = rhs.getDelim();
	  index = rhs.index;
	  keepDelim = rhs.keepDelim;
	  string = rhs.string;
	  tokenStream << string;
	  return *this;
	}

std::string utils::StringTokenizer::nextToken(bool keepD)
	{
	  if(tokenStream.eof())
	    throw TokenizerException("No more tokens in string tokenizer");

	  std::string returnMe;

          if(delim[0] == ' '  && !keepD)//optimized for ' ' delimiter by taking advantage of stringstream's operator>>
          {
              tokenStream >> returnMe;
              return returnMe;
          }
          
	  if(this->string[index] == delim[0] && keepD)
	    {
	      index++;
	      return delim;
	    }
	  
	  getline(tokenStream,returnMe,delim[0]);
	  index += returnMe.size();
	  return returnMe;
	}
	
bool utils::StringTokenizer::hasMoreTokens()
{
  return !tokenStream.eof();
}

std::string utils::StringTokenizer::peek()
{
  std::string bean = nextToken(this->keepDelim);
  this->index -= 1;
  tokenStream.seekg(-1*bean.size(),std::ios::cur);
  return bean;
}




