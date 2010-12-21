#include "XMLParser.h"

mmpbsa_utils::XMLParser::XMLParser() {
    head = 0;
}

mmpbsa_utils::XMLParser::XMLParser(const std::string& rootName, const std::map<std::string,std::string>& docMap)
{
    using std::string;
    head = new XMLNode(rootName,"");
    
    XMLNode* currNode = 0;
    for(std::map<string,string>::const_iterator it = docMap.begin();it != docMap.end();it++)
    {
        if(currNode == 0)
        {
            currNode = new XMLNode(it->first,it->second);
            head->children = currNode;
        }
        else
        {
            currNode->replaceSiblings(new XMLNode(it->first,it->second));
            currNode = currNode->siblings;
        }
    }
}

mmpbsa_utils::XMLParser::XMLParser(const std::map<std::string,std::string>& docMap)
{
    using std::string;
    head = new XMLNode("NULL","");

    XMLNode* currNode = 0;
    for(std::map<string,string>::const_iterator it = docMap.begin();it != docMap.end();it++)
    {
        if(currNode == 0)
        {
            currNode = new XMLNode(it->first,it->second);
            head->children = currNode;
        }
        else
        {
            currNode->replaceSiblings(new XMLNode(it->first,it->second));
            currNode = currNode->siblings;
        }
    }
}


mmpbsa_utils::XMLParser::~XMLParser() {
    delete head;
}

void mmpbsa_utils::XMLParser::parse(const std::string& xmlFilename) throw (mmpbsa::XMLParserException)
{
    delete head;head = 0;
    std::fstream xmlFile(xmlFilename.c_str(),std::ios::in);
    if(!xmlFile.good())
    {
        std::string error = "Could not open xml file: " + xmlFilename;
        throw mmpbsa::XMLParserException(error,mmpbsa::FILE_IO_ERROR);
    }

    using mmpbsa_io::getNextLine;
    using mmpbsa_utils::trimString;
    using std::string;
    using std::pair;
    using mmpbsa_utils::XMLNode;

    bool isClosedTag = false;
    pair<string,string> currLine = parseLine(getNextLine(xmlFile),isClosedTag);

    head = new XMLNode(currLine.first,currLine.second);
    XMLNode* currNode = head;
    if(isClosedTag)
        currNode = 0;
    while(!xmlFile.eof())
    {
        currLine = parseLine(getNextLine(xmlFile),isClosedTag);
        if(currLine.first.size() > 0)
        {
            if(currLine.first == currNode->getName())
            {
                currNode = currNode->parent;//in this case, we're closing what was until now an unclosed tag.
                continue;
            }
            XMLNode* newNode = new XMLNode(currLine.first,currLine.second);
            if(currNode == 0)
            {
                std::string error = "XML Error in " + xmlFilename +
                        ": Missing top level tag. " + head->getName() + "should "
                        "not have a sibling without a parent.";
                throw mmpbsa::XMLParserException(error,mmpbsa::BAD_XML_TAG);
            }
            currNode->insertChild(newNode);
            if(!isClosedTag)
                currNode = newNode;
        }
        else
        {
            if(currLine.second.size() == 0)
                continue;
            currNode->setText(currNode->getText() + currLine.second);
        }
    }
    if(currNode != 0)
    {
        std::string error = xmlFilename + " is a malformed XML file. Missing an end tag "
                "to go with " + head->getName();
        throw mmpbsa::XMLParserException(error,mmpbsa::BAD_XML_TAG);
    }
}

std::pair<std::string,std::string> mmpbsa_utils::XMLParser::parseLine(const std::string& _line,
        bool& isClosedTag) throw (mmpbsa::XMLParserException)
{
    std::pair<std::string,std::string> returnMe("","");
    using std::string;
    string line = mmpbsa_utils::trimString(_line);
    while(line.find(CR_CHAR) != string::npos)//CR char
    	line.erase(line.find(CR_CHAR),1);

    if(line.size() == 0)
        return returnMe;
    if(line[0] != '<' || line.find_first_of('>') == string::npos)
    {
        returnMe.second = line;
        isClosedTag = false;
        return returnMe;
    }

    size_t tagEnd = line.find_first_of('>');
    returnMe.first = line.substr(1,tagEnd-1);
    line.erase(0,tagEnd+1);
    if(line[0] == '?')
        return std::pair<string,string>("","");

    if(returnMe.first[returnMe.first.size()-1] == '/' || returnMe.first[0] == '/')
    {
        isClosedTag = true;
        returnMe.first.erase(0,1);
        return returnMe;
    }

    size_t endTagPos = line.find_first_of('<');
    if(endTagPos == string::npos)
    {
        isClosedTag = false;
        return returnMe;
    }

    returnMe.second = line.substr(0,endTagPos);
    endTagPos = line.find_first_of('/',endTagPos);
    line = line.substr(endTagPos+1);
    line.erase(line.end()-1);
    if(line.size() > 0 && line != returnMe.first)
    {
        string error = "XML Parser does not yet support multiple tags on one line.\n"
                + _line;
        throw mmpbsa::XMLParserException(error,mmpbsa::BAD_XML_TAG);
    }

    isClosedTag = true;
    return returnMe;
}

void mmpbsa_utils::XMLParser::write(const char* fileName)const throw (mmpbsa::XMLParserException)
{
    std::fstream outfile(fileName,std::ios::out);
    if(!outfile.good())
    {
        char error[256];
        sprintf(error,"Could not open: %s",fileName);
        throw mmpbsa::XMLParserException(error,mmpbsa::FILE_IO_ERROR);
    }
    outfile << toString();
    outfile.close();
}

std::string mmpbsa_utils::XMLParser::toString()const
{
    if(head == 0)
        return "";
    return head->toString();
}

std::map<std::string,std::string> mmpbsa_utils::XMLParser::mapNode(XMLNode const * const theNode) throw (mmpbsa::XMLParserException)
{
    if(theNode == 0)
        throw mmpbsa::XMLParserException("Null XMLNode pointer provided to "
                "mapNode",mmpbsa::DATA_FORMAT_ERROR);
    
    std::map<std::string,std::string> returnMe;
    if(theNode->children == 0)
        return returnMe;
    
    for(XMLNode* currNode = theNode->children;currNode != 0;currNode = currNode->siblings)
        returnMe[currNode->getName()] = currNode->getText();

    return returnMe;
}

std::string mmpbsa_utils::XMLParser::mainTag()
{
    if(this->head)
        return this->head->getName();

    return "";
}

mmpbsa_utils::XMLNode* mmpbsa_utils::XMLParser::detachHead()
{
	mmpbsa_utils::XMLNode* detached = head;
	head = 0;
	return detached;
}
