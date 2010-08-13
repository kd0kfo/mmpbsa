/* 
 * Custom XML Parser for use within the MMPBSA program.
 *
 * Created by David Coss (David.Coss@stjude.org) 2010
 * 
 */

#ifndef XMLPARSER_H
#define	XMLPARSER_H

#include <fstream>
#include <string>
#include <map>
#include "XMLNode.h"

#include "mmpbsa_exceptions.h"
#include "mmpbsa_io.h"


namespace mmpbsa{

class XMLParserException : public mmpbsa::MMPBSAException
{
    public:
    /**
     * Exception for when something goes wrong with reading the XML data
     *
     * @param error
     */
    XMLParserException(const std::string& error) : mmpbsa::MMPBSAException( error){}

    XMLParserException(const std::string& error, const mmpbsa::MMPBSAErrorTypes& errorType)
        : mmpbsa::MMPBSAException(error,errorType){}

    const char* identifier(){return "XML Parser Error";}
};

};//end namespace mmpbsa

namespace mmpbsa_utils{

class XMLParser {
public:
    XMLParser();

    XMLParser(const XMLParser& orig);

    /**
     * Creates an XML Document Tree with tags and content based on the values of
     * the provided map. The root tag is set to the value of rootName.
     * 
     * @param rootName
     * @param docMap
     */
    XMLParser(const std::string& rootName, const std::map<std::string,std::string>& docMap);

    /**
     * Creates an XML Document Tree with tags and content based on the values of
     * the provided map.
     * 
     * @param docMap
     */
    XMLParser(const std::map<std::string,std::string>& docMap);

    XMLParser(XMLNode* head){this->head = head;}

    /**
     * Reads the file and creates an xml document tree.
     * 
     * @param xmlFile
     */
    void parse(const std::string& xmlFile) throw (mmpbsa::XMLParserException);

    static std::pair<std::string,std::string> parseLine(const std::string& line, 
        bool& isClosedTag) throw (mmpbsa::XMLParserException);

    /**
     * Writes the XML Document.
     * 
     * @param fileName
     */
    void write(const std::string& fileName)const{write(fileName.c_str());}
    void write(const char* fileName)const;

    std::string toString()const;

    /**
     * Takes the children of the root node and make a map of the names and content.
     * If a tag has no content, the value of the map is set to a null string, "".
     * 
     * @return
     */
    static std::map<std::string,std::string> mapNode(XMLNode const * const theNode) throw (mmpbsa::XMLParserException);

    const XMLNode * getHead(){return head;}

    std::string mainTag();

    virtual ~XMLParser();

    
private:
    XMLNode* head;/* top node in the tree */
};

};//end namespace mmpbsa_utils

#endif	/* XMLPARSER_H */

