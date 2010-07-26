/* 
 * Wrapper for libxml2, customized for use within the MMPBSA program.
 *
 * Created by David Coss (David.Coss@stjude.org) 2010
 * 
 */

#ifndef XMLPARSER_H
#define	XMLPARSER_H

#include <fstream>
#include <string>
#include <map>
#include "libxml/parser.h"//include libxml2 in Includes path.
#include "libxml/tree.h"
#include "libxml/xmlstring.h"

#include "mmpbsa_exceptions.h"
#include "mmpbsa_io.h"

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

    /**
     * Reads the file and creates an xml document tree.
     * 
     * @param xmlFile
     */
    void parse(const std::string& xmlFile);

    /**
     * Writes the XML Document.
     * 
     * @param fileName
     */
    void write(const std::string& fileName){write(fileName.c_str());}
    void write(const char* fileName);

    /**
     * Takes the children of the root node and make a map of the names and content.
     * If a tag has no content, the value of the map is set to a null string, "".
     * 
     * @return
     */
    std::map<std::string,std::string> getChildren();


    virtual ~XMLParser();

    
private:
    xmlDocPtr doc; /* the document tree */
    xmlNodePtr head;/* top node in the tree */
};

class XMLParserException : public MMPBSAException
{
    public:
    /**
     * Exception for when something goes wrong with reading the XML data
     *
     * @param error
     */
    XMLParserException(const std::string& error) : MMPBSAException( error){}

    XMLParserException(const std::string& error, const MMPBSAErrorTypes& errorType)
        : MMPBSAException(error,errorType){}

    const char* identifier(){return "XML Parser Error";}
};

#endif	/* XMLPARSER_H */

