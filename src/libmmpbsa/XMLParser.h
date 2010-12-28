/* 
 * Custom XML Parser for use within the MMPBSA program.
 *
 * This class was made with its end use in mind. Therefore it is not made to be
 * a all purpose, general public use XML parser. It works with documents made for
 * MD/MMPBSA. It will not:
 *      --Pay attention to comment tags. IMPORTANT
 *      --Pay attention to character sets, i.e. UTF-8.
 *      --Allow multiple tags on one line (though tags can span multiple lines)
 *      --Use flags. For example <a href="test">...</a> would have the name 'a href="test"'
 *      --Output tags of the for <stuff/>
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

#include "mmpbsa_io.h"
#include "mmpbsa_utils.h"
#ifndef CR_CHAR
#define CR_CHAR 0xd
#endif

namespace mmpbsa_utils{

class XMLParser {
public:
    /**
     * Creates an empty XML class with "NULL" as its tag.
     * 
     */
    XMLParser();

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
     * Creates an XML document using the given XMLNode.
     * 
     * @param head
     */
    XMLParser(XMLNode* head){this->head = head;}

    /**
     * Reads the file and creates an xml document tree. If there is an unexpected
     * error in the document, an mmpbsa::XMLParserException is thrown.
     * 
     * @param xmlFile
     */
    void parse(const std::string& xmlFile) throw (mmpbsa::XMLParserException);

    /**
     * Parses a line and turns it into a pair, where the first component of the
     * pair is the tag name and the second component is the text value, if such
     * a thing exists. If the line contains both the beginning and the end of the
     * tag isClosedTag is set to true. Otherwise it is false.
     * 
     * @param line
     * @param isClosedTag
     * @return
     */
    static std::pair<std::string,std::string> parseLine(const std::string& line, 
        bool& isClosedTag) throw (mmpbsa::XMLParserException);

    /**
     * Writes the XML Document to a file using toString().
     * 
     * @param fileName
     */
    void write(const std::string& fileName)const throw (mmpbsa::XMLParserException){write(fileName.c_str());}
    void write(const char* fileName)const throw (mmpbsa::XMLParserException);

    /**
     * Produces a string representation of the XML tree using XMLNode::toString()
     * for each element.
     * 
     * @return String represntation of the object.
     */
    std::string toString()const;

    /**
     * Takes the children of the specified node and makes a map of the names and content.
     * If a tag has no content, the value of the map is set to a null string, "".
     * 
     * @return
     */
    static std::map<std::string,std::string> mapNode(XMLNode const * const theNode) throw (mmpbsa::XMLParserException);

    /**
     * Returns a pointer to the document head.
     * 
     * @return Document Head
     */
    XMLNode * getHead(){return head;}
    const XMLNode * const getHead() const {return head;}

    /**
     * Detaches (and returns) the head node from the parser object, after which
     * the parser object's head is set to null. This allows one to remove the
     * node from the parser, preventing it from being deleted after the
     * parser's destructor is called.
     *
     * @return mmpbsa_utils::XMLNode* detached head node of document.
     */
    XMLNode* detachHead();

    /**
     * Replaces the head node of the XMLParser. Previous node is deleted.
     * 
     * @param newHead
     */
    void setHead(mmpbsa_utils::XMLNode* newHead){delete head;head = newHead;}

    /**
     * Returns the name of the head.
     * 
     * @return
     */
    std::string mainTag();

    /**
     * Calls delete head
     * 
     */
    virtual ~XMLParser();

    
private:
    XMLNode* head;///<top node in the tree 
};

};//end namespace mmpbsa_utils

#endif	/* XMLPARSER_H */

