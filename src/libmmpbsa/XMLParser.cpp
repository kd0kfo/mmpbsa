#include "XMLParser.h"

mmpbsa_utils::XMLParser::XMLParser() {
    doc = 0;
    head = 0;
}

mmpbsa_utils::XMLParser::XMLParser(const mmpbsa_utils::XMLParser& orig) {
    doc = orig.doc;
    head = orig.head;
}

mmpbsa_utils::XMLParser::XMLParser(const std::string& rootName, const std::map<std::string,std::string>& docMap)
{
    using std::string;
    doc = xmlNewDoc(BAD_CAST "1.0");
    head = xmlNewNode(NULL, BAD_CAST rootName.c_str());
    xmlDocSetRootElement(doc, head);
    
    for(std::map<string,string>::const_iterator it = docMap.begin();it != docMap.end();it++)
        xmlNewChild(head,NULL,BAD_CAST it->first.c_str(),BAD_CAST it->second.c_str());
}

mmpbsa_utils::XMLParser::XMLParser(const std::map<std::string,std::string>& docMap)
{
    using std::string;
    if(!doc)
    {
        doc = xmlNewDoc(BAD_CAST "1.0");
        if(!head)
            delete head;
        head = xmlNewNode(NULL, BAD_CAST "xml_doc");
        xmlDocSetRootElement(doc, head);
    }

    for(std::map<string,string>::const_iterator it = docMap.begin();it != docMap.end();it++)
        xmlNewChild(head,NULL,BAD_CAST it->first.c_str(),BAD_CAST it->second.c_str());
}


mmpbsa_utils::XMLParser::~XMLParser() {
    xmlFreeDoc(doc);
    xmlCleanupParser();
}


void mmpbsa_utils::XMLParser::parse(const std::string& xmlFile)
{
    this->doc = xmlReadFile(xmlFile.c_str(), NULL, 0);
    if (doc == NULL) {
        char error[256];
        sprintf(error, "Failed to parse %s\n", xmlFile.c_str());
        throw mmpbsa::XMLParserException(error,mmpbsa::FILE_READ_ERROR);
    }
    this->head = xmlDocGetRootElement(doc);
}

void mmpbsa_utils::XMLParser::write(const char* fileName)
{
    xmlSaveFormatFileEnc(fileName, doc, "UTF-8", 1);
}

std::map<std::string,std::string> mmpbsa_utils::XMLParser::getChildren() const
{
    std::map<std::string,std::string> returnMe;
    std::pair<std::string,std::string> currPair;
    for(xmlNodePtr currNode = head->children;currNode;currNode = currNode->next)
    {

        if(currNode->type == ::XML_ELEMENT_NODE)
        {
            currPair.first = (char*) currNode->name;
            if(currNode->children && currNode->children->type == ::XML_TEXT_NODE)
                currPair.second = (char*) currNode->children->content;
            returnMe.insert(currPair);
                
        }
        else if(currNode->type == ::XML_TEXT_NODE)
            returnMe[(char*) currNode->content] = "";
    }
    return returnMe;
}


