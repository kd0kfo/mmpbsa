/* 
 * Node used by XMLParser. Stores data and pointers to other 
 * nodes in the document tree, which XMLParser represents.
 *
 * Created by David Coss <David.Coss@stjude.org> 2010
 */

#ifndef XMLNODE_H
#define	XMLNODE_H

#include <string>

namespace mmpbsa_utils{

class XMLNode {
public:
    XMLNode();
    XMLNode(const std::string& name,const std::string& text = "");
    virtual ~XMLNode();

    //Returns old children
    const std::string& getText()const{return text;}
    void setText(const std::string& newText){text = newText;}
    const std::string& getName()const{return name;}
    void setName(const std::string& newName){name = newName;}
    size_t hasChildren()const;
    size_t hasSiblings()const;
    XMLNode* attachChildren(XMLNode* newChildren);
    XMLNode* attachSiblings(XMLNode* newSiblings);
    void replaceChildren(XMLNode* newChildren){delete (attachChildren(newChildren));}
    void replaceSiblings(XMLNode* newSiblings){delete (attachSiblings(newSiblings));}
    void insertChild(XMLNode* newChild);
    void insertChild(const std::string& name, const std::string& text){insertChild(new XMLNode(name,text));}
    XMLNode* parent;
    XMLNode* children;//Pointer to first child (not an array, cf destructor)
    XMLNode* siblings;//pointer to first sibling (not an array, cf destructor)
    std::string toString(const std::string& offset = "")const;
private:
    std::string text;
    std::string name;
};

}//end namespace mmpbsa_utils

#endif	/* XMLNODE_H */

