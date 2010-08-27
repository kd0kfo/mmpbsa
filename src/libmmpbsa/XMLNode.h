/* 
 * Node used by XMLParser. Stores data and pointers to other 
 * nodes in the document tree, which XMLParser represents.
 *
 * Created by David Coss (David.Coss@stjude.org) 2010
 */

#ifndef XMLNODE_H
#define	XMLNODE_H

#include <string>

namespace mmpbsa_utils{

class XMLNode {
public:
    /**
     * Creates an empty node.
     */
    XMLNode();

    /**
     * Creates a childless node, with the specified name and text content.
     * 
     * @param name
     * @param text
     */
    XMLNode(const std::string& name,const std::string& text = "");

    /**
     * Deletes itself and all of its siblings and children. Therefore, only the
     * top most node should be deleted.
     */
    virtual ~XMLNode();

    /**
     * Gives the content of the node.
     * 
     * @return Constant reference to the content of the XMLNode.
     */
    const std::string& getText()const{return text;}
    
    /**
     * Sets the content of the tag.
     * @param newText
     */
    void setText(const std::string& newText){text = newText;}

    /**
     * Gives the name of the Node.
     *
     * @return Constant reference to the name string.
     */
    const std::string& getName()const{return name;}

    /**
     * Sets the name of the node.
     * 
     * @param newChildren
     * @return
     */
    void setName(const std::string& newName){name = newName;}

    /**
     * Returns the number of children the that particular node.
     * For a full depth count of all branches off of this node, see
     * hasRelatives()
     *
     * @return size_t number of children of the node.
     */
    size_t hasChildren()const;

    /**
     * Returns a count of all decendents
     * @return 
     */
    size_t hasRelatives()const;

    /**
     * Returns the number of siblings of the node.
     * 
     * @returns size_t number of siblings of the node.
     */
    size_t hasSiblings()const;

    /**
     * Attaches children to the node. Children that were previously attached are
     * returned.
     *
     * @param newChildren
     * @return pointer to previous children.
     */
    XMLNode* attachChildren(XMLNode* newChildren);

    /**
     * Attaches siblings to the node. Previous siblings are returned.
     *
     * @param newSiblings
     * @return pointer to previous siblings.
     */
    XMLNode* attachSiblings(XMLNode* newSiblings);

    /**
     * Replaces children of the node. Previous children are deleted.
     *
     * @param newChildren
     */
    void replaceChildren(XMLNode* newChildren){delete (attachChildren(newChildren));}

    /**
     * Replaces siblings of the node. Previous siblings are deleted.
     *
     * @param newSiblings
     */
    void replaceSiblings(XMLNode* newSiblings){delete (attachSiblings(newSiblings));}

    /**
     * Inserts a child at the end of the children.
     *
     * @param newChild
     */
    void insertChild(XMLNode* newChild);

    /**
     * Inserts a new node with the given name and text content at then end of the
     * children
     *
     * @param name
     * @param text
     */
    void insertChild(const std::string& name, const std::string& text){insertChild(new XMLNode(name,text));}
    
    XMLNode* parent;///<Parent of the node.
    XMLNode* children;///<Pointer to the first child (not an array, cf destructor)
    XMLNode* siblings;///<Pointer to first sibling (not an array, cf destructor)
    XMLNode* end(){return 0;}///<To reproduce use of iterator.

    /**
     * Produces an XML document string based using the name and text content
     * values of the node and all of its siblings and children.
     *
     * An optional offset may be used to deliniate the branches of the document
     * tree. Default is three blank spaces.
     * 
     * @param offset
     * @return string based on the XML document format.
     */
    std::string toString(const std::string& offset = "   ")const;

    /**
     * Produces a pair containing the name and text of the XMLNode.
     *
     * @return std::pair<std::string,std::string>(name,text)
     */
    std::pair<std::string, std::string> toPair()const
    {
        return std::pair<std::string,std::string>(name,text);
    }
private:
    std::string text;
    std::string name;
};

}//end namespace mmpbsa_utils

/**
 * Performs the given function on all siblings from beginning until end (exclusive).
 * Since the end point is exclusive, end can equal null, which will essentially
 * iteratorate through all siblings.
 * 
 * @param beginning
 * @param end
 * @param function
 */
void foreach(mmpbsa_utils::XMLNode* beginning, mmpbsa_utils::XMLNode* end,
        void function(mmpbsa_utils::XMLNode*));

#endif	/* XMLNODE_H */

