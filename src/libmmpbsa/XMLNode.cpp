
#include "XMLNode.h"
#include "XMLParser.h"

mmpbsa_utils::XMLNode::XMLNode() {
    parent = children = siblings = 0;
    name = "NULL";
    text = "";
}

mmpbsa_utils::XMLNode::XMLNode(const std::string& name,const std::string& text)
{
    parent = children = siblings = 0;
    this->name = name;
    this->text = text;
}

mmpbsa_utils::XMLNode::~XMLNode() {
    delete siblings;
    delete children;
}

size_t mmpbsa_utils::XMLNode::hasChildren()const
{
    if(children == 0)
        return 0;
    return children->hasSiblings() + 1;//one for children itself.
}

size_t mmpbsa_utils::XMLNode::hasRelatives()const
{
    if(children == 0)
        return 0;
    mmpbsa_utils::XMLNode * currNode = children;
    size_t returnMe = children->hasSiblings()+1;
    while(currNode)
    {
        returnMe += currNode->hasChildren();
        currNode = currNode->siblings;
    }
    return returnMe;
}

size_t mmpbsa_utils::XMLNode::hasSiblings()const
{
    if(siblings == 0)
        return 0;
    return 1+siblings->hasSiblings();
}

mmpbsa_utils::XMLNode* mmpbsa_utils::XMLNode::attachChildren(mmpbsa_utils::XMLNode* newChildren)
{
    XMLNode* returnMe = children;
    children = newChildren;
    XMLNode* currNode = children;
    for(;currNode != 0;currNode = currNode->siblings)
        currNode->parent = this;
    for(currNode = returnMe;currNode != 0;currNode = currNode->siblings)
        currNode->parent = 0;
    
    return returnMe;
}

mmpbsa_utils::XMLNode* mmpbsa_utils::XMLNode::attachSiblings(mmpbsa_utils::XMLNode* newSiblings)
{
    mmpbsa_utils::XMLNode* returnMe = siblings;
    siblings = newSiblings;
    XMLNode* currNode = siblings;
    for(;currNode != 0;currNode = currNode->siblings)
        currNode->parent = parent;
    for(currNode = returnMe;currNode != 0;currNode = currNode->siblings)
        currNode->parent = 0;
    return returnMe;
}

void mmpbsa_utils::XMLNode::insertChild(mmpbsa_utils::XMLNode* newChild)
{
    if(children == 0)
        children = newChild;
    else
    {
        XMLNode* currNode = children;
        while(currNode->siblings != 0)
            currNode = currNode->siblings;
        currNode->siblings = newChild;
    }
    newChild->parent = this;
}


std::string mmpbsa_utils::XMLNode::toString(const std::string& offset) const
{
    std::string returnMe = offset+"<" + name + ">";
    if(text.size() != 0)
    {
        if(children != 0)
            returnMe += "\n"+offset+offset;
        returnMe += text;
    }

    if(children != 0)
    {
        returnMe += "\n" + children->toString( (parent == 0) ? offset : offset + offset )
                + "\n"+offset;
    }
    
    returnMe += "</" + name + ">";
    if(siblings != 0)
        returnMe += "\n" + siblings->toString(offset);
    
    return returnMe;
}

void foreach(mmpbsa_utils::XMLNode* beginning, mmpbsa_utils::XMLNode* end,
        void function(mmpbsa_utils::XMLNode*))
{
    for(mmpbsa_utils::XMLNode* currNode = beginning;currNode != end;currNode = currNode->siblings)
    {
        if(currNode == 0)
            throw mmpbsa::XMLParserException("foreach was given a null pointer.",mmpbsa::DATA_FORMAT_ERROR);
        
        function(currNode);
    }
}