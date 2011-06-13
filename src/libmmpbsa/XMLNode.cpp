
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


mmpbsa_utils::XMLNode* mmpbsa_utils::XMLNode::removeChild(XMLNode* to_be_removed) throw (mmpbsa::XMLParserException)
{
	XMLNode* left_of_delete;
	if(this->children == 0)
		return 0;

	for(left_of_delete = this->children;left_of_delete != 0;left_of_delete = left_of_delete->siblings)
	{
		if(left_of_delete->siblings == 0)//this means we've reached the end and to_be_removed is not a child of "this"
			throw mmpbsa::XMLParserException("mmpbsa_utils::XMLNode::removeChild: attempted to remove a non-child node.",mmpbsa::INVALID_XML_REQUEST);

		if(left_of_delete->siblings == to_be_removed)
			break;
	}
	if(left_of_delete->siblings == 0)
		return left_of_delete;

	left_of_delete->siblings = to_be_removed->siblings;
	to_be_removed->siblings = 0;
	delete to_be_removed;
	return left_of_delete;
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

#if 0 //old, less flexible method of sorting
std::vector<mmpbsa_utils::XMLNode*> quick_sort(mmpbsa_utils::XMLNode* node_to_sort,
		int sort_algo(mmpbsa_utils::XMLNode* partition_element, mmpbsa_utils::XMLNode* test_element))
{
	std::vector<mmpbsa_utils::XMLNode*> returnMe;
	std::vector<int> stack;

	if(node_to_sort == 0)
		throw mmpbsa::XMLParserException("quick_sort: given a null pointer.");

	mmpbsa_utils::XMLNode *it,*swapper;
	int curr_index = 0,index_length;
	int ins_inner_loop,ins_sort_index,median;
	int GREATER = 0x2,LESS = 0x4, EQUAL = 0x1;//retvals of sort_algo. Used as a greater than, less than, equal to bit map.

	//index array starting point.
	for(it = node_to_sort->children;it != 0;it = it->siblings)
		returnMe.push_back(it);

	index_length = returnMe.size() - 1;
	//quick sort algo
	for(;;)
	{
		if(index_length - 1 < 7)//suggested by Numerical Recipes to be the cutoff between quicksort and insert sort
		{
			for(ins_sort_index = curr_index + 1;ins_sort_index <= index_length;ins_sort_index++)
			{
				it = returnMe[ins_sort_index];
				for(ins_inner_loop = ins_sort_index - 1;ins_inner_loop >= curr_index; ins_inner_loop--)
				{
					if(sort_algo(returnMe[ins_inner_loop],it) & (LESS+EQUAL))
						break;
					returnMe[ins_inner_loop+1] = returnMe[ins_inner_loop];
				}
				returnMe[ins_inner_loop+1] = it;
			}
			if(stack.size() == 0)
				break;
			ins_sort_index = stack.back();stack.pop_back();
			index_length = stack.back();stack.pop_back();
		}
		else
		{
			median = (curr_index + index_length) >> 1;//more efficient then division.
			swapper = returnMe[median];
			returnMe[median] = returnMe[curr_index+1];
			returnMe[curr_index+1] = swapper;
			if(sort_algo(returnMe[curr_index],returnMe[index_length]) & GREATER)
			{
				swapper = returnMe[curr_index];
				returnMe[curr_index] = returnMe[index_length];
				returnMe[index_length] = swapper;
			}
			if(sort_algo(returnMe[curr_index+1],returnMe[index_length]) & GREATER)
			{
				swapper = returnMe[curr_index+1];
				returnMe[curr_index+1] = returnMe[index_length];
				returnMe[index_length] = swapper;
			}
			if(sort_algo(returnMe[curr_index],returnMe[curr_index+1]) & GREATER)
			{
				swapper = returnMe[curr_index];
				returnMe[curr_index] = returnMe[curr_index+1];
				returnMe[curr_index + 1] = swapper;
			}
			ins_inner_loop = curr_index+1;
			ins_sort_index = index_length;
			it = returnMe[curr_index+1];
			for(;;)
			{
				do ins_inner_loop++;while(sort_algo(returnMe[ins_inner_loop],it) & LESS);
				do ins_sort_index--;while(sort_algo(returnMe[ins_sort_index],it) & GREATER);
				if(ins_inner_loop > ins_sort_index)
					break;
				swapper = returnMe[ins_inner_loop];
				returnMe[ins_inner_loop] = returnMe[ins_sort_index];
				returnMe[ins_sort_index] = swapper;
			}
			returnMe[curr_index+1] = returnMe[ins_sort_index];
			returnMe[ins_sort_index] = it;
			if(index_length-ins_inner_loop+1 >= ins_sort_index-1)
			{
				stack.push_back(index_length);
				stack.push_back(ins_inner_loop);
				index_length = ins_sort_index - 1;
			}
			else
			{
				stack.push_back(ins_sort_index - 1);
				stack.push_back(curr_index);
				curr_index = ins_inner_loop;
			}
		}
	}

	return returnMe;
}
#endif

