#include "Vector.h"


mmpbsa::Vector::Vector()
{
	push_back(0);
	push_back(0);
	push_back(0);
}


mmpbsa::Vector::Vector(const mmpbsa::Vector::value_type& xi, const mmpbsa::Vector::value_type& yi, const mmpbsa::Vector::value_type& zi)
{
	push_back(xi);
	push_back(yi);
	push_back(zi);
}

mmpbsa::Vector::Vector(const mmpbsa::Vector::value_type* coords,const mmpbsa::Vector::size_type& dim)
{
	for(mmpbsa::Vector::size_type i = 0; i < dim; i++)
		push_back(coords[i]);
}


mmpbsa::Vector& mmpbsa::Vector::operator+= (const mmpbsa::Vector& a)
{
	mmpbsa::Vector::size_type i,dim = std::min(a.size(),this->size());
	for(i = 0;i<dim;i++)
		this->at(i) += a.at(i);
	return *this;
}

mmpbsa::Vector& mmpbsa::Vector::operator-= (const mmpbsa::Vector& a)
{
	mmpbsa::Vector::size_type i,dim = std::min(a.size(),this->size());
	for(i = 0;i<dim;i++)
		this->at(i) -= a.at(i);
	return *this;
}

mmpbsa::Vector& mmpbsa::Vector::operator*= (const mmpbsa::Vector::value_type& a)
{
	mmpbsa::Vector::size_type i,dim = this->size();
	for(i = 0;i<dim;i++)
		this->at(i) *= a;
	return *this;
}

mmpbsa::Vector& mmpbsa::Vector::operator/= (const mmpbsa::Vector::value_type& a)
{
	mmpbsa::Vector::size_type i,dim = this->size();
	for(i = 0;i<dim;i++)
		this->at(i) /= a;
	return *this;
}

mmpbsa::Vector mmpbsa::Vector::operator- () const
{
	mmpbsa::Vector a;
	mmpbsa::Vector::size_type i,dim = this->size();
	for(i = 0;i<dim;i++)
		a.at(i) = -(this->at(i));
	return a;
}

bool mmpbsa::Vector::operator> (const Vector& a) const
{
	mmpbsa::Vector::size_type i,dim = std::min(a.size(),this->size());
	for(i = 0;i<dim;i++)
		if(this->at(i) <= a.at(i))
			return false;

	return true;
}

bool mmpbsa::Vector::operator< (const Vector& a) const
{
	mmpbsa::Vector::size_type i,dim = std::min(a.size(),this->size());
	for(i = 0;i<dim;i++)
		if(this->at(i) >= a.at(i))
			return false;

	return true;
}

mmpbsa::Vector operator+ (const mmpbsa::Vector& a, const mmpbsa::Vector& b)
{
	mmpbsa::Vector sum;
	mmpbsa::Vector::size_type i,dim = std::min(a.size(),b.size());
	for(i = 0;i<dim;i++)
		sum.at(i) = a.at(i) + b.at(i);
	return sum;
}

mmpbsa::Vector operator- (const mmpbsa::Vector& a, const mmpbsa::Vector& b)
{
	mmpbsa::Vector diff;
	mmpbsa::Vector::size_type i,dim = std::min(a.size(),b.size());
	for(i = 0;i<dim;i++)
		diff.at(i) = a.at(i) - b.at(i);
	return diff;
}


mmpbsa::Vector operator* (const mmpbsa::Vector& a, const mmpbsa::Vector::value_type& b)
{
	mmpbsa::Vector product;
	mmpbsa::Vector::size_type i,dim = a.size();
	for(i = 0;i<dim;i++)
		product.at(i) = a.at(i) * b;
	return product;
}

mmpbsa::Vector operator* (const mmpbsa::Vector::value_type& b, const mmpbsa::Vector& a)
{
	mmpbsa::Vector prod = a; prod*=b; return prod;
}

mmpbsa::Vector::value_type operator* (const mmpbsa::Vector& a, const mmpbsa::Vector& b)
{
	mmpbsa::Vector::value_type inner_prod = 0;
	mmpbsa::Vector::size_type i,dim = std::min(a.size(),b.size());
	for(i = 0;i<dim;i++)
		inner_prod += a.at(i) * b.at(i);
	return inner_prod;
}

mmpbsa::Vector cross (const mmpbsa::Vector& a, const mmpbsa::Vector& b)
{
	if(a.size() != b.size())
		throw mmpbsa::MMPBSAException("cross: the cross product requires Vectors of equal length.");
	if(a.size() != 3)
		throw mmpbsa::MMPBSAException("cross: the cross product has not yet been implemented for dimensions other than three.");

	mmpbsa::Vector prod(a.y()*b.z() - a.z()*b.y(), b.x()*a.z() - a.x()*b.z(), a.x()*b.y() - b.x()*a.y());
	return prod;
}

mmpbsa::Vector::value_type dot (const mmpbsa::Vector& a, const mmpbsa::Vector& b)
{
	return a * b;
}

mmpbsa::Vector operator/ (const mmpbsa::Vector& a, const mmpbsa::Vector::value_type& b)
{
	mmpbsa::Vector prod = a;
	prod /= b;
	return prod;
}

mmpbsa::Vector::value_type mmpbsa::Vector::distance(const mmpbsa::Vector& other_point)const
{
	mmpbsa::Vector diff = *this - other_point;
	return sqrt(diff * diff);
}

mmpbsa::Vector::value_type mmpbsa::Vector::modulus()const
{
	return sqrt(*this * *this);
}

bool operator == (const mmpbsa::Vector& a, const mmpbsa::Vector& b)
{
	mmpbsa::Vector::size_type i,dim = std::min(a.size(),b.size());
	for(i = 0;i<dim;i++)
		if(a.at(i) != b.at(i))
			return false;
	return true;
}

bool operator != (const mmpbsa::Vector& a, const mmpbsa::Vector& b)
{
	mmpbsa::Vector::size_type i,dim = std::min(a.size(),b.size());
	for(i = 0;i<dim;i++)
		if(a.at(i) != b.at(i))
			return true;
	return false;
}

std::ostream& operator<<(std::ostream& ost,const mmpbsa::Vector& v)
{
    ost << "(";
    mmpbsa::Vector::const_iterator val = v.begin();
    for(;val != v.end();val++)
    {
    	ost << *val;
    	if((val + 1) != v.end())
    		ost << ", ";
    }
    ost << ")";
    return ost;
  }




