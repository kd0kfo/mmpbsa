#ifndef MMPBSA_VECTOR
#define MMPBSA_VECTOR

#include <iostream>
#include <vector>
#include <algorithm>

namespace mmpbsa{

//Similar to Mead's Coord, but improved
class Vector : public std::vector<mmpbsa_t> {
public:
  Vector ();
  Vector (const value_type& xi, const value_type& yi, const value_type& zi);
  Vector(const value_type* coords,const Vector::size_type& dim);
  ~Vector () {}
  Vector& operator+= (const Vector& a);
  Vector& operator-= (const Vector& a);
  Vector& operator*= (const value_type& a);
  Vector& operator/= (const value_type& a);
  Vector operator- () const;
  bool operator> (const Vector& a) const;
  bool operator< (const Vector& a) const;

  const value_type& x()const{return this->at(0);}
  const value_type& y()const{return this->at(1);}
  const value_type& z()const{return this->at(2);}

  value_type distance(const mmpbsa::Vector& other_point)const;
  value_type modulus()const;

  friend std::ostream& operator<<(std::ostream& ost,const mmpbsa::Vector& v);

};

}//namespace mmpbsa

mmpbsa::Vector::value_type dot (const mmpbsa::Vector& a, const mmpbsa::Vector& b);
mmpbsa::Vector cross (const mmpbsa::Vector& a, const mmpbsa::Vector& b);

mmpbsa::Vector::value_type operator* (const mmpbsa::Vector& a, const mmpbsa::Vector& b);
mmpbsa::Vector operator* (const mmpbsa::Vector::value_type& b, const mmpbsa::Vector& a);
mmpbsa::Vector operator* (const mmpbsa::Vector& a, const mmpbsa::Vector::value_type& b);
mmpbsa::Vector operator- (const mmpbsa::Vector& a, const mmpbsa::Vector& b);
mmpbsa::Vector operator+ (const mmpbsa::Vector& a, const mmpbsa::Vector& b);

bool operator != (const mmpbsa::Vector& a, const mmpbsa::Vector& b);
bool operator == (const mmpbsa::Vector& a, const mmpbsa::Vector& b);
mmpbsa::Vector operator/ (const mmpbsa::Vector& a, const mmpbsa::Vector::value_type& b);


#endif
