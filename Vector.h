/*contains declaration and implementation for the Vector class
This class is used for defining vectors in the mathematical sense. i.e. that have an x,y, and z component and can be used for 
vector operations.
*/

#ifndef _Vector_h
#define _Vector_h

#include <math.h>
#include <iostream>
#include <boost/mpi.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>


using namespace std;

class Vector;

struct VaddV {
  const Vector & v1;
  const Vector & v2;
  VaddV(const Vector & _v1,const Vector & _v2) : v1(_v1), v2(_v2) {}
};

struct VsubV {
  const Vector & v1;
  const Vector & v2;
  VsubV(const Vector & _v1,const Vector & _v2) : v1(_v1), v2(_v2) {}
};
struct DmulV {
  const double c;
  const Vector & v;
  DmulV(const double _c,const Vector & _v) : c(_c), v(_v) {}
};

struct Dmul_VaddV {
  const double c;
  const VaddV & v;
  Dmul_VaddV(const double _c,const VaddV & _v) : c(_c), v(_v) {}
};

struct Dmul_VsubV {
  const double c;
  const VsubV & v;
  Dmul_VsubV(const double _c,const VsubV & _v) : c(_c), v(_v) {}
};

struct DmulV_addV {
  const DmulV & dv;
  const Vector & v;
  DmulV_addV(const DmulV & _dv, const Vector & _v) : dv(_dv), v(_v) {}
};

struct DmulV_sub_DmulV {
  const DmulV & dv1;
  const DmulV & dv2;
  DmulV_sub_DmulV(const DmulV & _dv1, const DmulV & _dv2) : dv1(_dv1), dv2(_dv2) {}
};

struct DmulV_add_DmulV {
  const DmulV & dv1;
  const DmulV & dv2;
  DmulV_add_DmulV(const DmulV & _dv1, const DmulV & _dv2) : dv1(_dv1), dv2(_dv2) {}
};

struct VtimesV {
  const Vector & v1;
  const Vector & v2;
  VtimesV(const Vector & _v1,const Vector & _v2) : v1(_v1), v2(_v2) {}
};

class Vector {
  friend istream & operator >> (istream & is, Vector & v)
  {
    is >> v._x >> v._y >> v._z;
    return is;
  }
  friend ostream & operator << (ostream & os, const Vector & v)
  {
    os << v.x() << " " << v.y() << " " << v.z();
    return os;
  }
  friend VaddV operator + (const Vector & v1, const Vector & v2)
  {
    return VaddV(v1,v2);
  }
  friend DmulV_addV operator + (const DmulV & dv, const Vector & v)
  {
    return DmulV_addV(dv,v);
  }
  friend DmulV_addV operator + (const Vector & v, const DmulV & dv)
  {
    return DmulV_addV(dv,v);
  }
  friend DmulV_add_DmulV operator + (const DmulV & dv1, const DmulV & dv2)
  {
    return DmulV_add_DmulV(dv1,dv2);
  }
  friend VsubV operator - (const Vector & v1, const Vector & v2)
  {
    return VsubV(v1,v2);
  }
  friend DmulV_sub_DmulV operator - (const DmulV & dv1, const DmulV & dv2)
  {
    return DmulV_sub_DmulV(dv1,dv2);
  }
  friend DmulV operator * (double c, const Vector & v){
    return DmulV(c,v);
  }
  friend Dmul_VaddV operator * (double c, const VaddV & v){
    return Dmul_VaddV(c,v);
  }
  friend Dmul_VsubV operator * (double c, const VsubV & v){
    return Dmul_VsubV(c,v);
  }

  friend Vector operator * (const Vector & p, double c){
    return c*p;
  }
  friend double Distance(const Vector & v1,const Vector & v2){
    double dx=v1._x-v2._x;
    double dy=v1._y-v2._y;
  double dz=v1._z-v2._z;
    
    return sqrt(dx*dx+dy*dy+dz*dz);
  }
  friend double scalprod(const Vector & v1, const Vector & v2)
  {
    return v1._x*v2._x + v1._y*v2._y + v1._z*v2._z;
  }
  friend double scalprod(const VsubV & vsv, const Vector & v2)
  {
    return (vsv.v1._x-vsv.v2._x)*v2._x + (vsv.v1._y-vsv.v2._y)*v2._y + (vsv.v1._z-vsv.v2._z)*v2._z;
  }

  friend Vector vecprod(const Vector & v1, const Vector & v2){
    return Vector(v1._y*v2._z-v1._z*v2._y,v1._z*v2._x-v1._x*v2._z,v1._x*v2._y-v1._y*v2._x);
  }

  friend void addmul(Vector & res,double d, const Vector & v){
    res._x+=d*v._x;
    res._y+=d*v._y;
	res._z+=d*v._z;
  }
public:
  explicit Vector(double x=0,double y=0,double z=0): _x(x), _y(y), _z(z){};


  Vector(const VaddV & vv){
    _x=vv.v1._x+vv.v2._x;
    _y=vv.v1._y+vv.v2._y;
    _z=vv.v1._z+vv.v2._z;
  }
  Vector(const VsubV & vv){
    _x=vv.v1._x-vv.v2._x;
    _y=vv.v1._y-vv.v2._y;
    _z=vv.v1._z-vv.v2._z;
  }
  Vector(const DmulV & dv){
    _x=dv.c*dv.v._x;
    _y=dv.c*dv.v._y;
    _z=dv.c*dv.v._z;
  }
  
  Vector(const Dmul_VaddV & dvv){
    _x=dvv.c*(dvv.v.v1._x+dvv.v.v2._x);
    _y=dvv.c*(dvv.v.v1._y+dvv.v.v2._y);
    _z=dvv.c*(dvv.v.v1._z+dvv.v.v2._z);
  }
  Vector(const DmulV_addV & dvv){
    _x=dvv.dv.c*dvv.dv.v._x+dvv.v._x;
    _y=dvv.dv.c*dvv.dv.v._y+dvv.v._y;
    _z=dvv.dv.c*dvv.dv.v._z+dvv.v._z;
  }
  Vector(const Dmul_VsubV & dvv){
    _x=dvv.c*(dvv.v.v1._x-dvv.v.v2._x);
    _y=dvv.c*(dvv.v.v1._y-dvv.v.v2._y);
    _z=dvv.c*(dvv.v.v1._z-dvv.v.v2._z);
  }

  Vector(const DmulV_sub_DmulV & dvsdv){
    _x=dvsdv.dv1.c*dvsdv.dv1.v._x - dvsdv.dv2.c*dvsdv.dv2.v._x;
    _y=dvsdv.dv1.c*dvsdv.dv1.v._y - dvsdv.dv2.c*dvsdv.dv2.v._y;
    _z=dvsdv.dv1.c*dvsdv.dv1.v._z - dvsdv.dv2.c*dvsdv.dv2.v._z;
  }
  Vector(const DmulV_add_DmulV & dvadv){
    _x=dvadv.dv1.c*dvadv.dv1.v._x + dvadv.dv2.c*dvadv.dv2.v._x;
    _y=dvadv.dv1.c*dvadv.dv1.v._y + dvadv.dv2.c*dvadv.dv2.v._y;
    _z=dvadv.dv1.c*dvadv.dv1.v._z + dvadv.dv2.c*dvadv.dv2.v._z;
  }

  double & x() {return _x;}
  double x() const {return _x;}
  double & y() {return _y;}
  double y() const {return _y;}
  double & z() {return _z;}
  double z() const {return _z;}

  double norm(){
	  return sqrt(_x*_x+_y*_y+_z*_z);
  }

  const Vector & operator += (const Vector & p){
    _x+=p._x; _y+=p._y; _z+=p._z;
    return *this;
  }
  const Vector & operator += (const DmulV & dv){
    _x+=dv.c*dv.v._x; _y+=dv.c*dv.v._y; _z+=dv.c*dv.v._z;
    return *this;
  }
  const Vector & operator -= (const Vector & p){
    _x-=p._x; _y-=p._y; _z-=p._z;
    return *this;
  }
  const Vector & operator -= (const DmulV & dv){
    _x-=dv.c*dv.v._x; _y-=dv.c*dv.v._y; _z-=dv.c*dv.v._z;
    return *this;
  }
  const Vector & operator *= (double c){
    _x*=c; _y*=c; _z*=c;
    return *this;
  }

  /*
  void normalize(const Vector & size){
    if(_x<-size.x()/2) _x+=size.x();    
    if(_x> size.x()/2) _x-=size.x();

    if(_y<-size.y()/2) _y+=size.y();    
    if(_y> size.y()/2) _y-=size.y();

    if(_z<-size.z()/2) _z+=size.z();    
    if(_z> size.z()/2) _z-=size.z();

  }
  */

private:
	friend class boost::serialization::access;

    template<class Archive>
    void serialize(Archive & ar, const unsigned int version)
    {
        ar & _x;
        ar & _y;
        ar & _z;
    }

  double _x,_y,_z;
};

const Vector null(0,0,0);
//BOOST_CLASS_IMPLEMENTATION(Vector,object_serializable);
//BOOST_IS_MPI_DATATYPE(Vector);
//BOOST_CLASS_TRACKING(Vector,track_never);
//BOOST_IS_BITWISE_SERIALIZABLE(Vector);
#endif
