/*Header file containing the particle class
default copy works
This class is used to represent a particle in the DEM simulations. It includes the following methods for manipulating them:
Constructors:
  Particle()->takes no inputs and initializes the particles with arbitrary parameters and position/velocity/higher order
		derivatives all to (0,0,0). This should only be used for initial declaration in the case of a vector. Before using
		the particles, all of the parameters and vectors should be set individually.
	Particle(Vector& velocity_in, Vector& position_in, double rho_in, double radius_in, Vector& omega_in, Vector& theta_in, int dimensions_in, int type_in) 
		-> takes several inputs and assigns appropriate vectors/parameters as follows:
			velocity_in(Vector) is the velocity of the particle
			position_in(Vector) is the position of the particle
			rho_in(double) is the density of the particle(the mass will automatically be assigned based on density and radius)
			radius_in(double) is the radius of the particle
			omega_in(Vector) is the initial rotational velocity of the particle
			theta_in(Vector) is the initial orientation of the particle
			dimensions_in is the number of dimensions in the system(currently does not affect anything)
			type_in is the type of particle(0 is a normal particle and 1 is a stationary wall particle)
			the other time derivatives of the position and orientation are automatically intialized as 0(support to change them will be added later)

Methods for accessing or directly changing particle properties
	methods to changing information about the particle
		Vector & position()-> used to change the particle's position
		Vector & velocity()-> used to change the particle's velocity
		Vector & omega()-> used to change the particle's angular velocity
		Vector & theta()-> used to change the particle's orientation
		double & rho()-> used to change the particle's density
		double & radius()-> used to change the particle's radius
		double & mass()-> used to change the particle's mass
		unsigned int & type()-> used to change the type of particle

	methods for returning information about the particle
		Vector position() const-> returns the particle's position
		Vector velocity() const-> return the particle's velocity
		Vector omega() const-> returns the particle's angular velocity
		Vector theta() const-> returns the particle's orientation
		double rho() const-> returns the particle's density
		double radius() const-> returns the particle's radius
		double mass() const-> returns the particle's mass
		unsigned int type() const->  returns the type of the particle	

Other methods:
	double rho()-> returns the density(double) of the particle and can also be used to assign it a new value
	double radius()-> returns the radius(double) of the particle and can also be used for assignment
	double type()-> returns the type(int) of the particle and can also be used for assignment
	double mass()-> returns the mass(double) of the particle and can also be used for assignment(NOTE: the mass will NOT automatically be updated if the radius or density change)
	Vector position()-> returns the position(Vector) of the particle
	Vector velocity()-> returns the velocity(Vector) of the particle
	Vector omega()-> returns the angular velocity(Vector) of the particle
	Vector theta()-> returns the orientation(Vector) of the particle
	void move(Vector&  incrpos)-> changes the position by adding the vector incrpos to the current position vector
	void accelerate(Vector& incrvel)-> changes the velocity by adding the vector incrvel to the current position vector
	void accelerate_omega(Vector& incromega)-> changes the angular velocity by adding the vector incromega to the current angular velocity vector
	void change_orientation(Vector& incrtheta)-> changes the orientation by adding the vector incrtheta to the current orientation vector
	void predict(User_input&)->takes in a User_input type and predicts the new positions, velocities etc.(Gear algorithm) of the particle
	void add_force(const Vector& force_in)-> adds the vector force_in to the current force vector
	void correct(User_input&)-> takes in a User_input type and corrects the predicted positions, velocities etc.(Gear algorithm)
	double kinetic energy() const->outputs the kinetic energy(double) of the particle
	void set_force_to_zero()-> sets the force vector to (0,0,0)
	void periodic_bc(User_input &)-> takes in a User_input type and applies appropriate periodic conditions to the particle if necessary
	void tangential
	void plume_force(User_input &)-> computes the plume force exerted on the particle based on its height above the anchoring plane



friend functions:
	friend std::istream & operator >>(std::istream &, Particle &)->overloads the >> operator to take direct input from istream in the following order:
		xposition   yposition   zposition   radius    density    xvelocity    yvelocity    zvelocity
	friend std::ostream & operator <<(std::ostream &, Particle &)-> overloads the << operator to give output to ostream in same order as above
	friend std::ifstream & operator >>(std::ifstream &, Particle &)-> same as 2 functions up but with ifstream
	friend std::ofstream & operator <<(std::ofstream &, Particle &)-> same as 2 functions up but with ofstream
	friend double Distance(const Particle & p1, const Particle & p2, User_input & input_data)-> finds the distance between particles p1 and p2(accounts for periodicity)
	friend void force(Particle&, Particle&, User_input&)-> calculates the force between two particles and adds them to the force vector of each particle

Inline functions:
	double normalize(double dx, double L, bool periodic_wall)-> selects the appropriate image of a particle if periodic_wall is true(dx is the 
		non-periodic distance between them and L is the length of the box in the direction of dx)

Implementations of these functions can be found in Particle.cpp
*/

#ifndef Particle_h
#define Particle_h

#include <iostream>
#include <math.h>
#include <map>
#include <vector>
#include <boost/mpi.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <boost/serialization/map.hpp>

#include "Vector.h"

			//class User_input;
			//class data_tracking;

			inline double normalize(double dx, double L, bool periodic_wall)
			{
				//selects the appropirate periodic image of a particle for calculating distances
				if(periodic_wall)
				{
					while(dx<-L/2) dx+=L;
					while(dx>=L/2) dx-=L;
				}
				return dx;
			}


			class Particle 
			{

				/*friend double Distance(const Particle & p1, const Particle & p2, User_input & input_data)
				{
					//make a call to normalize(for each direction) in case periodic image is necessary
					double dx=normalize(p1._position.x() - p2._position.x(),input_data.xlength(),input_data.periodic_walls_x());
					double dy=normalize(p1._position.y() - p2._position.y(),input_data.ylength(),input_data.periodic_walls_y());
					double dz=normalize(p1._position.z() - p2._position.z(),input_data.zlength(),input_data.periodic_walls_z());
					//return the distance between particles
					return sqrt(dx*dx+dy*dy+dz*dz);
				}*/

				//friend void force(Particle &, Particle &, User_input &, data_tracking &, std::vector<Particle_props> &, std::vector<Particle_forces> &, std::vector<int> &, const int);

				private:
				friend class boost::serialization::access;

				template<class Archive>
				void serialize(Archive & ar, const unsigned int version)
				{
					ar & _velocity;
					ar & _position;
					ar & _omega;
					ar & _global_part_num;
					ar & _tan_contact_hist;
					ar & _location_number;
				}
				Vector _velocity, _position, _omega;
				//velocity is a Vector representing the 
				//position is an array representing positions in x, y, and z in m
				//omega is the angular velocity of the particle in x, y, and z in 1/s
				//theta is the orientation of the particle in x, y, and z
				int _global_part_num, _location_number;
				std::map<int,Vector> _tan_contact_hist;
				//map<int,double> _last_contact_time;
		
				public:
				//constructors
				Particle(): _velocity(null), _position(null), _omega(null), _global_part_num(0), _location_number(-1) {};
		
				//methods to set particle properties
				Vector & position() {return _position;} //used to change the particle's position
				Vector & velocity() {return _velocity;} //return the particle's velocity
				Vector & omega() {return _omega;} //used to change the particle's angular velocity
				int & global_part_num() {return _global_part_num;}
				int & location_number() {return _location_number;}
				std::map<int,Vector> & tan_contact_hist() {return _tan_contact_hist;}
		
				//methods necessary for algorithm
				/*void read_particle(std::ifstream &);
				void correct(User_input & input_data, std::vector<Particle_props> &, std::vector<Particle_forces> &);
				double kinetic_energy(double) const;
				double rotational_energy(double) const;
				double momentum_x(double) const;
				double momentum_y(double) const;
				double momentum_z(double) const;
				void wall_collision(User_input &, std::string &, std::vector<Particle_props> &, std::vector<Particle_forces> &);
				void periodic_bc(User_input &);
				void tangential_displacement(const int &, Vector &, const Vector &, const Vector &, User_input &, std::map<int,Vector>&);
				void tangential_displacement_wall(const int &, Vector &, const Vector &, const Vector &, User_input &);
				void plume_force(User_input &, std::vector<Particle_props> &, std::vector<Particle_forces> &);*/
			};


//BOOST_CLASS_IMPLEMENTATION(Particle,object_serializable);
//BOOST_IS_MPI_DATATYPE(Particle);
//BOOST_CLASS_TRACKING(Particle,track_never);
#endif
