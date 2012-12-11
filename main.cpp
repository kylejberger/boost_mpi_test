//Outline for new DEM code
//Creation begun: 8/5/2011
//This program, when completed, will perform a DEM simulation of particles defined by the user
//using an input file.

//include directives here
#include <vector>
#include <boost/mpi.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include "Particle.h"
#include "Vector.h"
#include <boost/serialization/map.hpp>

//standard namespace and berger_DEM namespace
using namespace std;
namespace mpi=boost::mpi;



int main(int argc, char* argv[]) 
{
  //initialize the mpi communicator
	boost::mpi::environment env(argc, argv);
	boost::mpi::communicator world;
	
	vector<vector<Particle> >particles_to_be_sent;
	vector<Particle> particles_received;
	//vector<int> particles_received;
	//vector<vector<int> >particles_to_be_sent;
	particles_to_be_sent.resize(world.size());
	

	//send particles 
		Particle test; //test particle, will now fill all of its variables with something
		test.velocity()=null;//null is a global variable of the right type for these variables
		test.position()=null;
		test.omega()=null;
		test.global_part_num()=0;
		test.location_number()=1;
		test.tan_contact_hist().insert(pair<int,Vector>(1,null));
		int test_int=10;

		for(int i=0;i<world.size();i++)
		{
			particles_to_be_sent[i].assign(35,test); //number of particles sent is assigned here
			//particles_to_be_sent[i].assign(100000,test_int);
		}

		for(int i=0;i<world.size();i++)
		{
			world.isend(i,i,particles_to_be_sent[i]);
		}

		//receive particles
		int receive_counter=0;
		//boost::optional<boost::mpi::status> msginfo=world.iprobe(mpi::any_source, world.rank());
		while (receive_counter<world.size())
		{
			//if (msginfo)
			//{
				cout<<"rank "<<world.rank()<<" receiving message from rank "<<receive_counter<<endl;
				world.recv(receive_counter,world.rank(),particles_received);
				cout<<"message to rank "<<world.rank()<<" from rank "<<receive_counter<<" is being analyzed: consists of "<<particles_received.size()<<" particles"<<endl;
				cout<<"message to rank "<<world.rank()<<" from rank "<<receive_counter<<" successfully received"<<endl;
				receive_counter++;
			//}
			//msginfo=world.iprobe(mpi::any_source, world.rank());
		}
	return 0;
}
