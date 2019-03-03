#include <iostream>
#include <iomanip>
#include <math.h>
#include <plstream.h>
#include <thread>
#include <random>
#include <stdexcept>
#include <string>



/*
  Gravitational force vector on body i due to body j is
------------------------------------------------------------------------------------
    F_ij  =  C  *  m_i  *  m_j  *  ( pos_j  -  pos_i )  /  abs( pos_j  -  pos_i )^3
------------------------------------------------------------------------------------
  where C is the gravitational constant (a parameter in this program :), 
  m are masses, and pos are position vectors. 
  ================================================================================
  The force on body i 
------------------------------------------------------------------------------------
    F  =  m_i  *  acc_i  =  m_i  *  d vel_i  /  d t
------------------------------------------------------------------------------------
  is equal to the sum over force vectors due to other bodies, i.e. the 
  N-body problem involves the coupled ordinary differential equations
------------------------------------------------------------------------------------
    d vel_i  /  d t  =  sum( F_ij ) for i!=j  /  m_i
    d pos_i  /  d t  =  vel_i
  =================================================================================
*/

double const C=1.01;


class body {
  public:
    body () : m{rm(1.0)}, x{rp(0.0)}, y{rp(0.0)}, v{rv(0.0)}, w{rv(0.0)}, f{0.0}, g{0.0} {}
    double m, x, y, v, w, f, g;
  private:
    double rm (double offset_from_zero) { // random mass
      std::random_device rd; // random number from hardware
      std::mt19937 eng(rd()); // seed generator
      std::weibull_distribution<> distr(1,2); // distribution parameters
      return offset_from_zero+distr(eng);
    } 
    double rp (double mean) { // random position
      std::random_device rd;  // random number from hardware
      std::mt19937 eng(rd()); // seed generator
      std::normal_distribution<> distr(0.0, 10.0); // distribution parameters
      return mean+distr(eng);
    } 
    double rv (double mean) { // random velocity
      std::random_device rd;  // random number from hardware
      std::mt19937 eng(rd()); // seed generator
      std::normal_distribution<> distr(0.0, 0.5); // distribution parameters
      return mean+distr(eng);
    } 
};


void input_verify_stoi(std::vector<std::string> vec, int *verified) {
  for (unsigned i = 1; i < vec.size(); i++) {
    auto arg = vec[i];
    int x;
    try {
      std::size_t pos;
      verified[i] = std::stoi(arg, &pos);
      if (pos < arg.size()) {
        std::cerr << "Trailing characters after number: " << arg << '\n';
      }
    } catch (std::invalid_argument const &ex) {
      std::cerr << "Invalid number: " << arg << '\n';
    } catch (std::out_of_range const &ex) {
      std::cerr << "Number out of range: " << arg << '\n';
    }
  }
}


body fsum ( body computed, body affecting ) {
  double xvec = affecting.x - computed.x;
  double yvec = affecting.y - computed.y;
  double bodyd = std::sqrt(xvec*xvec + yvec*yvec);
  double fscal = C * computed.m * affecting.m / (bodyd * bodyd * bodyd);    // scalar part of force
  computed.f += fscal * xvec;
  computed.g += fscal * yvec;
  return computed;
}; 


body integrate ( body computed, double dt ) {  // Euler time-integration 
  computed.v += dt * computed.f / computed.m;
  computed.w += dt * computed.g / computed.m;
  computed.x += dt * computed.v;
  computed.y += dt * computed.w;
  return computed;
}

void print_output(body b)
{
  std::cout.precision(3);
  std::cout << "  m = " << std::setw(4) << b.m <<
               "  x = " << std::setw(6) << b.x <<
               "  y = " << std::setw(6) << b.y <<
               "  v = " << std::setw(6) << b.v << 
               "  w = " << std::setw(6) << b.w << "\n";
}


int main (int argc, char **argv) {


  if (argc < 2) {
    std::cout << "\n" <<
	         " ----------------------------------------------------------------" << "\n" <<
	         "  Optionally specify arguments, " << "\n" << 
	         "    ./main <number-of-bodies>" << "\n" <<
		 "  or" << "\n"
		 "    ./main <number-of-bodies> <number-of-iterations>" << "\n"
		 "  or" << "\n"
		 "    ./main <number-of-bodies> <number-of-iterations> <time-step>" << "\n"
	         " ----------------------------------------------------------------" << "\n" <<
		 "\n";
  }


  std::vector<std::string> inputs(argv, argv + argc);
  int input[argc]; 
  input_verify_stoi(inputs, input); 

  int num_bodies = (argc > 1) ? input[1] : 10;
  int iterations = (argc > 2) ? input[2] : 100;
  int delta_t    = (argc > 3) ? input[3] : 1.0;
  

  std::cout << "\n" <<
	       "            Bodies: " << num_bodies << "\n" <<
               "        Iterations: " << iterations << "\n" <<
               "  Integration step: " << delta_t << "\n" <<
	       "\n";


  body bodies[num_bodies];
  
  
  PLFLT xmin{-75}, ymin{-75}, xmax={75}, ymax={75},
        x[6] = {-75, -45, -15, 15, 45, 75}, 
        y[6] = {-75, -45, -15, 15, 45, 75};
  PLINT just=0, axis=0;
  PLFLT xcoords[num_bodies], ycoords[num_bodies];
  plstream *pls;

  
  pls = new plstream();
//  plscolbg(255, 255, 255);      // uncomment to make white background
  pls->init();           // start plplot object
//  plscol0(0, 0, 0, 0);          // uncomment to make white background
  pls->env(xmin, xmax, ymin, ymax, just, axis );
  pls->lab( "x", "y", "N-body problem");
 


// Calculate sum of forces
  for (int timestep = 0; timestep < iterations; timestep++) {
    for (int i = 0; i < num_bodies; i++) {
      bodies[i].f = 0.0;
      bodies[i].g = 0.0;
      for (int j = 0; j < num_bodies; j++) {
        if ( j == i ) continue;
        bodies[i] = fsum( bodies[i], bodies[j] );
      }
    }

// Solve ODEs
    for (int i = 0; i < num_bodies; i++) {
      bodies[i] = integrate( bodies[i], delta_t );
      xcoords[i] = bodies[i].x;
      ycoords[i] = bodies[i].y;
    }
    for (int i = 0; i < num_bodies; i++) {
      print_output( bodies[i] );
      pls->poin( num_bodies , xcoords, ycoords, 1);
      std::this_thread::sleep_for(std::chrono::milliseconds(100));
    }
  }
  delete pls; 
}
