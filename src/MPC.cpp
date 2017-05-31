#include "MPC.h"
#include <cppad/cppad.hpp>
#include <cppad/ipopt/solve.hpp>
#include "Eigen-3.3/Eigen/Core"

using CppAD::AD;

// This value assumes the model presented in the classroom is used.
//
// It was obtained by measuring the radius formed by running the vehicle in the
// simulator around in a circle with a constant steering angle and velocity on a
// flat terrain.
//
// Lf was tuned until the the radius formed by the simulating the model
// presented in the classroom matched the previous radius.
//
// This is the length from front to CoG that has a similar radius.
const double Lf = 2.67;

// Both the reference cross track and orientation errors are 0.
// The reference velocity is set to 100 mph.
double ref_cte = 0;
double ref_epsi = 0;
double v_max = 100;
double ref_psi = 0.1;
double psi_max = 0;

// The solver takes all the state variables and actuator
// variables in a singular vector. Thus, we should to establish
// when one variable starts and another ends to make our lifes easier.
size_t x_start = 0;
size_t y_start = x_start + N;
size_t psi_start = y_start + N;
size_t v_start = psi_start + N;
size_t cte_start = v_start + N;
size_t epsi_start = cte_start + N;
size_t delta_start = epsi_start + N;
size_t a_start = delta_start + N - 1;
int debug = 0;

class FG_eval {
 public:
  // Fitted polynomial coefficients
  Eigen::VectorXd coeffs;
  FG_eval(Eigen::VectorXd coeffs) { this->coeffs = coeffs; }

  typedef CPPAD_TESTVECTOR(AD<double>) ADvector;
  void operator()(ADvector& fg, const ADvector& vars) {
    
    //cout << "Vars = " << vars << endl;
    // The cost is stored is the first element of `fg`.
    // Any additions to the cost should be added to `fg[0]`.
    fg[0] = 0;
    
    // The part of the cost based on the reference state.
    for (int i = 0; i < N; i++) {
      CppAD::AD<double> cte_cost = 1000 * CppAD::pow((vars[cte_start + i] - ref_cte), 2);
      CppAD::AD<double> epsi_cost = 10000 * CppAD::pow((vars[epsi_start + i] - ref_epsi), 2);
      CppAD::AD<double> v_cost = CppAD::pow((vars[v_start + i] - v_max), 2);
      
      if(i == 0 && debug)
        cout << "Cost for i = " << i << ", " << cte_cost << ", " \
                << epsi_cost <<", " << v_cost <<", ";
      fg[0] += cte_cost + epsi_cost + v_cost;
    
    
    // Minimize the use of actuators.
      if( i < N-1) {
        CppAD::AD<double> delta_cost = 100*CppAD::pow(vars[delta_start + i], 2);
        CppAD::AD<double> a_cost = 10*CppAD::pow(vars[a_start + i], 2);
        fg[0] += delta_cost + a_cost;
        if(i == 0 && debug)
          cout << delta_cost << ", " << a_cost << ", ";
      }
      if(i < N-2) {
        // Minimize the value gap between sequential actuations.

        CppAD::AD<double> diff_delta_cost = CppAD::pow(vars[delta_start + i + 1] - vars[delta_start + i], 2);
        CppAD::AD<double> diff_a_cost = 10*CppAD::pow(vars[a_start + i + 1] - vars[a_start + i], 2);
        fg[0] += diff_delta_cost + diff_a_cost;
        if(i == 0 && debug)
          cout << diff_delta_cost << ", " << diff_a_cost;
      }
      if(i == 0 && debug)
        cout << endl;
    }
    //cout << "fg[0] = " << fg[0] << endl;
    //
    // Setup Constraints
    //
    // NOTE: In this section you'll setup the model constraints.
    
    // Initial constraints
    //
    // We add 1 to each of the starting indices due to cost being located at
    // index 0 of `fg`.
    // This bumps up the position of all the other values.
    fg[1 + x_start] = vars[x_start];
    fg[1 + y_start] = vars[y_start];
    fg[1 + psi_start] = vars[psi_start];
    fg[1 + v_start] = vars[v_start];
    fg[1 + cte_start] = vars[cte_start];
    fg[1 + epsi_start] = vars[epsi_start];
    //cout << "fg[1] = " << fg[1] << endl;
    // The rest of the constraints
    for (int i = 0; i < N - 1; i++) {
      // The state at time t+1 .
      AD<double> x1 = vars[x_start + i + 1];
      AD<double> y1 = vars[y_start + i + 1];
      AD<double> psi1 = vars[psi_start + i + 1];
      AD<double> v1 = vars[v_start + i + 1];
      AD<double> cte1 = vars[cte_start + i + 1];
      AD<double> epsi1 = vars[epsi_start + i + 1];
      
      // The state at time t.
      AD<double> x0 = vars[x_start + i];
      AD<double> y0 = vars[y_start + i];
      AD<double> psi0 = vars[psi_start + i];
      AD<double> v0 = vars[v_start + i];
      AD<double> cte0 = vars[cte_start + i];
      AD<double> epsi0 = vars[epsi_start + i];
      
      // Only consider the actuation at time t.
      AD<double> delta0 = vars[delta_start + i];
      AD<double> a0 = vars[a_start + i];
      //value of y based on polynomial
      AD<double> f0 = coeffs[0] + coeffs[1] * x0 + \
                    coeffs[2] * x0 * x0 + coeffs[3] * x0 * x0 * x0;
      //value of psi based on detivative of poly.
      AD<double> psides0 = CppAD::atan(3.0 * coeffs[3] * x0 * x0 + \
                                        2.0 * coeffs[2] * x0 + coeffs[1]);
      
      // Here's `x` to get you started.
      // The idea here is to constraint this value to be 0.
      //
      // Recall the equations for the model:
      // x_[t+1] = x[t] + v[t] * cos(psi[t]) * dt
      // y_[t+1] = y[t] + v[t] * sin(psi[t]) * dt
      // psi_[t+1] = psi[t] + v[t] / Lf * delta[t] * dt
      // v_[t+1] = v[t] + a[t] * dt
      // cte[t+1] = f(x[t]) - y[t] + v[t] * sin(epsi[t]) * dt
      // epsi[t+1] = psi[t] - psides[t] + v[t] * delta[t] / Lf * dt
      fg[2 + x_start + i] = x1 - (x0 + v0 * CppAD::cos(psi0) * dt);
      fg[2 + y_start + i] = y1 - (y0 + v0 * CppAD::sin(psi0) * dt);
      fg[2 + psi_start + i] = psi1 - (psi0 + v0 * (-delta0) / Lf * dt);
      fg[2 + v_start + i] = v1 - (v0 + a0 * dt);
      fg[2 + cte_start + i] =
      cte1 - ((f0 - y0) + (v0 * CppAD::sin(epsi0) * dt));
      fg[2 + epsi_start + i] =
                  epsi1 - ((psi0 - psides0) + v0 * (-delta0) / Lf * dt);
    }
    //cout << "max x " << vars[y_start - 1] << endl;
  }
};

//
// MPC class definition implementation.
//
MPC::MPC() {}
MPC::~MPC() {}

vector<double> MPC::Solve(Eigen::VectorXd state, Eigen::VectorXd coeffs) {
  bool ok = true;
  size_t i;
  typedef CPPAD_TESTVECTOR(double) Dvector;
  // As actuator has 100ms latency we set the initial variables to where they will be in 100ms.
  const double latency = 0.1;// in seconds
  double x = state[0];
  double y = state[1];
  double psi = state[2];
  double v = state[3];
  double steering_angle = state[4];
  double throttle = state[5];
  double x_max = state[6];
  double y_minmax = state[7];
  
  //The maximum distance we are safe to get to is 1/2 of x_max as we have no waypoints
  // after x_max and we want to keep some buffer (50% for now)
  // so the max speed we can get to is the speed we can reach in N*dt seconds?
  //we don't know the unit of x is it feet?
  //v_max = x_max/(N*latency);
  //cout << "v_max = " << v_max << endl;
  
  double cte = coeffs[0]; // current CTE is fitted polynomial (road curve) evaluated at px = 0.0
  // epsi = arctan(f') where f' is the derivative of the fitted polynomial
  // f' = 3.0 * K[3] * px0 * px0 + 2.0 * K[2] * px0 + K[1]
  double epsi = -atan(coeffs[1]);
  //cout << "cte, epsi = " << cte << ", " << epsi << endl;
  // TODO: Set the number of model variables (includes both states and inputs).
  // For example: If the state is a 4 element vector, the actuators is a 2
  // element vector and there are 10 timesteps. The number of variables is:
  //
  // 4 * 10 + 2 * 9
  // number of independent variables
  // N timesteps == N - 1 actuations
  size_t n_vars = N * 6 + (N - 1) * 2;
  // Number of constraints
  size_t n_constraints = N * 6;

  // Initial value of the independent variables.
  // SHOULD BE 0 besides initial state.
  Dvector vars(n_vars);
  for (int i = 0; i < n_vars; i++) {
    vars[i] = 0;
  }
  for (int i = y_start; i < n_vars; i++) {
    vars[i] = 0;
  }
  // Set the initial variable values -
  Dvector vars_lowerbound(n_vars);
  Dvector vars_upperbound(n_vars);
  // TODO: Set lower and upper limits for variables.
  //upper & lower bound for x is from 0 to max waypoint
  for (int i = x_start; i < x_start + N; i++) {
    vars_lowerbound[i] = 0;
    vars_upperbound[i] = x_max;
  }
  
  // Set all non-actuators upper and lowerlimits
  // to the max negative and positive values.
  //cout << "y_minmax = " << y_minmax << endl;
  for (int i = y_start; i < y_start + N; i++) {
    if(y_minmax < 0) {
      vars_lowerbound[i] = -1.0e19; //y_minmax;
      vars_upperbound[i] = 1.0e19; //-y_minmax;
      //cout << "lt " << vars_lowerbound[i] << ", " << vars_upperbound[i] << endl;
    } else {
      vars_lowerbound[i] = -1.0e19; //-y_minmax;
      vars_upperbound[i] = 1.0e19; //y_minmax;
      //cout << vars_lowerbound[i] << ", " << vars_upperbound[i] << endl;
    }
  }
  for (int i = psi_start; i < delta_start; i++) {
    vars_lowerbound[i] = -1.0e19;
    vars_upperbound[i] = 1.0e19;
  }
  // The upper and lower limits of delta are set to -25 and 25
  // degrees (values in radians).
  // NOTE: Feel free to change this to something else.
  for (int i = delta_start; i < a_start; i++) {
    vars_lowerbound[i] = -0.436332;
    vars_upperbound[i] = 0.436332;
  }

  // Acceleration/decceleration upper and lower limits.
  // NOTE: Feel free to change this to something else.
  for (int i = a_start; i < n_vars; i++) {
    vars_lowerbound[i] = -0.5; //If we don't want to break hard then we can reduce this.
    //upper limit of throttle adjusted based inversely with speed.
    vars_upperbound[i] = 1-v/v_max;
  }
  
  // Lower and upper limits for the constraints
  // Should be 0 besides initial state.
  Dvector constraints_lowerbound(n_constraints);
  Dvector constraints_upperbound(n_constraints);
  for (int i = 0; i < n_constraints; i++) {
    constraints_lowerbound[i] = 0;
    constraints_upperbound[i] = 0;
  }
  //We should take delay into account here.
 
  vars[x_start] = 0 + v*latency; //0 in vehicle coordinates
  //vars[y_start -1] = x_max;
  vars[y_start] = 0; //0 in vehicle coordinates
  vars[psi_start] = 0 + v * (-steering_angle) / Lf * latency; //psi is 0 in vehicle coordinates
  vars[v_start] = v + throttle*latency;
  vars[cte_start] = cte + v * sin(epsi) * dt;
  vars[epsi_start] = epsi + v * (-steering_angle) / Lf * dt;

  constraints_lowerbound[x_start] = vars[x_start];
  constraints_lowerbound[y_start] = vars[y_start];
  constraints_lowerbound[psi_start] = vars[psi_start];
  constraints_lowerbound[v_start] = vars[v_start];
  constraints_lowerbound[cte_start] = vars[cte_start];
  constraints_lowerbound[epsi_start] = vars[epsi_start];
  
  constraints_upperbound[x_start] = vars[x_start];
  constraints_upperbound[y_start] = vars[y_start];
  constraints_upperbound[psi_start] = vars[psi_start];
  constraints_upperbound[v_start] = vars[v_start];
  constraints_upperbound[cte_start] = vars[cte_start];
  constraints_upperbound[epsi_start] = vars[epsi_start];
  
  // object that computes objective and constraints
  FG_eval fg_eval(coeffs);

  //
  // NOTE: You don't have to worry about these options
  //
  // options for IPOPT solver
  std::string options;
  // Uncomment this if you'd like more print information
  options += "Integer print_level  0\n";
  // NOTE: Setting sparse to true allows the solver to take advantage
  // of sparse routines, this makes the computation MUCH FASTER. If you
  // can uncomment 1 of these and see if it makes a difference or not but
  // if you uncomment both the computation time should go up in orders of
  // magnitude.
  options += "Sparse  true        forward\n";
  options += "Sparse  true        reverse\n";
  // NOTE: Currently the solver has a maximum time limit of 0.5 seconds.
  // Change this as you see fit.
  options += "Numeric max_cpu_time          0.5\n";

  // place to return solution
  CppAD::ipopt::solve_result<Dvector> solution;
  
  //cout << "options = " << options << endl;
  //cout << "vars=" << vars << endl;
  //cout << "vars_lowerbound = " << vars_lowerbound << endl;
  //cout << "vars_upperbound = " << vars_upperbound << endl;
  //cout << "constraints_lowerbound = " << constraints_lowerbound << endl;
  //cout << "constraints_upperbound = " << constraints_upperbound << endl;
  // solve the problem
  CppAD::ipopt::solve<Dvector, FG_eval>(
      options, vars, vars_lowerbound, vars_upperbound, constraints_lowerbound,
      constraints_upperbound, fg_eval, solution);
  
  // Cost
  auto cost = solution.obj_value;
  
  // Check some of the solution values
  ok &= solution.status == CppAD::ipopt::solve_result<Dvector>::success;
 /* if (ok) {
    std::cout << "OK! Cost:" << cost << std::endl;
  } else {
    std::cout << "SOMETHING IS WRONG!" << cost << std::endl;
  }*/

  //std::cout <<  "Cost " << cost << std::endl;

  // TODO: Return the first actuator values. The variables can be accessed with
  // `solution.x[i]`.
  //
  // {...} is shorthand for creating a vector, so auto x1 = {1.0,2.0}
  // creates a 2 element double vector.
  //cout << "mpc_y_vals = ";
  mpc_x_vals.resize(N);
  mpc_y_vals.resize(N);
  for (int i = 0; i < N; ++i) {
    mpc_x_vals[i] = solution.x[x_start + i];
    mpc_y_vals[i] = solution.x[y_start + i];
    //cout << mpc_y_vals[i] << ", " ;
  }
  //cout << endl;
  
  //cout << "solution: " << solution.x[delta_start] << "," << solution.x[a_start] << endl;
  return {solution.x[delta_start], solution.x[a_start]};
}
