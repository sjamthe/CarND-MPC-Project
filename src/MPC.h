#ifndef MPC_H
#define MPC_H

#include <vector>
#include "Eigen-3.3/Eigen/Core"

using namespace std;

const size_t N = 15; //Predict 10 steps in future
const double dt = 0.08; // predict at duration of 0.1 secs

class MPC {

  public:
    vector<double> mpc_x_vals;
    vector<double> mpc_y_vals;
  MPC();

  virtual ~MPC();

  // Solve the model given an initial state and polynomial coefficients.
  // Return the first actuatotions.
  vector<double> Solve(Eigen::VectorXd state, Eigen::VectorXd coeffs);
};

#endif /* MPC_H */
