#ifndef HAMILTONIAN_HPP
#define HAMILTONIAN_HPP

#include <vector>

#include <boost/tuple/tuple.hpp>

#include "bcl.hpp"
#include "worms/weight.hpp"




template <int Spin>
class hamilonian {
  typedef bcl::markov<random01_t> markov_t;

public:
  hamiltonian(const double& _offset) : offset_(_offset) {}
  inline int n_types() const { return couplings_.size(); }

  void add_term(const heisenberg_operator<Spin>& op, const int& s0,
		const int& s1)
  {
    // Create weights and transition kernels if type is not registered yet
    double J = op.coupling;
    int type;
    if (std::find(couplings_.begin(), couplings_.end(), J) == couplings_.end())
      {
	type = couplings_.size();
	couplings_.push_back(J);
	weights_.push_back(weight<Spin>(op, offset));

	// Generate diagonal acceptance probabilities
	double lambda = weights_[type].max_diagonal_weight();
	std::vector<double> accept(1 << (2*n_bits<SPIN>::val));
	for (int c0 = 0; c0 < Spin; ++c0)
	  for (int c1 = 0; c1 < Spin; ++c1)
	    accept[spin_state::c2u<Spin>(c0, c1)] =
	      weights_[type][spin_state::c2p<Spin>(c0, c1, c0, c1)] / lambda;
	
	// Generate transition probabilities
	std::vector<markov_t> markov_raise_top;
	std::vector<markov_t> markov_raise_bottom;
	std::vector<std::vector<double> > raise_top_weights;
	std::vector<std::vector<double> > raise_bottom_weights;
	
	markov_raise_top.resize(operatorsize<Spin>::val);
	markov_raise_bottom.resize(operatorsize<Spin>::val);
	raise_top_weights.resize(operatorsize<Spin>::val);
	raise_bottom_weights.resize(operatorsize<Spin>::val);
	for (int c = 0; c < operatorsize<Spin>::val; ++c){
	  raise_top_weights[c].resize(4);
	  raise_bottom_weights[c].resize(4);

	  // Compute weights for raising worm
	  for (int ext_leg = 0; ext_leg < 4; ++ext_leg){
	    int direc = (ext_leg >> 1);
	    int state = c;

	    // Compute weights for raise top worms
	    if (direc==1){
	      if (spin_state::raise<Spin>(ext_leg, &state))
		raise_top_weights[c][ext_leg] = weights_[type][state];
	      else
		raise_top_weights[c][ext_leg] = 0;
	    }
	    else {
	      if (spin_state::lower<Spin>(ext_leg, &state))
		raise_top_weights[c][ext_leg] = weights_[type][state];
	      else
		raise_top_weights[c][ext_leg] = 0;
	    }
	    markov_raise_top[c] = markov_t(bcl::st2010(), raise_top_weights[c]);  
      

	    // Compute weights for raise bottom worm
	    state = c;
	    if (direc==1){
	      if (spin_state::lower<Spin>(ext_leg, &state))
		raise_bottom_weights[c][ext_leg] = weights_[type][state];
	      else
		raise_bottom_weights[c][ext_leg] = 0;
	    }
	    else {
	      if (spin_state::raise<Spin>(ext_leg, &state))
		raise_bottom_weights[c][ext_leg] = weights_[type][state];
	      else
		raise_bottom_weights[c][ext_leg] = 0;
	    }
	    markov_raise_bottom[c] = markov_t(bcl::st2010(), raise_bottom_weights[c]);  
	  }

	  markov_raise_top_.push_back(markov_raise_top);
	  markov_raise_bottom_.push_back(markov_raise_bottom);
	}
      }
    else
      type = std::find(couplings_.begin(), couplings_.end(), J) - couplings_.begin();
    
    sites_and_type_.push_back(boost::make_tuple(s0, s1, type));
  }
  
  template<class Engine>
  std::vector<boost::tuple<int, int, int> > random_operator(const Engine& eng) const
  { return sites_and_type_[sites_and_type_.size()*eng()]; }
  
    
private:
  const double offset_;   
      
  std::vector<double> couplings_;
  std::vector<boost::tuple<int, int, int> > sites_and_type_;
    
  std::vector<weight<Spin> > weights_;

  std::vector<std::vector<double> > accept_diag_;
  
  std::vector<std::vector<markov_t> > markov_raise_top_;
  std::vector<std::vector<markov_t> > markov_raise_bottom_;
  
};


#endif
