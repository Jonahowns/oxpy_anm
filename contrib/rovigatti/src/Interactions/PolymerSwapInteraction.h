#ifndef POLYMERSWAP_INTERACTION_H
#define POLYMERSWAP_INTERACTION_H

#include "Interactions/BaseInteraction.h"

#include <vector>
#include <array>

struct PSBond {
	BaseParticle *other;
	number energy;
	LR_vector force;
	LR_vector r;

	PSBond(BaseParticle *o, number e, LR_vector &nr) :
		other(o),
		energy(e),
		r(nr) {
			
	}
};

struct PSBondCompare {
	bool operator()(const PSBond &lhs, const PSBond &rhs) {
		if(lhs.other->index == rhs.other->index) {
			return false;
		}
		else {
			return true;
		}
	}
};

/**
 * @brief Handles interactions in microgel systems.
 *
 * This interaction is selected with
 * interaction_type = PolymerSwapInteraction
 *
 * @verbatim

 @endverbatim
 */
class PolymerSwapInteraction: public BaseInteraction<PolymerSwapInteraction> {
protected:
	std::array<number, 3> _Kfene = { {15., 15., 15.} };
	std::array<number, 3> _rfene = { {1.5, 1.5, 1.5} };
	std::array<number, 3> _sqr_rfene = { {2.25, 2.25, 2.25} };
	std::array<number, 3> _WCA_sigma = { {1.0, 1.0, 1.0} };
	std::array<number, 3> _PS_sqr_rep_rcut;

	number _PS_alpha = 0.;
	number _PS_beta = 0.;
	number _PS_gamma = 0.;
	int _PS_n = 6;

	/// Three-body potential stuff
	number _3b_rcut = -1.;
	number _sqr_3b_rcut = -1.;
	number _3b_sigma = 1.05;
	number _3b_range = 1.6;
	number _3b_lambda = 1.0;
	number _3b_epsilon = 10.;
	number _3b_prefactor = 0.1;
	number _3b_A_part = 0.;
	number _3b_B_part = 0.;

	std::vector<LR_vector> _inter_chain_forces;

	std::map<int, std::set<PSBond, PSBondCompare> > _bonds;

	std::string _bond_filename;

	bool _only_links_in_bondfile = true;
	int _chain_size = -1;
	int _N_chains = -1;
	number _T = 0.;

	void _update_inter_chain_forces(BaseParticle *p, BaseParticle *q, LR_vector p_force, const LR_vector &pq_r);

	number _fene(BaseParticle *p, BaseParticle *q, bool update_forces);
	number _WCA(BaseParticle *p, BaseParticle *q, bool update_forces);
	number _sticky(BaseParticle *p, BaseParticle *q, bool update_forces);
	number _three_body(BaseParticle *p, PSBond &new_bond, bool update_forces);

public:
	enum {
		BONDED = 0, NONBONDED = 1
	};
	enum {
		MONOMER = 0, STICKY = 1
	};

	bool no_three_body = false;

	PolymerSwapInteraction();
	virtual ~PolymerSwapInteraction();

	virtual void get_settings(input_file &inp);
	virtual void init();

	number P_inter_chain();

	virtual void allocate_particles(std::vector<BaseParticle *> &particles);

	void begin_energy_computation() override;

	virtual number pair_interaction(BaseParticle *p, BaseParticle *q, bool compute_r = true, bool update_forces = false);
	virtual number pair_interaction_bonded(BaseParticle *p, BaseParticle *q, bool compute_r = true, bool update_forces = false);
	virtual number pair_interaction_nonbonded(BaseParticle *p, BaseParticle *q, bool compute_r = true, bool update_forces = false);
	virtual number pair_interaction_term(int name, BaseParticle *p, BaseParticle *q, bool compute_r = true, bool update_forces = false) {
		return _pair_interaction_term_wrapper(this, name, p, q, compute_r, update_forces);
	}

	virtual void read_topology(int *N_stars, std::vector<BaseParticle *> &particles);
	virtual void check_input_sanity(std::vector<BaseParticle *> &particles);
};

extern "C" PolymerSwapInteraction *make_PolymerSwapInteraction();

#endif /* POLYMERSWAP_INTERACTION_H */
