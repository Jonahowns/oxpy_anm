#include "PolymerSwapInteraction.h"

#include "Particles/CustomParticle.h"
#include <fstream>

PolymerSwapInteraction::PolymerSwapInteraction() :
				BaseInteraction<PolymerSwapInteraction>() {
	_int_map[BONDED] = &PolymerSwapInteraction::pair_interaction_bonded;
	_int_map[NONBONDED] = &PolymerSwapInteraction::pair_interaction_nonbonded;

}

PolymerSwapInteraction::~PolymerSwapInteraction() {

}

void PolymerSwapInteraction::get_settings(input_file &inp) {
	IBaseInteraction::get_settings(inp);

	getInputInt(&inp, "PS_n", &_PS_n, 0);
	getInputString(&inp, "PS_bond_file", _bond_filename, 1);
	getInputBool(&inp, "PS_only_links_in_bondfile", &_only_links_in_bondfile, 0);
	getInputInt(&inp, "PS_chain_size", &_chain_size, 1);
	
	getInputNumber(&inp, "PS_alpha", &_PS_alpha, 0);

	if(_PS_alpha != 0.) {
		throw oxDNAException("PolymerSwap: PS_alpha != 0 is not supported yet!");
	}

	getInputNumber(&inp, "PS_3b_sigma", &_3b_sigma, 0);
	getInputNumber(&inp, "PS_3b_range", &_3b_range, 0);
	getInputNumber(&inp, "PS_3b_lambda", &_3b_lambda, 0);
	getInputNumber(&inp, "PS_3b_epsilon", &_3b_epsilon, 0);

	getInputNumber(&inp, "T", &_T, 1);

	for(int i = 0; i < 3; i++) {
		std::string key = Utils::sformat("PS_rfene[%d]", i);
		getInputNumber(&inp, key.c_str(), _rfene.data() + i, 0);

		key = Utils::sformat("PS_Kfene[%d]", i);
		getInputNumber(&inp, key.c_str(), _Kfene.data() + i, 0);

		key = Utils::sformat("PS_WCA_sigma[%d]", i);
		getInputNumber(&inp, key.c_str(), _WCA_sigma.data() + i, 0);
	}
}

void PolymerSwapInteraction::init() {
	std::transform(_rfene.begin(), _rfene.end(), _sqr_rfene.begin(), [](number r) {
		return SQR(r);
	});

	std::transform(_WCA_sigma.begin(), _WCA_sigma.end(), _PS_sqr_rep_rcut.begin(), [this](number sigma) {
		return pow(2. * sigma, 2. / this->_PS_n);
	});

	_rcut = sqrt(*std::max_element(_PS_sqr_rep_rcut.begin(), _PS_sqr_rep_rcut.end()));

	_3b_rcut = _3b_range * _3b_sigma;
	_sqr_3b_rcut = SQR(_3b_rcut);
	number B_ss = 1. / (1. + 4. * SQR(1. - _3b_range));
	_3b_A_part = -1. / (B_ss - 1.) / exp(1. / (1. - _3b_range));
	_3b_B_part = B_ss * pow(_3b_sigma, 4.);
	_3b_prefactor = _3b_lambda;
	if(_3b_epsilon > 0.) {
		_3b_prefactor /= _3b_epsilon;
	}

	if(_3b_rcut > _rcut) {
		_rcut = _3b_rcut;
	}

	if(_PS_alpha > 0.) {
		if(_rfene[0] > _rcut) {
			_rcut = _rfene[0];
		}
	}
	_sqr_rcut = SQR(_rcut);

	OX_LOG(Logger::LOG_INFO, "PolymerSwap: A_part: %lf, B_part: %lf, 3b_eps: %lf, total rcut: %lf (%lf)", _3b_A_part, _3b_B_part, _3b_epsilon, _rcut, _sqr_rcut);

	if(_PS_alpha != 0) {
		if(_PS_alpha < 0.) {
			throw oxDNAException("MG_alpha may not be negative");
		}
		_PS_gamma = M_PI / (_sqr_rfene[0] - pow(2., 1. / 3.));
		_PS_beta = 2 * M_PI - _sqr_rfene[0] * _PS_gamma;
		OX_LOG(Logger::LOG_INFO, "MG: alpha = %lf, beta = %lf, gamma = %lf", _PS_alpha, _PS_beta, _PS_gamma);
	}
}

void PolymerSwapInteraction::_update_inter_chain_stress_tensor(int chain, int ref_chain, LR_vector group_force) {
	LR_vector chain_pos = _chain_coms[ref_chain];
	if(chain != ref_chain) {
		chain_pos += _box->min_image(_chain_coms[ref_chain], _chain_coms[chain]);
	}

	_inter_chain_stress_tensor[0] += chain_pos[0] * group_force[0];
	_inter_chain_stress_tensor[1] += chain_pos[1] * group_force[1];
	_inter_chain_stress_tensor[2] += chain_pos[2] * group_force[2];
	_inter_chain_stress_tensor[3] += chain_pos[0] * group_force[1];
	_inter_chain_stress_tensor[4] += chain_pos[0] * group_force[2];
	_inter_chain_stress_tensor[5] += chain_pos[1] * group_force[2];
}

number PolymerSwapInteraction::P_inter_chain() {
	number V = CONFIG_INFO->box->V();
	return _T * (_N_chains / V) + (_inter_chain_stress_tensor[0] + _inter_chain_stress_tensor[1] + _inter_chain_stress_tensor[2]) / (3. * V);
}

void PolymerSwapInteraction::begin_energy_computation() {
	BaseInteraction<PolymerSwapInteraction>::begin_energy_computation();

	std::fill(_inter_chain_stress_tensor.begin(), _inter_chain_stress_tensor.end(), 0.);
	_chain_coms.resize(_N_chains, LR_vector(0., 0., 0.));

	int curr_idx = 0;
	for(int curr_chain = 0; curr_chain < _N_chains; curr_chain++) {
		LR_vector prev_r = CONFIG_INFO->particles()[curr_idx]->pos;
		LR_vector delta_r(0., 0., 0.);
		_chain_coms[curr_chain] = _chain_size * prev_r;
		curr_idx++;
		for(int i = 1; i < _chain_size; i++) {
			BaseParticle *pi = CONFIG_INFO->particles()[curr_idx];
			LR_vector ri = pi->pos;
			delta_r += _box->min_image(prev_r, ri);
			_chain_coms[curr_chain] += delta_r;
			prev_r = ri;

			if(pi->strand_id != curr_chain) {
				throw oxDNAException("mismatch detected between chain ids during the computation of P (%d != %d)", pi->strand_id, curr_chain);
			}

			curr_idx++;
		}
		_chain_coms[curr_chain] /= _chain_size;
	}

	_bonds.clear();
}

bool PolymerSwapInteraction::has_custom_stress_tensor() const {
	return true;
}

number PolymerSwapInteraction::_fene(BaseParticle *p, BaseParticle *q, bool update_forces) {
	number sqr_r = _computed_r.norm();

	int int_type = p->type + q->type;
	number sqr_rfene = _sqr_rfene[int_type];

	if(sqr_r > sqr_rfene) {
		if(update_forces) {
			throw oxDNAException("The distance between particles %d and %d (type: %d, r: %lf) exceeds the FENE distance (%lf)", p->index, q->index, int_type, sqrt(sqr_r), sqrt(sqr_rfene));
		}
		set_is_infinite(true);
		return 10e10;
	}

	number Kfene = _Kfene[int_type];
	number energy = -Kfene * sqr_rfene * log(1. - sqr_r / sqr_rfene);

	if(update_forces) {
		// this number is the module of the force over r, so we don't have to divide the distance vector by its module
		number force_mod = -2 * Kfene * sqr_rfene / (sqr_rfene - sqr_r);
		auto force = -_computed_r * force_mod;
		p->force += force;
		q->force -= force;

		if(p->strand_id != q->strand_id) {
			_update_inter_chain_stress_tensor(p->strand_id, p->strand_id, force);
			_update_inter_chain_stress_tensor(q->strand_id, p->strand_id, -force);
		}

		_update_stress_tensor(p->pos, force);
		_update_stress_tensor(p->pos + _computed_r, -force);
	}

	return energy;
}

number PolymerSwapInteraction::_WCA(BaseParticle *p, BaseParticle *q, bool update_forces) {
	number sqr_r = _computed_r.norm();
	int int_type = p->type + q->type;
	if(sqr_r > _PS_sqr_rep_rcut[int_type]) {
		return (number) 0.;
	}

	number energy = 0;
	// this number is the module of the force over r, so we don't have to divide the distance
	// vector for its module
	number force_mod = 0;

	// cut-off for all the repulsive interactions
	if(sqr_r < _PS_sqr_rep_rcut[int_type]) {
		number part = 1.;
		number ir2_scaled = SQR(_WCA_sigma[int_type]) / sqr_r;
		for(int i = 0; i < _PS_n / 2; i++) {
			part *= ir2_scaled;
		}
		energy += 4. * (part * (part - 1.)) + 1. - _PS_alpha;
		if(update_forces) {
			force_mod += 4. * _PS_n * part * (2. * part - 1.) / sqr_r;
		}
	}
	// attraction
	/*else if(int_type == 0 && _PS_alpha != 0.) {
		energy += 0.5 * _PS_alpha * (cos(_PS_gamma * sqr_r + _PS_beta) - 1.0);
		if(update_forces) {
			force_mod += _PS_alpha * _PS_gamma * sin(_PS_gamma * sqr_r + _PS_beta);
		}
	}*/

	if(update_forces) {
		auto force = -_computed_r * force_mod;
		p->force += force;
		q->force -= force;

		if(p->strand_id != q->strand_id) {
			_update_inter_chain_stress_tensor(p->strand_id, p->strand_id, force);
			_update_inter_chain_stress_tensor(q->strand_id, p->strand_id, -force);
		}

		_update_stress_tensor(p->pos, force);
		_update_stress_tensor(p->pos + _computed_r, -force);
	}

	return energy;
}

number PolymerSwapInteraction::_sticky(BaseParticle *p, BaseParticle *q, bool update_forces) {
	number sqr_r = _computed_r.norm();
	number energy = 0.;

	int int_type = p->type + q->type;
	// sticky-sticky
	if(int_type == 2) {
		if(sqr_r < _sqr_3b_rcut) {
			number r_mod = sqrt(sqr_r);
			number exp_part = exp(_3b_sigma / (r_mod - _3b_rcut));
			number tmp_energy = _3b_epsilon * _3b_A_part * exp_part * (_3b_B_part / SQR(sqr_r) - 1.);

			energy += tmp_energy;

			number tb_energy = (r_mod < _3b_sigma) ? _3b_epsilon : -tmp_energy;
			
			PSBond p_bond(q, tb_energy, _computed_r);
			PSBond q_bond(p, tb_energy, -_computed_r);

			if(update_forces) {
				number force_mod = _3b_epsilon * _3b_A_part * exp_part * (4. * _3b_B_part / (SQR(sqr_r) * r_mod)) + _3b_sigma * tmp_energy / SQR(r_mod - _3b_rcut);
				LR_vector tmp_force = -_computed_r * (force_mod / r_mod);

				p->force += tmp_force;
				q->force -= tmp_force;

				p_bond.force = -tmp_force;
				q_bond.force = +tmp_force;

				if(p->strand_id != q->strand_id) {
					_update_inter_chain_stress_tensor(p->strand_id, p->strand_id, tmp_force);
					_update_inter_chain_stress_tensor(q->strand_id, p->strand_id, -tmp_force);
				}

				_update_stress_tensor(p->pos, tmp_force);
				_update_stress_tensor(p->pos + _computed_r, -tmp_force);
			}

			if(!no_three_body) {
				energy += _three_body(p, p_bond, update_forces);
				energy += _three_body(q, q_bond, update_forces);

				_bonds[p->index].insert(p_bond);
				_bonds[q->index].insert(q_bond);
			}
		}
	}

	return energy;
}

number PolymerSwapInteraction::_three_body(BaseParticle *p, PSBond &new_bond, bool update_forces) {
	number energy = 0.;

	number curr_energy = new_bond.energy;
	for(auto &other_bond : _bonds[p->index]) {
		if(other_bond.other != new_bond.other) {
			number other_energy = other_bond.energy;

			energy += _3b_prefactor * curr_energy * other_energy;

			if(update_forces) {
				if(curr_energy != _3b_epsilon) {
					BaseParticle *other = new_bond.other;

					number factor = -_3b_prefactor * other_energy;
					LR_vector tmp_force = -factor * new_bond.force;

					p->force += tmp_force;
					other->force -= tmp_force;

					if(p->strand_id != other->strand_id) {
						_update_inter_chain_stress_tensor(p->strand_id, p->strand_id, tmp_force);
						_update_inter_chain_stress_tensor(other->strand_id, p->strand_id, -tmp_force);
					}

					_update_stress_tensor(p->pos, tmp_force);
					_update_stress_tensor(p->pos + new_bond.r, -tmp_force);
				}

				if(other_energy != _3b_epsilon) {
					BaseParticle *other = other_bond.other;

					number factor = -_3b_prefactor * curr_energy;
					LR_vector tmp_force = -factor * other_bond.force;

					p->force += tmp_force;
					other->force -= tmp_force;

					if(p->strand_id != other->strand_id) {
						_update_inter_chain_stress_tensor(p->strand_id, p->strand_id, tmp_force);
						_update_inter_chain_stress_tensor(other->strand_id, p->strand_id, -tmp_force);
					}

					_update_stress_tensor(p->pos, tmp_force);
					_update_stress_tensor(p->pos + other_bond.r, -tmp_force);
				}
			}
		}
	}

	return energy;
}

number PolymerSwapInteraction::pair_interaction(BaseParticle *p, BaseParticle *q, bool compute_r, bool update_forces) {
	if(p->is_bonded(q)) {
		return pair_interaction_bonded(p, q, compute_r, update_forces);
	}
	else {
		return pair_interaction_nonbonded(p, q, compute_r, update_forces);
	}
}

number PolymerSwapInteraction::pair_interaction_bonded(BaseParticle *p, BaseParticle *q, bool compute_r, bool update_forces) {
	number energy = (number) 0.f;

	if(p->is_bonded(q)) {
		if(compute_r) {
			if(q != P_VIRTUAL && p != P_VIRTUAL) {
				_computed_r = _box->min_image(p->pos, q->pos);
			}
		}

		energy = _fene(p, q, update_forces);
		energy += _WCA(p, q, update_forces);
	}

	return energy;
}

number PolymerSwapInteraction::pair_interaction_nonbonded(BaseParticle *p, BaseParticle *q, bool compute_r, bool update_forces) {
	if(p->is_bonded(q)) {
		return (number) 0.f;
	}

	if(compute_r) {
		_computed_r = _box->min_image(p->pos, q->pos);
	}

	int int_type = p->type + q->type;
	number energy = _WCA(p, q, update_forces);
	if(int_type == 2) {
		energy += _sticky(p, q, update_forces);
	}

	return energy;
}

void PolymerSwapInteraction::check_input_sanity(std::vector<BaseParticle*> &particles) {
	for(uint i = 0; i < particles.size(); i++) {
		/*CustomParticle *p = static_cast<CustomParticle*>(particles[i]);
		for(auto q : p->bonded_neighs) {

		}*/
	}
}

void PolymerSwapInteraction::allocate_particles(std::vector<BaseParticle*> &particles) {
	for(uint i = 0; i < particles.size(); i++) {
		particles[i] = new CustomParticle();
	}
}

void PolymerSwapInteraction::read_topology(int *N_strands, std::vector<BaseParticle*> &particles) {
	std::string line;
	unsigned int N_from_conf = particles.size();
	IBaseInteraction::read_topology(N_strands, particles);

	std::ifstream topology;
	topology.open(_topology_filename, std::ios::in);
	if(!topology.good()) {
		throw oxDNAException("Can't read topology file '%s'. Aborting", _topology_filename);
	}
	std::getline(topology, line);
	std::vector<int> sticky_particles;
	while(topology.good()) {
		std::getline(topology, line);
		auto ps = Utils::get_particles_from_string(particles, line, "PolymerSwapInteraction");
		sticky_particles.insert(sticky_particles.begin(), ps.begin(), ps.end());
	}
	topology.close();

	OX_LOG(Logger::LOG_INFO, "PolymerSwap: %d sticky particles found in the topology", sticky_particles.size());

	std::ifstream bond_file(_bond_filename.c_str());

	if(!bond_file.good()) {
		throw oxDNAException("Can't read bond file '%s'. Aborting", _bond_filename.c_str());
	}
	// skip the headers and the particle positions if the user specified that the bondfile does not contain bonds only
	if(!_only_links_in_bondfile) {
		int to_skip = N_from_conf + 2;
		for(int i = 0; i < to_skip; i++) {
			std::getline(bond_file, line);
		}
		if(!bond_file.good()) {
			throw oxDNAException("The bond file '%s' does not contain the right number of lines. Aborting", _bond_filename.c_str());
		}
	}

	for(unsigned int i = 0; i < N_from_conf; i++) {
		std::getline(bond_file, line);

		if(!bond_file.good()) {
			throw oxDNAException("The bond file should contain two lines per particle, but it seems there are info for only %d particles\n", i);
		}

		auto line_spl = Utils::split(line, ' ');
		CustomParticle *p;
		int n_bonds;

		unsigned int p_idx;
		if(line_spl.size() == 2) {
			p_idx = std::stoi(line_spl[0]);
			p_idx--;

			n_bonds = std::stoi(line_spl[1]);

			if(i != p_idx) {
				throw oxDNAException("There is something wrong with the bond file. Expected index %d, found %d\n", i, p_idx);
			}
		}
		else {
			n_bonds = std::stoi(line_spl[0]);
			p_idx = i;
		}
		if(p_idx >= N_from_conf) {
			throw oxDNAException("There is a mismatch between the configuration and links files: the latter refers to particle %d, which is larger than the largest possible index (%d)", p_idx, N_from_conf - 1);
		}
		
		p = static_cast<CustomParticle*>(particles[p_idx]);

		p->type = MONOMER;
		if(std::find(sticky_particles.begin(), sticky_particles.end(), p->index) != sticky_particles.end()) {
			p->type = STICKY;
		}
		p->btype = p->type;
		p->strand_id = p_idx / _chain_size;
		p->n3 = p->n5 = P_VIRTUAL;
		for(int j = 0; j < n_bonds; j++) {
			unsigned int n_idx;
			bond_file >> n_idx;
			// the >> operator always leaves '\n', which is subsequently picked up by getline if we don't explicitly ignore it
			bond_file.ignore();
			n_idx--;

			if(n_idx >= N_from_conf) {
				throw oxDNAException("There is a mismatch between the configuration and links files: the latter contains a link between particles %d and %d, which is not possible since the largest possible index in the configuration is %d", p->index, n_idx, N_from_conf - 1);
			}
			
			CustomParticle *q = static_cast<CustomParticle*>(particles[n_idx]);
			p->add_bonded_neigh(q);
		}
	}

	*N_strands = N_from_conf / _chain_size;
	if(*N_strands == 0) {
		*N_strands = 1;
	}
	_N_chains = *N_strands;
}

extern "C" PolymerSwapInteraction* make_PolymerSwapInteraction() {
	return new PolymerSwapInteraction();
}

