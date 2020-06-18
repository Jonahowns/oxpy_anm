/*
 * Pressure.cpp
 *
 *  Created on: 25/ott/2013
 *      Author: lorenzo
 */

#include "Pressure.h"

Pressure::Pressure() :
				BaseObservable(),
				_custom_stress_tensor(true),
				_with_stress_tensor(false),
				_PV_only(false),
				_P(0.),
				_shear_rate(0.) {

}

Pressure::~Pressure() {

}

void Pressure::get_settings(input_file &my_inp, input_file &sim_inp) {
	BaseObservable::get_settings(my_inp, sim_inp);

	char raw_T[256];
	getInputString(&sim_inp, "T", raw_T, 1);
	_T = Utils::get_temperature(raw_T);

	getInputBool(&my_inp, "custom_stress_tensor", &_custom_stress_tensor, 0);
	getInputBool(&my_inp, "stress_tensor", &_with_stress_tensor, 0);
	getInputBool(&my_inp, "PV_only", &_PV_only, 0);

	bool lees_edwards = false;
	getInputBool(&sim_inp, "lees_edwards", &lees_edwards, 0);
	if(lees_edwards) getInputNumber(&sim_inp, "lees_edwards_shear_rate", &_shear_rate, 0);
}

void Pressure::update_pressure() {
	int N = _config_info->N();
	std::vector<ParticlePair> pairs = _config_info->lists->get_potential_interactions();

	_config_info->interaction->begin_energy_computation();

	double virial = 0;
	_stress_tensor = { 0., 0., 0., 0., 0., 0. };
	double energy = 0.;
	// we loop over all the pairs in order to update the forces
	for(auto &pair : pairs) {
		BaseParticle *p = pair.first;
		BaseParticle *q = pair.second;
		LR_vector r = _config_info->box->min_image(p->pos, q->pos);

		// pair_interaction will change these vectors, but we still need them in the next
		// first integration step. For this reason we copy and then restore their values
		// after the calculation
		LR_vector old_p_force(p->force);
		LR_vector old_q_force(q->force);
		LR_vector old_p_torque(p->torque);
		LR_vector old_q_torque(q->torque);

		p->force = q->force = p->torque = q->torque = LR_vector();

		_config_info->interaction->set_computed_r(r);
		energy += (double) _config_info->interaction->pair_interaction(p, q, false, true);

		_stress_tensor[0] -= r.x * p->force.x;
		_stress_tensor[1] -= r.y * p->force.y;
		_stress_tensor[2] -= r.z * p->force.z;
		_stress_tensor[3] -= r.x * p->force.y;
		_stress_tensor[4] -= r.x * p->force.z;
		_stress_tensor[5] -= r.y * p->force.z;

		virial -= (r * p->force);

		p->force = old_p_force;
		q->force = old_q_force;
		p->torque = old_p_torque;
		q->torque = old_q_torque;
	}

	/*
	for(auto p : _config_info->particles()) {
		LR_vector vel = p->vel;
		if(_shear_rate > 0.) {
			number Ly = CONFIG_INFO->box->box_sides().y;
			number y_in_box = p->pos.y - floor(p->pos.y / Ly) * Ly - 0.5 * Ly;
			number flow_vx = y_in_box * _shear_rate;
			vel.x -= flow_vx;
		}
		_stress_tensor.v1.x += SQR(vel.x);
		_stress_tensor.v1.y += vel.x * vel.y;
		_stress_tensor.v1.z += vel.x * vel.z;
		_stress_tensor.v2.x += vel.y * vel.x;
		_stress_tensor.v2.y += SQR(vel.y);
		_stress_tensor.v2.z += vel.y * vel.z;
		_stress_tensor.v3.x += vel.z * vel.x;
		_stress_tensor.v3.y += vel.z * vel.y;
		_stress_tensor.v3.z += SQR(vel.z);
	}
	*/

	double V = (_PV_only) ? 1 : _config_info->box->V();
	_P = _T * (N / V) + virial / (3. * V);
	for(auto &v : _stress_tensor) {
		v /= 3. * V;
	}
}

void Pressure::update_pressure_with_custom_stress_tensor() {
	int N = _config_info->N();
	std::vector<ParticlePair> pairs = _config_info->lists->get_potential_interactions();

	_config_info->interaction->begin_energy_computation();

	for(auto p : _config_info->particles()) {
		p->force = p->torque = LR_vector();
	}

	for(auto &pair : pairs) {
		BaseParticle *p = pair.first;
		BaseParticle *q = pair.second;

		_config_info->interaction->pair_interaction(p, q, true, true);
	}

	double V = (_PV_only) ? 1 : _config_info->box->V();
	_stress_tensor = _config_info->interaction->stress_tensor();
	for(auto &v : _stress_tensor) {
		v /= 3. * V;
	}
	_P = _T * (N / V) + (_stress_tensor[0] + _stress_tensor[1] + _stress_tensor[2]);
}

std::string Pressure::get_output_string(llint curr_step) {
	if(_custom_stress_tensor && _config_info->interaction->has_custom_stress_tensor()) {
		update_pressure_with_custom_stress_tensor();
	}
	else {
		update_pressure();
	}

	std::string to_ret;

	if(_with_stress_tensor) {
		to_ret += Utils::sformat("% .8e % .8e % .8e % .8e % .8e % .8e % .8e", _P, _stress_tensor[0], _stress_tensor[1], _stress_tensor[2], _stress_tensor[3], _stress_tensor[4], _stress_tensor[5]);
	}
	else {
		to_ret += Utils::sformat("% .8e", _P);
	}

	return to_ret;
}
