/*
 * SkewTrap.cpp
 *
 *  Created on: 20/oct/2021
 *      Author: Jonah
 */

#include "SkewTrap.h"
#include "../Particles/BaseParticle.h"
#include "../Boxes/BaseBox.h"

SkewTrap::SkewTrap() :
				BaseForce() {
	_ref_id = -2;
	_particle = -2;
	_p_ptr = NULL;
	_r0 = -1.;
	PBC = false;
	_box_ptr = NULL;
	_val1 = 0.f;
	_val2 = 0.f;
	_val3 = 0.f;
	_val4 = 0.f;
}

std::tuple<std::vector<int>, std::string> SkewTrap::init(input_file &inp, BaseBox * box_ptr) {
	getInputInt(&inp, "particle", &_particle, 1);
	getInputInt(&inp, "ref_particle", &_ref_id, 1);
	getInputNumber(&inp, "r0", &_r0, 1);
	getInputNumber(&inp, "stdev", &_s, 1);
	getInputNumber(&inp, "shape", &_a, 1);
	getInputBool(&inp, "PBC", &PBC, 0);
	_rate = 0.f; //default rate is 0
	getInputNumber(&inp, "rate", &_rate, 0);

	// set precomputed values
	_val1 = 1.f / (2*pow(_s, 2)); // 1/(2s^2)
	_val2 = sqrt(2*PI) * _s; // Sqrt(2*Pi) * s
	_val3 = -1.f/ pow(_s, 2); // -1/(2s^2)
	_val4 = sqrt(2/PI);

	int N = CONFIG_INFO->particles().size();
	if(_ref_id < 0 || _ref_id >= N) {
		throw oxDNAException("Invalid reference particle %d for Skew Trap", _ref_id);
	}
	_p_ptr = CONFIG_INFO->particles()[_ref_id];

	_box_ptr = box_ptr;

	if(_particle >= N || N < -1) {
		throw oxDNAException("Trying to add a SkewTrap on non-existent particle %d. Aborting", _particle);
	}
	if(_particle == -1) {
		throw oxDNAException("Cannot apply SkewTrap to all particles. Aborting");
	}

	std::string description = Utils::sformat("SkewTrap (stdev=%g, shape=%g, rate=%g, r0=%g, ref_particle=%d, PBC=%d)", _s, _a, _rate, _r0, _ref_id, PBC);

	return std::make_tuple(std::vector<int>{_particle}, description);
}

LR_vector SkewTrap::_distance(LR_vector u, LR_vector v) {
	if(PBC)
		return _box_ptr->min_image(u, v);
	else
		return v - u;
}

LR_vector SkewTrap::value(llint step, LR_vector &pos) {
    // Calculates: - x/ (s^2) + (a e^((-1/2) *a^2 * x^2) * Sqrt(2/Pi)) / (1 + Erf[(a*x)/Sqrt(2)])]
	LR_vector dr = _distance(pos, _box_ptr->get_abs_pos(_p_ptr));
    number dx = dr.module() - (_r0 + (_rate * step));
    number numerator = _a*exp(pow(dx, 2) * pow(_a, 2) * -0.5f) * _val3; // (a e^((-1/2) *a^2 * x^2) * Sqrt(2/Pi))
    number denominator = 1 + erf(_a*dx*0.7071067811865475f); // 1 +Erf[ (a*x)/ Sqrt(2) ]
    return (dr / dr.module()) * (dx * _val4 + numerator/denominator);
}

number SkewTrap::potential(llint step, LR_vector &pos) {
    // Calculates: Log[(e^(x^2/(2s^2))*Sqrt(2Pi)*s)  /  (1 + Erf[(a*x)/Sqrt(2)])]
	LR_vector dr = _distance(pos, _box_ptr->get_abs_pos(_p_ptr));
    number dx = dr.module() - (_r0 + (_rate * step));
    number numerator = exp(pow(dx, 2) * _val1) * _val2;
    number denominator = 1 + erf(_a*dx*0.7071067811865475f); // 1 +Erf[ (a*x)/ Sqrt(2) ]
	return log(numerator/denominator);
}
