/*
 * RNANMInteraction.cpp
 *
 *  Created on: Apr 17, 2019
 *      Author: jonah
 *  Inherits from the RNAInteraction2 Class
 *  Uses RNA2 Model and ANMProtein Model
 *  ANMProtein Functions are implemented here as to avoid a multitude of multiple inheritance
 *  problems of which I am unsure are truly solvable
 */


#include "RNANMInteraction.h"
#include <sstream>
#include <fstream>

#include "../Particles/RNANucleotide.h"
#include "../Particles/ANMParticle.h"


RNANMInteraction::RNANMInteraction() : RNA2Interaction() { // @suppress("Class members should be properly initialized")
    // TODO: Re-examine These
    //Protein Methods Function Pointers
    _int_map[SPRING] = (number (RNAInteraction::*)(BaseParticle *p, BaseParticle *q, bool compute_r, bool update_forces)) &RNANMInteraction::_protein_spring;
    _int_map[PRO_EXC_VOL] = (number (RNAInteraction::*)(BaseParticle *p, BaseParticle *q, bool compute_r, bool update_forces)) &RNANMInteraction::_protein_exc_volume;

    //Protein-RNA Function Pointers
    _int_map[PRO_RNA_EXC_VOL] = (number (RNAInteraction::*)(BaseParticle *p, BaseParticle *q, bool compute_r, bool update_forces)) &RNANMInteraction::_protein_rna_exc_volume;

    //RNA Methods Function Pointers
    _int_map[BACKBONE] = (number (RNAInteraction::*)(BaseParticle *p, BaseParticle *q, bool compute_r, bool update_forces)) &RNANMInteraction::_backbone;
    _int_map[COAXIAL_STACKING] = (number (RNAInteraction::*)(BaseParticle *p, BaseParticle *q, bool compute_r, bool update_forces))  &RNANMInteraction::_coaxial_stacking;
    _int_map[CROSS_STACKING] = (number (RNAInteraction::*)(BaseParticle *p, BaseParticle *q, bool compute_r, bool update_forces))  &RNANMInteraction::_cross_stacking;
    _int_map[BONDED_EXCLUDED_VOLUME] = (number (RNAInteraction::*)(BaseParticle *p, BaseParticle *q, bool compute_r, bool update_forces))  &RNANMInteraction::_bonded_excluded_volume;
    _int_map[STACKING] = (number (RNAInteraction::*)(BaseParticle *p, BaseParticle *q, bool compute_r, bool update_forces))  &RNANMInteraction::_stacking;
    _int_map[HYDROGEN_BONDING] = (number (RNAInteraction::*)(BaseParticle *p, BaseParticle *q, bool compute_r, bool update_forces))  &RNANMInteraction::_hydrogen_bonding;
    _int_map[NONBONDED_EXCLUDED_VOLUME] = (number (RNAInteraction::*)(BaseParticle *p, BaseParticle *q, bool compute_r, bool update_forces))  &RNANMInteraction::_nonbonded_excluded_volume;
    _int_map[DEBYE_HUCKEL] = (number (RNAInteraction::*)(BaseParticle *p, BaseParticle *q, bool compute_r, bool update_forces))  &RNANMInteraction::_debye_huckel;
}


void RNANMInteraction::get_settings(input_file &inp){
    this->RNA2Interaction::get_settings(inp);

    getInputString(&inp, "parfile", _parameterfile, 0);
    //Addition of Reading Parameter File
    char n[5] = "none";

    auto valid_spring_params = [](int N, int x, int y, double d, char s, double k){
        if(x < 0 || x > N) throw oxDNAException("Invalid Particle ID %d in Parameter File", x);
        if(y < 0 || y > N) throw oxDNAException("Invalid Particle ID %d in Parameter File", y);
        if(d < 0) throw oxDNAException("Invalid Eq Distance %d in Parameter File", d);
        if(s != 's') throw oxDNAException("Potential Type %c Not Supported", s);
        if(k < 0) throw oxDNAException("Spring Constant %f Not Supported", k);
    };

    if(strcmp(_parameterfile, n) != 0) {
        int key1, key2;
        char potswitch;
        double potential, dist;
        int N;
        std::fstream parameters;
        parameters.open(_parameterfile, std::ios::in);
        parameters >> N;
        int spring_connection_num = 0;
        if (parameters.is_open())
        {
            while (parameters >> key1 >> key2 >> dist >> potswitch >> potential)
            {
                valid_spring_params(N, key1, key2, dist, potswitch, potential);
                spring_connection_num += 1;
                std::pair <int, int> lkeys (key1, key2);
                std::pair <char, double> pot (potswitch, potential);
                _rknot[lkeys] = dist;
                _potential[lkeys] = pot;
            }
        }
        else
        {
            throw oxDNAException("ParameterFile Could Not Be Opened");
        }
        parameters.close();
        if(spring_connection_num == 1 && N > 2) throw oxDNAException("Invalid Parameter File Format, cannot use a RNACT Parameter File");
    } else {
        OX_LOG(Logger::LOG_INFO, "Parfile: NONE, No protein parameters were filled");
    }
}


void RNANMInteraction::check_input_sanity(std::vector<BaseParticle*> &particles){
    //this->RNA2Interaction::check_input_sanity(particles,N);
    //Need to make own function that checks the input sanity
}



void RNANMInteraction::allocate_particles(std::vector<BaseParticle*> &particles) {
    if (nrna==0 || nrnas==0) {
        OX_LOG(Logger::LOG_INFO, "No RNA Particles Specified, Continuing with just Protein Particles");
        for (int i = 0; i < npro; i++) particles[i] = new ANMParticle();
    } else if (npro == 0) {
        OX_LOG(Logger::LOG_INFO, "No Protein Particles Specified, Continuing with just RNA Particles");
        for (int i = 0; i < nrna; i++) particles[i] = new RNANucleotide(this->_grooving);
    } else {
        if (_firststrand > 0){
            for (int i = 0; i < nrna; i++) particles[i] = new RNANucleotide(this->_grooving);
            for (uint i = nrna; i < particles.size(); i++) particles[i] = new ANMParticle();
        } else {
            for (int i = 0; i < npro; i++) particles[i] = new ANMParticle();
            for (uint i = npro; i < particles.size(); i++) particles[i] = new RNANucleotide(this->_grooving);
        }
    }
}


void RNANMInteraction::read_topology(int *N_strands, std::vector<BaseParticle*> &particles) {
    int my_N, my_N_strands;

    char line[5120];
    std::ifstream topology;
    topology.open(this->_topology_filename, std::ios::in);

    if (!topology.good()) throw oxDNAException("Can't read topology file '%s'. Aborting",this->_topology_filename);

    topology.getline(line, 5120);
    std::stringstream head(line);

    head >> my_N >> my_N_strands >> nrna >> npro >>nrnas;
    if (head.fail()) throw oxDNAException("Problem with header make sure the format is correct for RNANM Interaction");

    if(my_N_strands < 0 || my_N_strands > my_N || nrna > my_N || nrna < 0 || npro > my_N || npro < 0 || nrnas < 0 || nrnas > my_N) {
        throw oxDNAException("Problem with header make sure the format is correct for RNANM Interaction");
    }


    int strand, i = 0;
    while (topology.good()) {
        topology.getline(line, 5120);
        if (strlen(line) == 0 || line[0] == '#')
            continue;
        if (i == my_N)
            throw oxDNAException("Too many particles found in the topology file (should be %d), aborting", my_N);

        std::stringstream ss(line);
        ss >> strand;
        if (i == 0) {
            _firststrand = strand; //Must be set prior to allocation of particles
            allocate_particles(particles);
            for (int j = 0; j < my_N; j++) {
                particles[j]->index = j;
                particles[j]->type = 26; //A_INVALID
                particles[j]->strand_id = 0;
            }
        }

        // Amino Acid
        if (strand < 0) {
            char aminoacid[256];
            int nside, cside;
            ss >> aminoacid >> nside >> cside;

            int x;
            std::set<int> myneighs;
            if (nside >= 0) myneighs.insert(nside);
            if (cside >= 0) myneighs.insert(cside);
            while (ss.good()) {
                ss >> x;
                if (x < 0 || x >= my_N) {
                    throw oxDNAException("Line %d of the topology file has an invalid syntax, neighbor has invalid id", i + 2);
                }
                myneighs.insert(x);
            }

            auto *p = dynamic_cast< ANMParticle * > (particles[i]);

            if (strlen(aminoacid) == 1) {
                p->type = Utils::decode_aa(aminoacid[0]);
                p->btype = Utils::decode_aa(aminoacid[0]);
            }

            p->strand_id = abs(strand) + nrnas - 1;
            p->index = i;

            for(auto & k : myneighs){
                if (p->index < k) {
                    p->add_bonded_neighbor(dynamic_cast<ANMParticle *> (particles[k]));
                }
            }

            i++;
        }
        if (strand > 0) {
            char base[256];
            int tmpn3, tmpn5;
            ss >> base >> tmpn3 >> tmpn5;

            auto *p = dynamic_cast<RNANucleotide *> (particles[i]);

            if (tmpn3 < 0) p->n3 = P_VIRTUAL;
            else p->n3 = particles[tmpn3];
            if (tmpn5 < 0) p->n5 = P_VIRTUAL;
            else p->n5 = particles[tmpn5];


            p->strand_id = strand - 1;

            // the base can be either a char or an integer
            if (strlen(base) == 1) {
                p->type = Utils::decode_base(base[0]);
                p->btype = Utils::decode_base(base[0]);

            } else {
                if (atoi(base) > 0) p->type = atoi(base) % 4;
                else p->type = 3 - ((3 - atoi(base)) % 4);
                p->btype = atoi(base);
            }

            if (p->type == P_INVALID) throw oxDNAException("Particle #%d in strand #%d contains a non valid base '%c'. Aborting", i, strand, base);

            p->index = i;
            i++;

            // here we fill the affected vector
            if (p->n3 != P_VIRTUAL) p->affected.emplace_back(ParticlePair(p->n3, p));
            if (p->n5 != P_VIRTUAL) p->affected.emplace_back(ParticlePair(p, p->n5));

        }
        if (strand == 0) throw oxDNAException("No strand 0 should be present please check topology file");

    }


    if (i < my_N)
        throw oxDNAException("Not enough particles found in the topology file (should be %d). Aborting", my_N);

    topology.close();

    if (my_N != (int) particles.size())
        throw oxDNAException("Number of lines in the configuration file and number of particles in the topology files don't match. Aborting");

    *N_strands = my_N_strands;

}



number RNANMInteraction::pair_interaction(BaseParticle *p, BaseParticle *q, bool compute_r, bool update_forces){
    if (p->btype < 4 && q->btype < 4){
        if(p->is_bonded(q)) return pair_interaction_bonded(p, q, compute_r, update_forces);
        else return pair_interaction_nonbonded(p, q, compute_r, update_forces);
    }
    if ((p->btype < 4 && q->btype > 4) || (p->btype > 4 && q->btype < 4)) return this->pair_interaction_nonbonded(p, q, compute_r, update_forces);

    if (p->btype > 4 && q->btype > 4){
        auto *cp = dynamic_cast< ANMParticle * > (p);
        if ((*cp).is_bonded(q)) return pair_interaction_bonded(p, q, compute_r, update_forces);
        else return pair_interaction_nonbonded(p, q, compute_r, update_forces);
    }
    return 0.f;
}


number RNANMInteraction::pair_interaction_bonded(BaseParticle *p, BaseParticle *q, bool compute_r, bool update_forces) {
    if(compute_r)
        if (q != P_VIRTUAL && p != P_VIRTUAL)
            _computed_r = this->_box->min_image(p->pos, q->pos);

    if (p->btype < 4 && q->btype < 4){
        if(!this->_check_bonded_neighbour(&p, &q, compute_r)) return (number) 0;
        number energy = _backbone(p, q, false, update_forces);
        energy += _bonded_excluded_volume(p, q, false, update_forces);
        energy += _stacking(p, q, false, update_forces);
        return energy;
    } else if (p->btype > 4 && q->btype > 4){
        auto *cp = dynamic_cast< ANMParticle * > (p);
        if (!(*cp).is_bonded(q)) return 0.f;
        number energy = _protein_spring(p, q, compute_r, update_forces);
        energy += _protein_exc_volume(p, q, compute_r, update_forces);
        return energy;
    } else
        return 0.f;

    // Expect Topology to not contain Bonds b/t protein & RNA
    //if ((p->btype >= 0 && q->btype <0) ||(p->btype <0 && q->btype >=0)) return 0.f;
}


number RNANMInteraction::pair_interaction_nonbonded(BaseParticle *p, BaseParticle *q, bool compute_r, bool update_forces) {
    if (compute_r)
        _computed_r = this->_box->min_image(p->pos, q->pos);

    number rnorm = _computed_r.norm();
    if (p->btype < 4 && q->btype < 4) { //RNA-RNA Interaction
        if (rnorm >= _sqr_rcut) return (number) 0.f;
        number energy = RNAInteraction::pair_interaction_nonbonded(p, q, false, update_forces);
        energy += _debye_huckel(p, q, false, update_forces);
        return energy;
    }

    if ((p->btype < 4 && q->btype > 4) || (p->btype > 4 && q->btype < 4)) {
        if (rnorm >= _pro_rna_sqr_rcut) return (number) 0.f;
        number energy = _protein_rna_exc_volume(p, q, compute_r, update_forces);
        return energy;
    }

    if (p->btype > 4 && q->btype > 4) {
        if (rnorm >= _pro_sqr_rcut) return (number) 0.f;
        number energy = _protein_exc_volume(p, q, compute_r, update_forces);
        return energy;
    }

    return 0.f;
}


number RNANMInteraction::_protein_rna_exc_volume(BaseParticle *p, BaseParticle *q, bool compute_r, bool update_forces)
{
    BaseParticle *protein;
    BaseParticle *nuc;

    LR_vector force(0, 0, 0);
    LR_vector rcenter = _computed_r;

    if(p->btype < 4 && q->btype > 4)
    {
        //rcenter = -rcenter;
        protein = q;
        nuc = p;
    }
    else if (p->btype > 4 && q->btype < 4)
    {
        rcenter = -rcenter;
        protein = p;
        nuc = q;
    }
    else
        return 0.f;

    LR_vector r_to_back = rcenter  - nuc->int_centers[RNANucleotide::BACK];
    LR_vector r_to_base = rcenter  - nuc->int_centers[RNANucleotide::BASE];

    LR_vector torquenuc(0,0,0);
    number energy;

    if(r_to_back.module() < _pro_backbone_sqr_rcut) {
        energy = _protein_rna_repulsive_lj(r_to_back, force, update_forces, _pro_backbone_sigma, _pro_backbone_b, _pro_backbone_rstar,_pro_backbone_rcut,_pro_backbone_stiffness);
        //printf("back-pro %d %d %f\n",p->index,q->index,energy);
        if (update_forces) {
            torquenuc  -= nuc->int_centers[RNANucleotide::BACK].cross(force);
            nuc->force -= force;
            protein->force += force;
        }
    }

    if(r_to_base.module() > _pro_base_sqr_rcut) return energy;

    energy += _protein_rna_repulsive_lj(r_to_base, force, update_forces, _pro_base_sigma, _pro_base_b, _pro_base_rstar, _pro_base_rcut, _pro_base_stiffness);
    if(update_forces) {
        torquenuc  -= nuc->int_centers[RNANucleotide::BASE].cross(force);
        nuc->torque += nuc->orientationT * torquenuc;

        nuc->force -= force;
        protein->force += force;
    }


    return energy;
}


number RNANMInteraction::_protein_rna_repulsive_lj(const LR_vector &r, LR_vector &force, bool update_forces, number &sigma, number &b, number &rstar, number &rcut, number &stiffness) {
    // this is a bit faster than calling r.norm()
    number rnorm = SQR(r.x) + SQR(r.y) + SQR(r.z);
    number energy = (number) 0;

    if(rnorm < SQR(rcut)) {
        if(rnorm > SQR(rstar)) {
            number rmod = sqrt(rnorm);
            number rrc = rmod - rcut;
            energy = stiffness * b * SQR(SQR(rrc));
            if(update_forces) force = -r * (4 * stiffness * b * CUB(rrc) / rmod);
        }
        else {
            number tmp = SQR(sigma) / rnorm;
            number lj_part = tmp * tmp * tmp;
            energy = 4 * stiffness * (SQR(lj_part) - lj_part);
            if(update_forces) force = -r * (24 * stiffness * (lj_part - 2*SQR(lj_part)) / rnorm);
        }
    }

    if(update_forces && energy == (number) 0) force.x = force.y = force.z = (number) 0;

    return energy;
}


void RNANMInteraction::init() {
    this->RNA2Interaction::init();
    nrna=0, npro=0, nrnas =0;
    //let's try this
    _pro_backbone_sigma = 0.57f;
    _pro_backbone_rstar= 0.569f;
    _pro_backbone_b = 178699253.5f;
    _pro_backbone_rcut = 0.572934f;
    _pro_backbone_stiffness = 1.0f;
    _pro_backbone_sqr_rcut = 0.3283f;
    //Base-Protein Excluded Volume Parameters
    _pro_base_sigma = 0.36f;
    _pro_base_rstar= 0.359f;
    _pro_base_b = 296866090.f;
    _pro_base_rcut = 0.362897f;
    _pro_base_stiffness = 1.0f;
    _pro_base_sqr_rcut = 0.1317f;
    //Protein-Protein Excluded Volume Parameters
    _pro_sigma = 0.35f;
    _pro_rstar = 0.349f;
    _pro_b = 306484596.421f;
    _pro_rcut = 0.352894;
    _pro_sqr_rcut = 0.12454f; //_rc ^2

    _pro_rna_sqr_rcut = 3.0625f; //placeholder test value 1.75^2
}

//Functions from ANMInteraction.h
//Stolen due to inheritance issues


number RNANMInteraction::_protein_repulsive_lj(const LR_vector &r, LR_vector &force, bool update_forces) {
    // this is a bit faster than calling r.norm()
    //changed to a quartic form
    number rnorm = SQR(r.x) + SQR(r.y) + SQR(r.z);
    number energy = (number) 0;
    if(rnorm < SQR(_pro_rcut)) {
        if(rnorm > SQR(_pro_rstar)) {
            number rmod = sqrt(rnorm);
            number rrc = rmod - _pro_rcut;
            energy = EXCL_EPS * _pro_b * SQR(SQR(rrc));
            if(update_forces) force = -r * (4 * EXCL_EPS * _pro_b * CUB(rrc)/ rmod);
        }
        else {
            number tmp = SQR(_pro_sigma) / rnorm;
            number lj_part = tmp * tmp * tmp;
            energy = 4 * EXCL_EPS * (SQR(lj_part) - lj_part);
            if(update_forces) force = -r* (24 * EXCL_EPS * (lj_part - 2*SQR(lj_part))/rnorm);
        }
    }

    if(update_forces && energy == (number) 0) force.x = force.y = force.z = (number) 0;

    return energy;
}


number RNANMInteraction::_protein_exc_volume(BaseParticle *p, BaseParticle *q, bool compute_r, bool update_forces) {
    LR_vector force(0,0,0);
    number energy =  _protein_repulsive_lj(_computed_r, force, update_forces);

    if(update_forces)
    {
        p->force -= force;
        q->force += force;
    }

    return energy;
}


number RNANMInteraction::_protein_spring(BaseParticle *p, BaseParticle *q, bool compute_r, bool update_forces) {

    std::pair <int,int> keys (std::min(p->index, q->index), std::max(p->index, q->index));

    //Harmonic Spring Potential
    number _k = _potential[keys].second; //stiffness of the spring

    number rinsta = _computed_r.module(); // distance b/t p and q
    number disp = rinsta - _rknot[keys]; // current distance - eqdistance
    number energy = 0.5 * _k * SQR(disp);

    if (update_forces) {
        LR_vector force(_computed_r);
        force *= (-1.0f * _k ) * disp/rinsta;

        p->force -= force;
        q->force += force;
    }

    return energy;
}

RNANMInteraction::~RNANMInteraction() {
}





