/**
 * @brief Hybrid DNA/ANM Model
 *
 *
 *
 * Jonah Feb. 2021
 * 
 * To use the oxDNA2 model with ANM Protein, set
 *
 * interaction_type = CGDNA
 *
 * in the input file
 *
 * Input options:
 *
 * @verbatim
 parfile = string (set parameter file for protein component, set to none for DNA only sims)
 massfile = string (set massfile for simulations w/ different massed items)
 */

#ifndef CGDNA_INTERACTION_H
#define CGDNA_INTERACTION_H

#include "DNA2Interaction.h"


class CGDNAInteraction: public DNA2Interaction {

protected:
    int ndna;//How many particles of DNA type: Used in allocate_particles
    int npro;//How many bonds b/t different particle types
    int ndnas;//Number of Strands that are dna, rest assumed to be protein
    int _firststrand; //+ for dna, - for protein or gs
    int ngs; // number of generic sphere particles
    int npep; // number of protein strands
    int ngstrands; // number of gs strands

    // excluded volume parameters for quartic LJ on backbone/protein + base/protein
    number _pro_backbone_sigma, _pro_backbone_rstar, _pro_backbone_b, _pro_backbone_rcut, _pro_backbone_stiffness;
    number _pro_base_sigma,_pro_base_rstar, _pro_base_b, _pro_base_rcut, _pro_base_stiffness;
    std::map<std::pair<int, int>, double> _rknot; //Both maps used just as they are in ACInteraction
    std::map<std::pair<int, int>, std::pair<char, double> > _potential;
    std::map<int, number> masses;
    bool _parameter_kbkt; //Controls whether kb/kt values are global or read from parameter file
    std::map<int, std::pair<double, double> > _ang_stiff; // Stores per particle pair, kb kt values
    number _k_bend, _k_tor;
    std::map<int, std::vector <double> > _ang_vals;
    number _pro_sigma, _pro_rstar, _pro_b, _pro_rcut; // Protein-protein quartic LJ params same as in ANM
    number _pro_base_sqr_rcut, _pro_backbone_sqr_rcut, _pro_sqr_rcut, _pro_dna_sqr_rcut;

    bool _fill_in_constants(number interaction_radius, number& sigma, number& rstar, number &b, number &rc)
    {
        rstar = interaction_radius;
        sigma = rstar + 0.05f;
        b = _calc_b(sigma,rstar);
        rc = _calc_rc(sigma,rstar);
        return _check_repulsion_smoothness(sigma,rstar,b,rc);

    }

    number _calc_b(number sigma,number R)
    {
        return 36.0f * POW6(sigma)  * SQR((- POW6(R) + 2.0f * POW6(sigma)))   / (POW14(R) *(-POW6(R)+ POW6(sigma)));
    }

    number _calc_rc(number sigma,number R)
    {
        return R *(4.0f *  POW6(R) -7.0f * POW6(sigma)  ) /(3.0f * (POW6(R) -2.0f *  POW6(sigma) ));
    }

    bool _check_repulsion_smoothness(number sigma,number rstar, number b, number rc)
    {
        /*TODO IMPLEMENT ME NICER!*/
        return true;
        // printf
        // printf("##### STARTING checking for %f\n",rstar);
        double tolerance = 0.2;
        LR_vector<number> rr(0,0,rstar-0.2);
        LR_vector<number> forcer(0,0,0);
        number old_energy =  this->_crowder_repulsive_lj(rr,forcer, 1.,sigma, rstar,  b, rc, true);
        number oldforce = forcer.norm();
        for (double x = rstar-0.2; x < rc+0.01; x += 0.0001)
        {
            LR_vector<number> r(0,0,x);
            LR_vector<number> force(0,0,0);
            number energy_new = this->_crowder_repulsive_lj(r,force, 1.,sigma, rstar,  b, rc, true) ;
            number force_new = force.norm();
            //printf("#### %f %f %f \n",x,energy_new,force_new);
            //if( fabs(old_energy - energy_new) > tolerance || fabs(force_new - oldforce) > tolerance || old_energy < energy_new)
            if(  old_energy < energy_new || oldforce < force_new)
            {
                throw oxDNAException("Non-monotonous change in the repulsion potential for distance  r = %f",x);
                return false;
            }
            else
            {
                old_energy = energy_new;
                oldforce = force_new;
            }
        }
        return true;
    };


public:
    enum {
        SPRING = 8,
        PRO_EXC_VOL = 9,
        PRO_DNA_EXC_VOL = 10,
        PRO_ANG_POT = 11,
        GS_DNA_EXC_VOL = 12,
        GS_PRO_EXC_VOL = 13
        //Assigned 8 9 and 10 so it won't overwrite the already existing DNA function pointers in the _int_map
    };

    char _parameterfile[512];
    std::string _massfile;
    bool _angular;

    CGDNAInteraction(bool btp); //btn controls whether bending/torsion potential is applied
    void load_massfile(std::string &filename);
    virtual ~CGDNAInteraction();
    virtual void get_settings(input_file &inp);
    virtual void allocate_particles(std::vector<BaseParticle *> &particles);
    virtual void read_topology(int *N_strands, std::vector<BaseParticle *> &particles);
    void load_masses(std::string &massfile);
    virtual number pair_interaction(BaseParticle *p, BaseParticle *q, bool compute_r, bool update_forces);
    virtual number pair_interaction_bonded(BaseParticle *p, BaseParticle *q, bool compute_r, bool update_forces);
    virtual number pair_interaction_nonbonded(BaseParticle *p, BaseParticle *q, bool compute_r, bool update_forces);
    number _protein_dna_exc_volume(BaseParticle *p,BaseParticle *q, bool compute_r, bool update_forces);
    number _protein_dna_repulsive_lj(const LR_vector &r, LR_vector &force, bool update_forces, number &sigma, number &b, number &rstar, number &rcut, number &stiffness);
    virtual void check_input_sanity(std::vector<BaseParticle *> &particles);
    virtual void init();
    virtual number _protein_repulsive_lj(const LR_vector &r, LR_vector &force, bool update_forces);
    virtual number _protein_exc_volume(BaseParticle *p, BaseParticle *q, bool compute_r, bool update_forces);
    virtual number _protein_spring(BaseParticle *p, BaseParticle *q, bool compute_r, bool update_forces);
    virtual number _protein_ang_pot(BaseParticle *p, BaseParticle*q, bool compute_r, bool update_forces);
};

#endif /* CGDNA_INTERACTION_H */
