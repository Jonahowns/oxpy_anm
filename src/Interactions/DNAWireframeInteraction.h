/**
 * @brief Hybrid DNA/ANM Model
 *
 *
 *
 * Jonah Feb. 2021
 * 
 * To use the oxDNA2 model with ANM Protein, set
 *
 * interaction_type = DNANM
 *f
 * in the input file
 *
 * Input options:
 *
 * @verbatim
 parfile = string (set parameter file for protein component, set to none for DNA only sims)
 massfile = string (set massfile for simulations w/ different massed items)
 */

#ifndef DNAWireframe_INTERACTION_H
#define DNAWireframe_INTERACTION_H

#include "DNA2Interaction.h"


class DNAWireframeInteraction: public DNA2Interaction {

protected:
    int ndna;//How many particles of DNA type: Used in allocate_particles
    int npro;//How many bonds b/t different particle types
    int ndnas;//Number of Strands that are dna, rest assumed to be protein
    int _firststrand; //+ for dna, - for protein

    // excluded volume parameters for quartic LJ on backbone/protein + base/protein
    number _pro_backbone_sigma, _pro_backbone_rstar, _pro_backbone_b, _pro_backbone_rcut, _pro_backbone_stiffness;
    number _pro_base_sigma,_pro_base_rstar, _pro_base_b, _pro_base_rcut, _pro_base_stiffness;
    std::map<std::pair<int, int>, double> _rknot; //Both maps used just as they are in ACInteraction
    std::map<std::pair<int, int>, std::pair<char, double> > _potential;
    number _pro_sigma, _pro_rstar, _pro_b, _pro_rcut; // Protein-protein quartic LJ params same as in ANM
    number _pro_base_sqr_rcut, _pro_backbone_sqr_rcut, _pro_sqr_rcut, _pro_dna_sqr_rcut;

public:
    enum {
        SPRING = 8,
        PRO_EXC_VOL = 9,
        PRO_DNA_EXC_VOL = 10
        //Assigned 8 9 and 10 so it won't overwrite the already existing DNA function pointers in the _int_map
    };

    char _parameterfile[500];

    DNAWireframeInteraction();
    virtual ~DNAWireframeInteraction();
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
    number _protein_repulsive_lj(const LR_vector &r, LR_vector &force, bool update_forces);
    number _protein_exc_volume(BaseParticle *p, BaseParticle *q, bool compute_r, bool update_forces);
    virtual number _protein_spring(BaseParticle *p, BaseParticle *q, bool compute_r, bool update_forces);
};

#endif /* DNAWireframe_INTERACTION_H */
