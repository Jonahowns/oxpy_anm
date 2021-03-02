//
// Created by jonah on 2/25/21.
//

#include "ANMTParticle.h"


LR_vector const ANMTParticle::principal_axis(1, 0, 0);
LR_vector const ANMTParticle::second_axis(0, 0, 1);
LR_vector const ANMTParticle::third_axis(0, 1, 0);

ANMTParticle::ANMTParticle() : BaseParticle()  {

}


ANMTParticle::~ANMTParticle() {

}


void ANMTParticle::add_bonded_neighbor(ANMTParticle *nn) {
    if(!is_bonded(nn)) {
        bonded_neighs.insert(nn);
        nn->bonded_neighs.insert(this);

        ParticlePair new_pair(this, nn);
        this->affected.push_back(new_pair);
        nn->affected.push_back(new_pair);
    }
}



bool ANMTParticle::is_bonded(BaseParticle *q) {
    auto *Cq = dynamic_cast<ANMTParticle *>(q);
    return !(bonded_neighs.find(Cq) == bonded_neighs.end());
}
