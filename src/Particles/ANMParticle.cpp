//
// Created by jonah on 2/25/21.
//

#include "ANMParticle.h"


ANMParticle::ANMParticle() : BaseParticle()  {

}


ANMParticle::~ANMParticle() {

}


void ANMParticle::add_bonded_neighbor(ANMParticle *nn) {
    if(!is_bonded(nn)) {
        bonded_neighs.insert(nn);
        nn->bonded_neighs.insert(this);

        ParticlePair new_pair(this, nn);
        this->affected.push_back(new_pair);
        nn->affected.push_back(new_pair);
    }
}



bool ANMParticle::is_bonded(BaseParticle *q) {
    auto *Cq = dynamic_cast<ANMParticle *>(q);
    return !(bonded_neighs.find(Cq) == bonded_neighs.end());
}
