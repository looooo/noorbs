#include "unwrap.h"
#include "nurbs.h"


namespace nurbs{

NurbsFlat::NurbsFlat(RowMat<double, 3> poles, 
					 NurbsBase2D nurbs_base, 
        			 std::vector<long> fixed_pins={})
{
	this->poles = poles;
	this->fixed_pins = fixed_pins;
	this->base = nurbs_base;
	// we have to make a map from u_i, v_i to matrix entry and the other direction
	// find the straight direction!

    this->set_fixed_pins();

    
    int fixed_count = 0;
    for (long i=0; i < this->vertices.cols(); i++)
    {   
        if (fixed_count < this->fixed_pins.size())
        {
            if (i == this->fixed_pins[fixed_count])
                fixed_count ++;
            else
                this->old_order.push_back(i);
        }
        else
            this->old_order.push_back(i);
    }

    for (auto fixed_index: this->fixed_pins)
        this->old_order.push_back(fixed_index);

    // get the reversed map:
    this->new_order.resize(this->old_order.size());
    long j = 0;
    for (auto index: this->old_order)
    {
        this->new_order[index] = j;
        j++;
    }
}

void NurbsFlat::lscm()
{



}

}