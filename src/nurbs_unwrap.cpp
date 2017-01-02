#include "unwrap.h"
#include "nurbs.h"


namespace nurbs{

    NurbsFlat::NurbsFlat(RowMat<double, 3> poles, 
              NurbsBase2D* nurbs_base, 
              std::vector<long> fixed_pins)
{
    // for a first implementation fix the last two poles

    // 1: define the collocation points in u-v

    // 2: for every collocation point:
            // 1: get the derivate vector of all controllpoints
            // 2: construct the local coordinatensystem H
            // 3: compute the jacobimatrix
            // 4: fill matrix qith equation (use 3 + 1)

    // 3: use last two collumns to compute rhs knowen values are: (0, 0, 1, 0)

    // 4: solve system and set flat vertices


}

void NurbsFlat::lscm()
{



}

}