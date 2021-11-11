//
//

#ifndef SLICEDOPTIM_IOPOINTSET_H
#define SLICEDOPTIM_IOPOINTSET_H

#include <iostream>
#include <fstream>
#include <vector>
#include "../Math/VecX.h"

template <class VECTYPE>
inline void savePointsetND(std::ostream& out, const std::vector<VECTYPE>& points){
    VECTYPE min,max;
    for(int i=0;i<points.size();++i){
        for (int d = 0; d < min.dim();++d)
        {
            if(points[i][d]<min[d]){
                min[d] = points[i][d];
            }
            if(points[i][d]>max[d]){
                max[d] = points[i][d];
            }
        }
    }
    
    out.precision(11);
    for (auto &v : points) {
        out << (v-min)/(max-min) << std::endl;
    }
}

#endif //SLICEDOPTIM_IOPOINTSET_H
