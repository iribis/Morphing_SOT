//
//

#ifndef SLICEDOPTIM_NBALLSLICEDOPTIMALTRANSPORT_H
#define SLICEDOPTIM_NBALLSLICEDOPTIMALTRANSPORT_H

#include <vector>
#include <iostream>
#include <random>
#include <algorithm>
#include <cstring>
#include "NBallRadonManager.h"
#include "../Math/VecX.h"
#include "../Math/myMath.h"
#include "../Tools/iopointset.h"
#include "../Tools/mapping.h"
#include "../Tools/my_utility.h"

NBallRadonManager nbrm(3);
const double PI  =3.141592653589793238463;

template <class VECTYPE>
void project(const std::vector<VECTYPE>& points, std::vector<std::pair<double, int>> &pointsProject, const VECTYPE& dir,int N){
    for(int s = 0; s < points.size(); ++s){
        VECTYPE p = points[s];
        double proj = (p * dir);
        std::pair<double, int> val_indice = std::pair<double, int>(proj, s);
        pointsProject.push_back(val_indice);
    }

}

template <class VECTYPE>
inline void getInverseRadonNCube(int D, int nbSamples, std::vector<double>& pos,VECTYPE dir,std::vector<VECTYPE>& pointsTarget){

    pos.resize(nbSamples);
    double mul = pointsTarget.size()/nbSamples;
    
    std::vector<double> obj_projected;
    obj_projected.resize(pointsTarget.size());
    for (int i = 0; i < pointsTarget.size(); i++)
    {
        obj_projected[i] = pointsTarget[i]*dir;
        //std::cout<<pointsTarget[i]<<std::endl;
    }
    

    std::sort(obj_projected.begin(),obj_projected.end());
    for (int i = 0; i < nbSamples; i++)
    {
        pos[i] = obj_projected[int(i*mul + mul/2)];
        //std::cout<<pos[i]<<std::endl;
    }
}


template <class VECTYPE>
inline bool testPoissonND(const VECTYPE& v, typename std::vector<VECTYPE>::const_iterator begin, typename std::vector<VECTYPE>::const_iterator end, double radius){
    bool test = true;
    for (typename std::vector<VECTYPE>::const_iterator it = begin; it != end && test; ++it){
        test = (v * (*it) < radius);
    }
    return test;
}

/**
 * Choose \p m directions in N dimension N being defined by the dimention of the content of directions.
 * Two selection methods are available. Either the direction are uniformly selected in ND
 * or one can force half of the them to lie in 2D planes to optimize the projection repartition as well.
 *
 * @param directions Table of directions to output. Must be initialized with \p m VECTYPE of the disired dimension
 * @param m Number of directions to pick
 * @param seed Seed for the random generator. Only applied once
 * @param projective If true half the directions will lie in 2D planes.
 */
template <class VECTYPE>
inline void chooseDirectionsND(std::vector<VECTYPE>& directions, int m, int seed){

    static std::mt19937 generatorND(seed);
    static std::normal_distribution<>normalND;
    static std::uniform_real_distribution<double> unif(0, 1);

    int dim = directions.front().dim();

    double pradius = std::cos(0.99 * std::pow(nsphereArea(1, dim) / (m * nballVolume(1, dim - 1)), 1. / (dim - 1)));
    for (int k = 0; k < m; ++k){
        if(dim == 2){
            double theta =  (float(k)/m + unif(generatorND)/float(m))*2*PI; // stratified 2D directions
            directions[k][0] = cos(theta);
            directions[k][1] = sin(theta);

        }else{
            for (int j = 0; j < dim; ++j){
                directions[k][j] = normalND(generatorND);
            }
            
        }   
        directions[k].normalize();
    }
}

/**
 * Compute optimal transport in 1D for direction \f$ \theta \f$ and \f$d_{j}\f$ being the 1D displacement of \f$\x^j\f$
 * that minimize the 1D sliced optimal transport along \f$ \theta \f$.
 *
 * Denoting $\sigma$ the permutations of the indices \f$\{j\}_{j=1..N}\f$ such that
 * \f$\bigl(\x^{\sigma(j)} \cdot \theta \bigr)_j\f$, is a sorted sequence of increasing values,
 * one can compute \f$d_{j}\f$ via
 * \f$ d_{j} = C_{\theta}^{-1}\left(\frac{\sigma(j)-\frac12}{N}\right)\,. \vspace*{-1mm}\f$
 *
 * @param dir Direction \f$ \theta \f$
 * @param points Table containing the points \f$ x_j \f$
 * @param pos Table containing the optimal solution in 1D
 * @param shift Output the 1D shift to apply to minimize transport cost
 * @param pointsProject Memory buffer used to store the 1D projection of points. Must be the same size as \p points
 */
template<class VECTYPE>
inline void slicedStepNBall(const VECTYPE& dir,
                            const std::vector<VECTYPE>& points,
                            std::vector<VECTYPE>& pointsTarget,
                            std::vector<double>& shift)
{
    int N = points.front().dim();
    
    float rnd = ((float)rand() / RAND_MAX);

    std::vector<std::pair<double, int>> pointsProject;
    project(points, pointsProject, dir, N);

    std::vector<double> pos(pointsProject.size());
    getInverseRadonNCube(N,pointsProject.size(), pos,dir, pointsTarget);

    std::sort(pointsProject.begin(), pointsProject.end());

    //Computes required shift to optimize 1D optimal transport
    for (int i = 0; i < pointsProject.size(); ++i) {
        //Compute shifting
        double s = pos[i] - pointsProject[i].first;
        shift[pointsProject[i].second] += s;
    }

}

/**
 * Compute optimal transport in 1D for the \p directions and displace \f$ x_j \f$ by
 * \f$\pmb{\delta}^j \EqDef \frac{1}{K} \sum_{i=1}^K d_{i,j}\, \theta_i \vspace*{-1mm}\f$ with
 * with \f$ d_{i,j} \f$ being the displacement of \f$\x^j\f$ that minimize the 1D sliced optimal
 * transport along direction \f$ \theta_i \f$
 *
 * @param pointsOut Table containing the points \f$ x_j \f$
 * @param directions Table containing the \f$\theta_i\f$
 * @param pos Table containing the 1D positions optimizing transport in 1D
 * @param shift Used to avoid having to allocate uge chunks of memory. Must be a vector of size m containing vectors of same size as \p pointsOut.
 * @param finalShift Used to avoid having to allocate huge chunks of memory Must be a vector of same size as \p pointsOut.
 * @param pointsProject Used to avoid having to allocate uge chunks of memory. Must be a vector of size m containing vectors of same size as \p pointsOut.
 * @return the Wasserstein cost of the current iteration.
 */
template <class VECTYPE>
inline void slicedOptimalTransportBatch(std::vector<VECTYPE>& pointsOut,
                                 std::vector<VECTYPE>& pointsTarget,
                                 const std::vector<VECTYPE>& directions,
                                 std::vector<std::vector<double>>& shift,
                                 std::vector<VECTYPE>& finalShift)
 {

    int m = directions.size();
    int nbPoints = pointsOut.size();

    //Compute the shift along each direction
#pragma omp parallel for shared(directions, shift)
    for (int k = 0; k < m; ++k){
        for(double& v : shift[k]){
			v = 0.;
        }
        const VECTYPE& dir = directions[k];

        slicedStepNBall(dir, pointsOut,pointsTarget, shift[k]);
    }

        //Accumulate shift from all directions
#pragma omp parallel for
    for (int i = 0; i < nbPoints; ++i) {
        VECTYPE sh(finalShift[i].dim());
        memset(&sh[0], 0, finalShift[i].dim() * sizeof(sh[0]));
        for (int k = 0; k < m; ++k) {
            sh += shift[k][i] * directions[k];
        }
        finalShift[i] = sh;
        finalShift[i] /= m;
    }

    //Displace points according to accumulated shift
#pragma omp parallel for
    for (int i = 0; i < nbPoints; ++i) {
        pointsOut[i] += finalShift[i];
    }

}

void print_progress(double ratio){
    int barWidth = 70;

    std::cout << "[";
    int pos = barWidth * ratio;
    for (int i = 0; i < barWidth; ++i) {
        if (i < pos) std::cout << "=";
        else if (i == pos) std::cout << ">";
        else std::cout << " ";
    }
    std::cout << "] " << int(ratio * 100.0) << " %\r";
    std::cout.flush();
}

/**
 *  \brief Computes an optimized point set to uniformly sample the unit N-Ball using sliced optimal transport
 *
 * @param pointsIn Contains input ND points
 * @param pointsOut Contains optimized ND points
 * @param nbIter Number of iterations
 * @param m Number of slice per iteration
 * @param seed random seed
 * @return the Sliced Wasserstein distance between the samples and the uniform distribution.
 */

template <class VECTYPE>
inline void slicedOptimalTransportNBall(const std::vector<VECTYPE>& pointsIn,
                                        std::vector<VECTYPE>& pointsOut,
                                        std::vector<VECTYPE>& pointsTarget,
                                        int nbIter,
                                        int m,
                                        int seed, bool silent = false)
{
    
    int N = pointsIn.front().dim();
    pointsOut = pointsIn;
    nbrm = NBallRadonManager(N);
    //Accumulation shift to be applied later on
    std::vector<std::vector<double>> shift(m, std::vector<double>(pointsOut.size()));
    std::vector<VECTYPE> finalShift(pointsOut.size(), VECTYPE(N));

    std::vector<VECTYPE> directions(m, VECTYPE(N));

    //Iterate 1D Optimal Transport
    for (int i = 0; i < nbIter; i += 1){
        if(!silent){
            print_progress(double(i)/nbIter);
        }
        chooseDirectionsND(directions, m, seed);

        slicedOptimalTransportBatch(pointsOut,pointsTarget, directions, shift, finalShift);
    }
    print_progress(1.0);

}

/**
 * Compute the 1D transport cost of the from given ND \param points
 * to the distribution whose optimal 1D positions along dir are given in \param pos
 *
 * @tparam VECTYPE any vector type with * being the dot product and [] granting access to the coordinates value
 * @param points ND points
 * @param dir Projection direction
 * @param pos Sorted optimal 1D positions for target distribution along dir
 * @param pointsProject Pre allocated buffer used to sort the points's projection
 */
template <class VECTYPE>
inline double computeSlicedOTCost(const std::vector<VECTYPE>& points,
                                  const VECTYPE& dir,
                                  const std::vector<double>& pos,
                                  std::vector<double>& pointsProject
                                  ){
    double cost = 0.;

    for (int i = 0; i < pointsProject.size(); ++i){
        pointsProject[i] = points[i] * dir;
    }
    sort(pointsProject.begin(), pointsProject.end());

    //Place them at optimal places
    for (int i = 0; i < pointsProject.size(); ++i) {
        cost += (pointsProject[i] - pos[i]) * (pointsProject[i] - pos[i]);
    }
    return cost / points.size();
}

#endif //SLICEDOPTIM_NBALLSLICEDOPTIMALTRANSPORT_H
