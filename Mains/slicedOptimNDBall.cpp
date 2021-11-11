//
//

#include <iostream>
#include <iomanip>
#include <vector>
#include <cmath>
#include <string>
#include <random>
#include <cstdlib>
#include "../Math/VecX.h"
#include "../Tools/iopointset.h"
#include "../Tools/mapping.h"
#include "../Transport/slicedOptimalTransportNBall.h"

#define DIM 3

using namespace std;

void usage(const char **argv){
    cerr << argv[0] << " [-o <OutputFileName>] "
                       "[-n <nbPoints>] [-f <objFileName>] [-m <nbRealisations>] [-p <nbIteration>]"
                       "[--step <nbDirectionPerStep>] [-s <seed>] [-d <dimension>] [-c]" << endl;
}

void handleParameters(int argc,
                      const char** argv,
                      string& outPrefix,
                      int& nbIter,
                      int& m,
                      int& p,
                      int& seed,
                      int& nbPointsets,
                      int& dim,
                      bool& cubify,
                      bool& silent,
                      std::string& filename){
    int i = 1;
    while (i < argc){
        if (!strncmp(argv[i], "-o", 2)) {
            outPrefix = (argv[i+1]);
            ++i;
        } else if (!strncmp(argv[i], "-n", 2)) {
            p = atoi(argv[i+1]);
            ++i;
        } else if (!strncmp(argv[i], "--step", 6)) {
            m = atoi(argv[i+1]);
            ++i;
        } else if (!strncmp(argv[i], "-m", 2)) {
            nbPointsets = atoi(argv[i+1]);
            ++i;
        } else if (!strncmp(argv[i], "-p", 2)) {
            nbIter = atoi(argv[i+1]);
            ++i;
        } else if (!strncmp(argv[i], "-s", 2)) {
            seed = atoi(argv[i+1]);
            ++i;
        }else if (!strncmp(argv[i], "-f", 2)) {
            filename = argv[i+1];
            ++i;
        } else if (!strncmp(argv[i], "-d", 2)) {
            dim = atoi(argv[i+1]);
            ++i;
        } else if (!strncmp(argv[i], "-c", 2)) {
            cubify = true;
        } else if (!strncmp(argv[i], "-h", 2) || !strncmp(argv[i], "--help", 6)) {
            cerr << "Help: " << endl;
            cerr << "Option list:" << endl;
            cerr << "\t-o <OutputFileName> (optional): Specifies an output file in which points will be written."
                 << "If unset standard output will be used" << endl;
            cerr << "\t-n <nbPoints> (default 1024): Specifies the number of points to generate" << endl;
            cerr << "\t-m <nbRealisations> (default 1): Specifies the number of generated pointsets" << endl;
            cerr << "\t-p <nbIteration> (default 4096): Specifies the number of batches in the optimization process" << endl;
            cerr << "\t--step <nbDirectionPerStep> (default 32): Specifies the number of slices per batch in the optimization process" << endl;
            cerr << "\t-s <seed> (default 133742): Specifies the random seed" << endl;
            cerr << "\t-d <dimension> (default 2): Specifies samples dimension" << endl;
            cerr << "\t-c (optional): If unset points will be given in the unit ball. Else they will be in the unit cube [0,1)^d" << endl;
            cerr << "\t--silent (optional): Cancels all outputs other than the points and errors" << endl;
            cerr << "\t" << endl;
            usage(argv);
            exit(2);
        } else if (!strncmp(argv[i], "--silent", 8)) {
            silent = true;
        } else {
            cerr << "Unknown option " << argv[i] << endl;
            exit(1);
        }
        ++i;
    }
}

template <class VECTYPE>
void loadOBJ(std::string& filename,vector<VECTYPE>& points){
    FILE * file = fopen(filename.c_str(), "r");
    if( file == NULL ){
        printf("Impossible to open the file !\n");
        return;
    }

    while( 1 ){

        char lineHeader[128];
        // read the first word of the line
        int res = fscanf(file, "%s", lineHeader);
        if (res == EOF)
            break; // EOF = End Of File. Quit the loop.

        // else : parse lineHeader
        if ( strcmp( lineHeader, "v" ) == 0 ){
            VECTYPE vertex;
            float f1,f2,f3;
            fscanf(file, "%f %f %f\n", &f1, &f2, &f3 );
            vertex[0] = f1;
            vertex[1] = f2;
            vertex[2] = f3;
            points.push_back(vertex);
        }
    }
    std::cout << points.size() << std::endl;
}

template <class VECTYPE>
int main_template(int argc, const char **argv) {

    int nbIter = 4096;
    int dim = DIM;
    int p = 1024;
    int m = 64;
    int nbPointsets = 1;
    //Default parameters value
    string outPrefix = "";
    int seed = 133742;
    bool cubify = false;
    bool silent = false;
    std::string filename = "./test.obj";

    handleParameters(argc, argv, outPrefix, nbIter, m, p, seed, nbPointsets, dim, cubify, silent,filename);

    //If file name ends in .bin then the output will be written in binary
    ostream* out = &cout;
    if (outPrefix != "") {
        out = new ofstream(outPrefix);
        if (out->fail()) {
            cerr << "Can't open output file \"" + outPrefix + "\"" << endl;
            exit(3);
        }
        cerr << "Output file: " + outPrefix + ".dat" << endl;
    }

    mt19937 generator(seed);
    if (!silent) {
        cerr << "Generating " << nbPointsets << " sets of " << p << " points in " << dim << "D using " << nbIter << " batches of " << m
             << " slices" << endl;
        cerr << "Map points to cube: " << (cubify ? "True" : "False") << endl;
    }

    for (int indPointset = 0; indPointset < nbPointsets; ++indPointset) {

        vector<VECTYPE> points(p, VECTYPE(dim));
        vector<VECTYPE> obj_point_list;
        
        loadOBJ(filename,obj_point_list);

        //Init from whitenoise in ball
        uniform_real_distribution<> unif(0., 1.0);
        normal_distribution<> norm(0., 1.0);

        for (size_t i = 0; i < points.size(); ++i) {
            points[i] = randomVectorInCube<VECTYPE>(dim, generator);
        }

        vector<VECTYPE> result;

        slicedOptimalTransportNBall(points, result,obj_point_list, nbIter, m, seed);

        if (indPointset != 0)
            *out << "#" << endl;

        savePointsetND(*out, result);
    }

    return 0;
}

#include "../Tools/dimensionsInstantiation.hpp"
