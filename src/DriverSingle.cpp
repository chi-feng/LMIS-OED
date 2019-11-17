#include <memory>
#include <string>
#include <iostream>
#include <iomanip>
#include <fstream>

#include <Eigen/Core>
#include "MultivariateT.h"
#include "ExpectedInformationEstimator.h"
#include "SimpleExperiment.h"
#include "Utilities.h"
#include "RandomGenerator.h"

using namespace std;
using namespace Eigen;
using namespace Utilities;

int main(int argc, char **argv)
{
    RandomGenerator::Initialize();

    int dim         = GetOption<int>(argc, argv, "-dim", 4);
    int poi         = GetOption<int>(argc, argv, "-poi", 1);
    double sigeps   = GetOption<double>(argc, argv, "-sigeps", 0.4);
    double k        = GetOption<double>(argc, argv, "-k", 5.0);
    auto experiment = make_shared<SimpleExperiment>(dim, poi, sigeps, k);

    string   designStr = GetOption<string>(argc, argv, "-design", "0.8");
    VectorXd design    = StringToVectorXd(designStr);

    auto estimator = make_shared<ExpectedInformationEstimator>(experiment);
    estimator->N                       = GetOption<int>(argc, argv, "-N", 1000);
    estimator->M1                      = GetOption<int>(argc, argv, "-M1", 100);
    estimator->M2                      = GetOption<int>(argc, argv, "-M2", 100);
    estimator->maxComponents           = GetOption<int>(argc, argv, "-maxComponents", estimator->N);
    estimator->useIS                   = GetOption<bool>(argc, argv, "-useIS", false);
    estimator->useMIS                  = GetOption<bool>(argc, argv, "-useMIS", true);
    estimator->useExactPosterior       = GetOption<bool>(argc, argv, "-useExactPosterior", false);
    estimator->useSortHeuristic        = GetOption<bool>(argc, argv, "-useSortHeuristic", true);
    estimator->useReverseLikelihood    = GetOption<bool>(argc, argv, "-useReverseLikelihood", true);
    estimator->useMinSampleDistance    = GetOption<bool>(argc, argv, "-useMinSampleDistance", false);
    estimator->useBiasCorrection       = GetOption<bool>(argc, argv, "-useBiasCorrection", false);
    estimator->verbose                 = GetOption<bool>(argc, argv, "-verbose", false);
    estimator->biasingDistributionType = GetOption<string>(argc, argv, "-biasingDistributionType", "MVT");
    estimator->dof                     = GetOption<double>(argc, argv, "-dof", 2.5);
    estimator->nugget                  = GetOption<double>(argc, argv, "-nugget", 1e-3);
    estimator->debugCondLikelihood     = GetOption<bool>(argc, argv, "-debugCondLikelihood", false);
    estimator->useMarginal             = GetOption<bool>(argc, argv, "-useMarginal", false);

    cout << estimator << endl;

    double   EIG = estimator->Evaluate(design);

    string   dumpfile = GetOption<string>(argc, argv, "-dumpfile", "none");
    if (dumpfile != "none") {
      estimator->WriteToFile(dumpfile, design);
    }
    
    return 0;
}

