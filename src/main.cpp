#include "solver.h"
#include "GitSHA1.h"

#include <iostream>

#include <vector>

#include <string>

#include <sys/time.h>
#include <sys/resource.h>

using namespace std;


int main(int argc, char *argv[]) {

  string input_file;
  Solver theSolver;

  cout << "c Outputting solution to console" << endl;
  cout << "c GANAK version 1.0.0" << endl;

  if (argc <= 1) {
    cout << "Usage: ganak [options] [CNF_File]" << endl;
    cout << "Options: " << endl;
    cout << "\t -noPP  \t\t turn off preprocessing" << endl;
    cout << "\t -q     \t\t quiet mode" << endl;
    cout << "\t -t [s] \t\t set time bound to s seconds" << endl;
    cout << "\t -noCC  \t\t turn off component caching" << endl;
    cout << "\t -cs [n]\t\t set max cache size to n MB" << endl;
    cout << "\t -noIBCP\t\t turn off implicit BCP" << endl;
    cout << "\t -noPCC\t\t\t turn off probabilistic component caching" << endl;
    cout << "\t -seed [n]\t\t set random seed to n (Default: 1000)" << endl;
    cout << "\t -m [n] \t\t set the range of hash function (= 64 x n) (Default: 1) "<<endl;
    cout << "\t -delta [n] \t\t set the confidence parameter to n (Default: 0.05) "<<endl;
    cout << "\t -noCSVSADS\t\t turn off CSVSADS variable branching heuristic" << endl;
    cout << "\t -pol [Polarity]\t Polarity: true, false, default, polaritycache (Default: polaritycache)" << endl;
    cout << "\t -EDR\t\t\t turn on EDR variable branching heuristic" << endl;
    cout << "\t -LSO [n]\t\t learn and start over after n decisions (Default: 5000)" << endl;
    cout << "\t -noPMC\t\t\t turn off projected model counting " << endl;
    cout << "\t -maxdec [n] [m] \t terminate after n decision if conflict is less than m "<<endl;
    cout << "\t" << endl;
    return -1;
  }

  for (int i = 1; i < argc; i++) {
    if (strcmp(argv[i], "-noCC") == 0)
      theSolver.config().perform_component_caching = false;
    else if (strcmp(argv[i], "-noIBCP") == 0)
      theSolver.config().perform_failed_lit_test = false;
    else if (strcmp(argv[i], "-noPP") == 0)
      theSolver.config().perform_pre_processing = false;
    else if (strcmp(argv[i], "-q") == 0)
      theSolver.config().quiet = true;
    else if (strcmp(argv[i], "-v") == 0)
      theSolver.config().verbose = true;
    else if (strcmp(argv[i], "-noPCC") == 0)
      theSolver.config().perform_pcc = false;
    else if (strcmp(argv[i], "-noCSVSADS") == 0)
      theSolver.config().use_csvsads = false;
    else if (strcmp(argv[i], "-noPMC") == 0)
      theSolver.config().perform_projectedmodelcounting = false;
    else if (strcmp(argv[i], "-EDR") == 0){
      theSolver.config().use_csvsads = false;
      theSolver.config().use_edr = true;
    }
    else if (strcmp(argv[i], "-LSO") == 0){
      if (argc <= i + 1){
        cout << "ERROR: wrong parameters" << endl;
        return -1;
      }
      else{
        theSolver.config().lsoafterdecisions = atol(argv[i + 1]);
        if (!theSolver.config().lsoafterdecisions)
          theSolver.config().use_lso = false;
      }
    }
    else if (strcmp(argv[i], "-seed") == 0){
      if (argc <= i + 1){
        cout << "ERROR: wrong parameters" << endl;
        return -1;
      }
      else{
        theSolver.config().randomseed = atol(argv[i + 1]);
      }
    }
    else if (strcmp(argv[i], "-m") == 0){
      if (argc <= i + 1){
        cout << "ERROR: wrong parameters" << endl;
        return -1;
      }
      else{
        theSolver.config().hashrange = atol(argv[i + 1]);
        cout << "c The value of hashrange is 64x"<<theSolver.config().hashrange<< endl;
      }
    }
    else if (strcmp(argv[i], "-delta") == 0){
      if (argc <= i + 1){
        cout << "ERROR: wrong parameters" << endl;
        return -1;
      }
      else{
        theSolver.config().delta = stof(argv[i + 1]);
        cout << "c The value of delta is "<< theSolver.config().delta << endl; 
      }
    }
    else if (strcmp(argv[i], "-pol") == 0) {
      if (argc <= i + 1) {
        cout << "ERROR: must give polarity type" << endl;
        return -1;
      }
      bool found = false;
      if (strcmp(argv[i + 1], "true") == 0) {
          theSolver.config().polarity_config = polar_true;
          found = true;
      }
      if (strcmp(argv[i + 1], "false") == 0) {
          theSolver.config().polarity_config = polar_false;
          found = true;
      }
      if (strcmp(argv[i + 1], "default") == 0) {
          theSolver.config().polarity_config = polar_default;
          found = true;
      }
      if (strcmp(argv[i + 1], "polaritycache") == 0) {
          theSolver.config().polarity_config = polaritycache;
          found = true;
      }
      if (!found) {
          cout << "ERROR: The option to '-pol' you gave, '" << argv[i + 1] << "' cannot be parsed" << endl;
          exit(-1);
      }
    }
    else if (strcmp(argv[i], "-t") == 0) {
      if (argc <= i + 1) {
        cout << "ERROR: wrong parameters" << endl;
        return -1;
      }
      theSolver.config().time_bound_seconds = atol(argv[i + 1]);
      if (theSolver.config().verbose)
        cout << "c time bound set to" << theSolver.config().time_bound_seconds << "s" << endl;
     }
    else if (strcmp(argv[i], "-cs") == 0) {
      if (argc <= i + 1) {
        cout << "ERROR: wrong parameters" << endl;
        return -1;
      }
      theSolver.statistics().maximum_cache_size_bytes_ = atol(argv[i + 1]) * (uint64_t) 1000000;
    }
    else if (strcmp(argv[i], "-maxdec") == 0){
      if (argc <= i + 2) {
        cout << "ERROR: wrong parameters" << endl;
        return -1;
      }
      else{
        theSolver.config().maxdec = atol(argv[i + 1]);
        theSolver.config().minconflicts_ = atol(argv[i + 2]);
        theSolver.config().maxdecterminate = true;
      }
    }else
      input_file = argv[i];
  }


  cout << "c ganak GIT revision: " << Ganak::get_version_sha1() << endl;
  cout << "c ganak build env: " << Ganak::get_compilation_env() << endl;
  theSolver.solve(input_file);
  return 0;
}
