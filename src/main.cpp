#include "solver.h"

#include <iostream>

#include <vector>

//#include <malloc.h>
#include <string>

#include <sys/time.h>
#include <sys/resource.h>

using namespace std;


int main(int argc, char *argv[]) {

  string input_file;
  Solver theSolver;


  if (argc <= 1) {
    cout << "Usage: ganak [options] [CNF_File]" << endl;
    cout << "Options: " << endl;
    cout << "\t -noPP  \t turn off preprocessing" << endl;
    cout << "\t -q     \t quiet mode" << endl;
    cout << "\t -t [s] \t set time bound to s seconds" << endl;
    cout << "\t -noCC  \t turn off component caching" << endl;
    cout << "\t -cs [n]\t set max cache size to n MB" << endl;
    cout << "\t -noIBCP\t turn off implicit BCP" << endl;
    cout << "\t -pol   \t Polarity: true, false, default, cache, reversecache, random, newrandom" << endl;
    cout << "\t -res [n]  \t do restart after n decisions" << endl;
    cout << "\t -act \t\t allow activity decrease CSVSADS" << endl;
    cout << "\t -seed [n] \t set random seed to n " << endl;
    cout << "\t -rand \t\t allow randomness "<< endl;
    cout << "\t -is file \t use independent support provided in file" << endl;
    cout << "\t -decc [n] \t use independent support after n decisions "<<endl;
    cout << "\t -maxdec [n] \t terminate after n decision "<<endl;
    cout << "\t -pcc \t\t use probablistic component caching "<< endl;
    cout << "\t -p  \t\t calculate projected model count on the set provided with -is" << endl;
    cout << "\t" << endl;

    return -1;
  }

  for (int i = 1; i < argc; i++) {
    if (strcmp(argv[i], "-noCC") == 0)
      theSolver.config().perform_component_caching = false;
    if (strcmp(argv[i], "-noIBCP") == 0)
      theSolver.config().perform_failed_lit_test = false;
    if (strcmp(argv[i], "-noPP") == 0)
      theSolver.config().perform_pre_processing = false;
    if (strcmp(argv[i], "-usegraphdecom") == 0)
      theSolver.config().usegraphdecom = true;
    if (strcmp(argv[i], "-pcc") == 0){
      cout << "Using probablistic component caching "<< endl;
      theSolver.config().usepcc = true;
    }
    if (strcmp(argv[i], "-res") == 0){
      if (argc <= i + 1) {
        cout << " wrong parameters" << endl;
        return -1;
      }
      else{
        theSolver.config().allowrestart = true;
        theSolver.config().restartafterdecisions = atol(argv[i + 1]);
      }
    }
    if (strcmp(argv[i], "-decc") == 0){
      if (argc <= i + 1) {
        cout << " wrong parameters" << endl;
        return -1;
      }
      else{
        theSolver.config().decc = atol(argv[i + 1]);
      }
    }
    if (strcmp(argv[i], "-maxdec") == 0){
      if (argc <= i + 2) {
        cout << " wrong parameters" << endl;
        return -1;
      }
      else{
        theSolver.config().maxdec = atol(argv[i + 1]);
        theSolver.config().minconflicts_ = atol(argv[i + 2]);
        theSolver.config().maxdecterminate = true;
        cout << "Using Max dec termination with maximum decisions "<< theSolver.config().maxdec 
        << " and number of conflicts should be less than " <<theSolver.config().minconflicts_<< endl;
      }
    }
    if (strcmp(argv[i], "-act") == 0){
      theSolver.config().allowactivitydecrease = true;
    }
    if (strcmp(argv[i], "-2random") == 0){
      theSolver.config().tworandom = true;
    }
    if (strcmp(argv[i], "-rand") == 0){
      theSolver.config().allowrandomness = true;
    }
    if (strcmp(argv[i], "-seed") == 0){
      if (argc <= i + 1){
        cout << " wrong parameters" << endl;
        return -1;
      }
      else{
        theSolver.config().randomseed = atol(argv[i + 1]);
      }
    }    
    else if (strcmp(argv[i], "-q") == 0)
      theSolver.config().quiet = true;
    else if (strcmp(argv[i], "-p") == 0)
      theSolver.config().projected = true;
    else if (strcmp(argv[i], "-v") == 0)
      theSolver.config().verbose = true;
    else if (strcmp(argv[i], "-t") == 0) {
      if (argc <= i + 1) {
        cout << " wrong parameters" << endl;
        return -1;
      }
      theSolver.config().time_bound_seconds = atol(argv[i + 1]);
      if (theSolver.config().verbose)
        cout << "time bound set to" << theSolver.config().time_bound_seconds << "s\n";
     } else if (strcmp(argv[i], "-pol") == 0) {
      if (argc <= i + 1) {
        cout << " must give polarity type" << endl;
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
      if (strcmp(argv[i + 1], "cache") == 0) {
          theSolver.config().polarity_config = polar_cached;
          found = true;
      }
      if (strcmp(argv[i + 1], "reversecache") == 0) {
          theSolver.config().polarity_config = reverse_cached;
          found = true;
      }
      if (strcmp(argv[i + 1], "random") == 0) {
          theSolver.config().polarity_config = random_pol;
          found = true;
      }
      if (strcmp(argv[i + 1], "newrandom") == 0) {
          theSolver.config().polarity_config = new_random_pol;
          found = true;
      }
      if (!found) {
          cout << "ERROR: The option to '-pol' you gave, '" << argv[i + 1] << "' cannot be parsed" << endl;
          exit(-1);
      }

     } else if (strcmp(argv[i], "-cs") == 0) {
      if (argc <= i + 1) {
        cout << " wrong parameters" << endl;
        return -1;
      }
      theSolver.statistics().maximum_cache_size_bytes_ = atol(argv[i + 1]) * (uint64_t) 1000000;
    }
    else if (strcmp(argv[i], "-cdb") == 0) {
      if (argc <= i + 1) {
        cout << " wrong parameters" << endl;
        return -1;
      }
      theSolver.config().cleardbafter = atol(argv[i + 1]);
    }
    else if (strcmp(argv[i], "-is") == 0){
      if (argc <= i+1){
        cout << " wrong parameters" << endl;
        return -1;
      }
      theSolver.config().useindependentsupport = true;
      theSolver.config().independent_support_file = argv[i+1];
    } 
    else if (strcmp(argv[i], "-poset") == 0){
      // if (argc <= i+1){
      //   cout << " wrong parameters" << endl;
      //   return -1;
      // }
      theSolver.config().useposet = true;
      // theSolver.config().posetfile = argv[i+1];
    }
    else if (strcmp(argv[i], "-cacheten") == 0) {
      theSolver.config().usecachetencoding = true;
    }
    else if (strcmp(argv[i], "-isoCC") == 0) {
      if (argc <= i+1){
        cout << " wrong parameters" << endl;
        return -1;
      }
      theSolver.config().useIsomorphicComponentCaching = true;
      theSolver.config().isomorphicComponentCachingn = atoi(argv[i+1]);
    } 
    else
      input_file = argv[i];
  }

  theSolver.solve(input_file);

//  cout << sizeof(LiteralID)<<"MALLOC_STATS:" << endl;
//  malloc_stats();

//  rusage ru;
//  getrusage(RUSAGE_SELF,&ru);
//
//   cout << "\nRus: " <<  ru.ru_maxrss*1024 << endl;
//  cout << "\nMALLINFO:" << endl;
//
//  cout << "total " << mallinfo().arena + mallinfo().hblkhd << endl;
//  cout <<  mallinfo().arena << "non-mmapped space allocated from system " << endl;
//  cout <<  mallinfo().ordblks << "number of free chunks " << endl;
//  cout <<  mallinfo().smblks<< "number of fastbin blocks " << endl;
//  cout <<  mallinfo().hblks<< " number of mmapped regions " << endl;
//  cout <<  mallinfo().hblkhd<< "space in mmapped regions " << endl;
//  cout <<  mallinfo().usmblks<< " maximum total allocated space " << endl;
//  cout <<  mallinfo().fsmblks<< "space available in freed fastbin blocks " << endl;
//  cout <<  mallinfo().uordblks<< " total allocated space " << endl;
//  cout <<  mallinfo().fordblks<< "total free space " << endl;
//  cout <<  mallinfo().keepcost<< " top-most, releasable (via malloc_trim) space " << endl;
  return 0;
}
