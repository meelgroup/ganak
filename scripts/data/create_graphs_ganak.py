#!/bin/python3

import os
import sqlite3
import re
import nbformat as nbf


def convert_to_cactus(fname, fname2):
    # print("fname:" , fname)
    f2 = open(fname2, "w")
    f = open(fname, "r")
    text = f.read()
    mylines = text.splitlines()
    i = 0;
    time = []
    for line in mylines:
      time.append(float(line.split()[0]))
      i += 1

    lastnum = -1
    for a in range(0, 3600, 1):
      num = 0
      for t in time:
        #print "t: %f a: %d" %(t, a)
        if (t < a) :
          num += 1

      if (lastnum != num):
          f2.write("%d \t%d\n" %(num, a))
      lastnum = num
    f.close()
    f2.close()
    return len(mylines)


def get_versions():
    vers = []
    con = sqlite3.connect("mydb.sql")
    cur = con.cursor()
    res = cur.execute("""
                      SELECT ganak_ver
                      FROM data
                      where ganak_ver is not NULL and ganak_ver != '' group by ganak_ver""")
    for a in res:
        vers.append(a[0])

    con.close()
    return vers

def get_dirs(ver : str):
    ret = []
    con = sqlite3.connect("mydb.sql")
    cur = con.cursor()
    res = cur.execute("SELECT dirname, ganak_call FROM data where ganak_ver='"+ver+"' group by dirname")
    for a in res:
        call = a[1]
        call = re.sub("././ganak", "", call)
        call = re.sub(" mc2022.*cnf.*", "", call)
        ret.append([a[0], call])

    con.close()
    return ret

def gnuplot_name_cleanup(name: str) -> str:
    # remove all non-alphanumeric characters except for underscores and dashes
    name = re.sub(r'\"', '', name)
    # replace multiple underscores or dashes with a single one
    name = re.sub(r'_', '=', name)
    print(name)
    return name

versions = get_versions()
fname2_s = []
# not_calls = ["ExactMC"]
not_calls = []
# not_versions = ["sharpsat", "gpmc", "6368237b"]
# only_calls = ["--ignore 1 --arjun 1 --maxcache 3500 --vivif 1 --decide 2 --sbva 1000"]
# not_versions = ["ganak"]
# not_calls = ["forcebranch", "target"] # "cachetime"
# exactm: out-ganak-6318929.pbs101-5/
# exactmc2: out-ganak-6328707.pbs101-7
# sharpsat: out-ganak-6318929.pbs101-7

# note: sharpsat thinks new track2 is all UNSAT, likely due to weights. maybe run unweighted.

only_dirs = [
            # "out-ganak-6683080.pbs101-1/", # best of GANAK so far, lbd 1
            # "out-ganak-6396805.pbs101-1/",  # best, but no cache
            # "out-ganak-6396805.pbs101-12/", # used to be best
            # "out-ganak-6318929.pbs101-4/", # approxmc
            # "out-ganak-6749880.pbs101-0/", # like out-ganak-6683080.pbs101-1, but cluster is slower(!!)
            # "out-ganak-6683080.pbs101-4/", # last run, without refactored Arjun+fixed SAT
            # "out-ganak-6841576.pbs101-0/", # trying TD exponent, turning off gates for Arjun
            # "out-ganak-6853697.pbs101-4/", # backbone first, if proj is smaller than TD, td=0.1, turn off gates -- ran with CMS 58286e58ef64156ed2c16ec4104820c7fb4aeda5
            # "out-ganak-6867174.pbs101-0", # fix with max grow 6, also try no freq
            # "out-ganak-6870225.pbs101-0", # max grow 6, maxtd 5
            # "out-ganak-6896385.pbs101-0", # supposedly new, not really exciting
            # "out-ganak-6907885.pbs101-", # BAD variable activity merged idea
            # "out-ganak-6909164.pbs101-", # BAD var activity merged, polarity changed
            # "out-ganak-6916780.pbs101-", # variable activity
            # "out-ganak-6910211.pbs101-2", # check diferent component sortings
            # "out-ganak-6917417.pbs101-", # sortings, lbd, polarity
            # "out-ganak-6918163.pbs101-0", # restart on mccomp 2023, 20k restarts
            # "out-ganak-6965806.pbs101", # mccomp 2022+3, restarts base etc
            # "out-ganak-7015279.pbs101", # new extension and symm run
            #######
            # "out-ganak-7021521.pbs101-0", # 2023 mccomp, restarts, fixed blocking lit, cube extend only no contract
            # "out-ganak-7028560.pbs101-", # 2022 + 2023
            # "out-ganak-7041485.", # kr-24 instances from Arijit
            # "out-ganak-7041554.pbs101-", # rapid restarts
            # "out-ganak-7047972.pbs101-0",
            # "out-ganak-7049225.pbs101-",
            # "out-ganak-7022833.pbs101-4", # best ever
            # "out-ganak-7048280.pbs101-0", # best ever now

            # unproj, i.e. track 1
            # "out-ganak-6318929.pbs101-5/", # exactmc
            # "out-ganak-6328707.pbs101-7/", # exactmc
            # "out-ganak-6318929.pbs101-7/", # sharpsat
            # no point in combining out-ganak-7178422.pbs101-0 with out-ganak-7184237.pbs101-2, 171 either way
            # "out-ganak-7178422.pbs101-2", # unproj-2023 16 GB d4
            # "out-ganak-7178422.pbs101-3", # unproj-2023 16 GB gpmc
            # "out-ganak-7178422.pbs101-", # unproj-2023 16 GB
            # "out-ganak-7180395.pbs101-", #unproj 2023 4GB, freq-128 (vsads)
            # "out-ganak-7184237.pbs101-2", # unproj-2023 4GB, vsads, checking vsads params
            # "out-ganak-7189905.pbs101-" # unproj 2023, 4GB, cache elem simplification, TD back to 0.1, SAT restarts
            # "out-ganak-7206355.pbs101", # unproj, new freq, not terrible with high TW weight
            # "out-ganak-7225470.pbs101-6",
            # "out-ganak-7225470.pbs101-7",
            # "out-ganak-7162995.pbs101-0", # new run, good. -- old freq setup
            # "out-ganak-7246958", # fixed cadiback, high tdmaxw
            # "out-ganak-7255018.pbs101-7", # tdmaxw 90 + higher tdminw. 12 minw seems to work well -- BEST (maxw 90, minw 12)
            # "out-ganak-7257776.pbs101-9", # higher tdmaxw, higher tdminw -- BEST (maxw 100, minw 15)
            #"out-ganak-7295009.pbs101-0", #fixed memout from contraction, fixed too much extend, distill-bin, BIG before backbone
            # "out-ganak-7306993.pbs101-0",# fixed cadical, fixed resolv-subs, fixed oracle, 2 new options for extra oracle & resolv-subs
                                         # best is "--tdminw 15 --tdmaxw 100 --resolvsub 0 --extraoracle 1"
            # BAD PARAMS "out-ganak-7316065.pbs101", # also appmc, play around with some opitions
            # "out-ganak-7320968.pbs101-1", # best: --tdminw 15 --tdmaxw 100, --appmct 2000
            # "out-ganak-7320968.pbs101-7", # best: --tdminw 15 --tdmaxw 100, --appmct 2000
            # "out-ganak-7320968.pbs101-7", # best: --tdminw 15 --tdmaxw 100 --varfreqdiv 30, --appmct 2000 (tdexpmul 1?)
            #"out-ganak-7334807.pbs101-7", # new run, some parameter tuning, fixed SAT restart bug (again, this time good) BEST appmc
            # "out-ganak-7377878.pbs101", # try different backbone setups, only oracle vivif, arjun probe, etc
            # old BEST: --tdmaxw 100 --varfreqdiv 25 --tdminw 15 --backbone 0 --arjunprobe 1
            # also: --tdexpmult 1.25
            # TODO: do a bit more probing, but not full backbone before indep minim
            # "out-ganak-7377878.pbs101-5", # best of the above
            # "out-ganak-7334807.pbs101-2", # new run, some parameter tuning, fixed SAT restart bug (again, this time good) BEST non-appmc
            # "out-ganak-7320968.pbs101-1", # best: --tdminw 15 --tdmaxw 100 --varfreqdiv 30, --appmct 2000 (tdexpmul 1?)
            # "out-ganak-7377878.pbs101-5", # best of the above
            # "out-ganak-7401120.pbs101-", # timeout in backbone 30s
            # "out-ganak-7419157.pbs101-3", # more simplification pre-backward, lower backward confl
            # "out-ganak-7419157.pbs101-", # more simplification pre-backward, lower backward confl
            # old BEST: --tdminw 15 --backbone 0 --arjunprobe 1 --arjunsimplev 2 --arjunbackwmaxc 20000
            # "out-ganak-7422239.pbs101-2" #unit cls to duplicated CNF, different cutoff for gates,
                    # some param tuning best: --tdminw 15 --backbone 0 --arjunprobe 1 --arjunsimplev 2 --arjunbackwmaxc 10000 --
                    # HOWEVER, duplication caused cadiback to go bad
            # both below are nice, with and without appmc
            # "out-ganak-7433320.pbs101-0", # default config
            # "out-ganak-7433320.pbs101-1", # --appmct 2000
            #SHITTY"out-ganak-7435410.pbs101-"# try all tdexp, go back to before mess
            # "out-ganak-7451423.pbs101-1", # gaank 32GB, --maxcache 26000 BEST. Max mem useage 30GB
            # "out-ganak-7466346.pbs101-2", # d4 32GB mem
            # "out-ganak-7466346.pbs101-", # diffocc vs non-diffocc
            # "out-ganak-7468615.pbs101-0", # without vivif
            # "out-ganak-7468615.pbs101-1", # with vivifevery 40k, BEST: --maxcache 26000 --vivifevery 40000
             # c o CMS revision: e897d58afa4caa767d9b4e6fe9eab244ffcdb97b
             # c o Arjun SHA revision: f8dfd3b824b2c4b404016fa7ed7fad0460dfc7ae
             # c o Arjun SBVA SHA revision: 4736d909ce8fe17b72429293d8285606a8e26925
             # c o GANAK SHA revision bf5867815d357f30911c0011b9d47a5941119107
            # "out-ganak-7549763.pbs101" # new run, with/without clause activities -- slightly slower. clact does nothing. GOOD, without clact


            # weighted wmc, i.e. track 2
            #"out-ganak-7484977.pbs101-", # 32GB, first run, segfault, 256B float
            # "out-ganak-7491797.pbs101-1", # ran out of memory float
            # "out-ganak-7491797.pbs101-2", # d4
            # "out-ganak-7491797.pbs101-3", # gpmc
            # "out-ganak-7503909.pbs101-0", # 40k vivif is slower
            # "out-ganak-7491797.pbs101-0", # no-memout, no-segfault mpq -- BEST: --maxcache 18000 --precise 1
            # "out-ganak-7503909.pbs101-1", # no-memout, no-segfault mpf

            # weighted, projected i.e. track 4
            # "out-ganak-7505064.pbs101-0", # BEST ganak: --maxcache 18000 --tdminw 0.05 --tdexpmult 0.3 --precise 1
            # "out-ganak-7505064.pbs101-4", # 32GB gpmc
            # "out-ganak-7505064.pbs101-3", # 32GB d4


            # unweighted, projected, i.e. track 3
            # "out-ganak-7173534.pbs101-0", # proj-2023 first run, 4GB
            # "out-ganak-7180435.pbs101-0", # proj-2023 4B, vsads
            # "out-ganak-7184202.pbs101-" #proj 2023, 4GB, freq-128 (vsads), checking vsads params
            # "out-ganak-7189867.pbs101-1" # proj 2023, 4GB, cache elem simplification, TD back to 0.1, SAT restarts
            # "out-ganak7197939",
            # "out-ganak-7206369.pbs101-",
            # TODO: higher maxw than out-ganak-7205692.pbs101-6, rerun fixed cadiback
            # "out-ganak-7178163.pbs101-0", # proj-2023 16GB, ganak
            # "out-ganak-7178163.pbs101-2", # proj-2023 16GB, d4
            # "out-ganak-7178163.pbs101-3", # proj-2023 16GB, gpmc
            # "out-ganak-7205692.pbs101-6", # old best
            # #"out-ganak-7247003.pbs101-", # fixed cadiback, --tdmaxw 40 seems best, but 15 is fine too actually # BEST
            # "out-ganak-7255014.pbs101-0", # different tdminw-s, also try lbd 2. LBD is indifferent, tdminw high is BAD
            # "out-ganak-7266814.pbs101-1", # contract over TDW, higher tdmaxw, higher tdminw. Was buggy in a few ways. Let's re-run. Fixed contraction in the meanwhile.
            # "out-ganak-7294423.pbs101-0", # fixed memout from contraction, fixed too much extend, distill-bin, BIG before backbone
            # "out-ganak-7307327.pbs101-4", # fixed cadical, fixed resolv-subs, fixed oracle, 2 new options for extra oracle & resolv-subs
            # "out-ganak-7308235.pbs101-", # newest cadical -- ALL very good
                                            # best is "--tdmaxw 10 --varfreqdiv 30"
            # BAD PARAMS.... "out-ganak-7316071.pbs101-", # also appmc, play around with some options
            # "out-ganak-7320764.pbs101-0", #best: --tdmaxw 10 --varfreqdiv 25 , --appmct 2000
            # "out-ganak-7320764.pbs101-9", #best: --tdmaxw 10 --varfreqdiv 25 , --appmct 2000
            # "out-ganak-7320764.pbs101-0", #best: --tdmaxw 10 --varfreqdiv 25 , --appmct 2000 (could try varfreqdiv 20)
            # "out-ganak-7334726.pbs101-1", #  new run, some parameter tuning, fixed SAT restart bug (again, this time good) -- BEST, non-appmc, --tdmaxw 10 --varfreqdiv 25
            # "out-ganak-7334726.pbs101-8", #  new run, some parameter tuning, fixed SAT restart bug (again, this time good) -- BEST, appmc
            #"out-ganak-7348395.pbs101-", # better extend -- maybe extend is not good?
            # todo: no sbva, no extend, --compsort 5 & 4,
            # "out-ganak-7377549.pbs101", # release version, try different bacdkbone setups (probe, only oracle vivif, etc), tdminw 0
            # "out-ganak-7366311.pbs101-4", # --backbone 0 is good, --tdexpmult 0.3 + --tdmaxw 100 is good
            # "out-ganak-7366311.pbs101-6", # --backbone 0 is good, --tdexpmult 0.3 + --tdmaxw 100 is good
            # "out-ganak-7377549.pbs101-", # best: --tdmaxw 100 --varfreqdiv 25 --tdexpmult 0.3 --backbone 0 --arjunprobe 1
            # SHITTY "out-ganak-7400920.pbs101-0", # timeout in backbone 30s, went really bad
            # "out-ganak-7419164.pbs101-3", # more simplification pre-backward, lower backward confl
                # BEST: --tdexpmult 0.3 --backbone 0 --arjunprobe 1 --arjunsimplev 2 --arjunbackwmaxc 20000
            # SHITTY "out-ganak-7422213.pbs101-" # unit cls to duplicated CNF, different cutoff for gates, some param tuning Wow, wrong.
                # Seemingly to do with cadiback what is going on....
            # "out-ganak-7435405.pbs101-", # try all tdexp, go back to before mess. Gates all turned off.
            # "out-ganak-7433580.pbs101-2", #d4 32G
            # "out-ganak-7511977.pbs101-0", # gpmc 32G
            # best: --tdminw 0.05 --tdexpmult 0.3
            # "out-ganak-7435405.pbs101-", # try all tdexp, go back to before mess. Gates all turned off.
            # best: --tdminw 0.05 --tdexpmult 0.3
            # "out-ganak-7451414.pbs101-1", # ganak 32G mem, --maxcache 26000 --tdminw 0.05 --tdexpmult 0.3 BEST. Max mem usage 29.82GB
            # "out-ganak-7466418.pbs101-0", # diffocc -- but still slower
            # "out-ganak-7466418.pbs101-1", # good, but with clause/variable early-exit, it's even faster below
            # "out-ganak-7468556.pbs101-1", # (32GB ganak) BEST: --maxcache 26000 --tdminw 0.05 --tdexpmult 0.3
             # c o CMS revision: e897d58afa4caa767d9b4e6fe9eab244ffcdb97b
             # c o Arjun SHA revision: f8dfd3b824b2c4b404016fa7ed7fad0460dfc7ae
             # c o Arjun SBVA SHA revision: 4736d909ce8fe17b72429293d8285606a8e26925
             # c o GANAK SHA revision bf5867815d357f30911c0011b9d47a5941119107
            # --> above also proves --arjunextend 0 is BAD.
            #"out-ganak-7482756.pbs101-", # now with more vivif -- bad, default is OK
            # BAD version "out-ganak-7515096.pbs101-" # try TD once more... BAAAD due to clause keeping
            # BAD version "out-ganak-7513732.pbs101-", # total use BAAAAD
            # "out-ganak-7513092.pbs101-", # different viviv cutoff totaluse
            # "out-ganak-7531643.pbs101-" # contract and new prime graph and TD
            # "out-ganak-7549745.pbs101-2", # clause activities/no cl activities (plus minor changes)
                                        # GOOD, without clact
            # "out-ganak-7559399.pbs101-1" # now with/without backbone only over opt-indep-set


            # arijit's experiment
            # "out-ganak-7511833.pbs101-", # pmc -> mc and run with all systems
            # "out-ganak-7433580.pbs101-2", #d4 32G     # pmc
            # "out-ganak-7511977.pbs101-0", # gpmc 32G  # pmc
            # "out-ganak-7468556.pbs101-1", # ganak 32G  # pmc

            ######################
            # 2024 track 1 (mc) public instances
            ######################
            # "out-ganak-7558933.pbs101-", # all very good
            # forget this: good:  --maxcache=26000 --maxcache=26000 --vivifevery 40000
            # BEST: --maxcache=26000 --tdminw 8 --tdmaxw 70
            # had 2GB left over.
            # "out-ganak-7602594.pbs101-", # with arjunoraclefindbins
            # "out-ganak-7601364.pbs101-", ## arjunoraclefindbins (ran it twice, oops)
            # "out-ganak-7608982.pbs101-", # more arjunoraclefindbins
            # "out-ganak-7608400.pbs101-", # more arjunoraclefindbins (or maybe sharpsat?)
            # "out-ganak-7608400.pbs101-5", # sstd
            # BEST: --maxcache=26000 --tdminw 8 --tdmaxw 70 --arjunoraclefindbins 4
            # TODO: WOW, sstd is doing very well. Check out: out-ganak-7608982.pbs101-0
            # "out-ganak-7622771.pbs101-", # changed clause deletion
            # "out-ganak-7622834.pbs101-", # changed clause deletion, update LBD
            # "out-ganak-7623566.pbs101-", # new cldel, sbva, etc
            # "out-ganak-7623566.pbs101-4", # BEST
            # BEST: call: --maxcache=24000 --sbva 1 --tdminw 5 --tdmaxw 50 --arjunoraclefindbins 4 --rdbclstarget 10000 # 2.8GB left over
            # TODO larger tdmin/max, different tdexp -- running: out-ganak-7625957.pbs101-0
            # "out-ganak-7625957.pbs101-", # new sbva configs
            # "out-ganak-7625957.pbs101-5", # new sbva configs (best)
            # BEST:
            # dir: out-ganak-7625957.pbs101-5   # 1GB left only!!
            # call: --maxcache=23000 --arjunverb 2 --sbva 1 --tdexpmult 1.1 --tdminw 7 --tdmaxw 60 --arjunoraclefindbins 6 --rdbclstarget 10000
            # "out-ganak-7628982.pbs101-", # getting bins from cadiback, more bin finding
            # "out-ganak-7631793.pbs101-", # arjun total timeout checks (4GB mem)
            # "out-ganak-7632152.pbs101-", # go back to old, restart, options, slightly more bins (4GB mem)
            # "out-ganak-7632478.pbs101-", # back to 32GB, restarts, go back to old as before
            # "out-ganak-7637537.pbs101-", # get bins from cadiback, do some bin pruning, only get <= 3 lbd learnt clauses
            # DAMN!! the issue was that --sbva 1 effectively TURNSOFF SBVA!!!
            # "out-ganak-7637617.pbs101-", # like above, but --sbva 1. Kinda OK.
            # "out-ganak-7637716.pbs101-", #  sbva configs BAD, --sbva 1 is likely best
            # "out-ganak-7639585.pbs101-", # different arjun setup
            # "out-ganak-7643061.pbs101-", # double-arjun in case it's small, no bve expand on extra BVE
            # "out-ganak-7644331.pbs101-", # FINAL
            # "out-ganak-7657683.pbs101-", # FINAL, fixed counting, basically the same
            # "out-ganak-7664145.pbs101-", # FINAL, different bve size limits -- no real advantage
            # "out-ganak-7664153.pbs101-", # FINAL, different bve size limits -- accidental re-run, no real advantage
            # "out-ganak-7665630.pbs101-", # FINAL, bvemaxres 10-12-14
              # nice: --bveresolvmaxsz 12
              # --bveresolvmaxsz 12 --maxcache=22000 --sbva 1 --tdexpmult 1.1 --tdminw 7 --tdmaxw 60 --arjunoraclefindbins 6 --rdbclstarget 10000
            # "out-ganak-7669369.pbs101-", # FINAL, double CNF minim -- NICE -- 1.2GB remain



            ######################
            # 2024 track 2 (wmc) public instances
            ######################
            # "out-ganak-7559210.pbs101-", # all are very good
            # BEST: --maxcache=18000
            # had 2GB left over.
            # "out-ganak-7601766.pbs101-", # with arjunoraclefindbins
            # BEST: --maxcache=18000 --arjunoraclefindbins 1
            # "out-ganak-7608600.pbs101", # more arjunoraclefindbins
            # BEST: --maxcache=18000 --arjunoraclefindbins 8
            # NOTE: sstd is very good, but it's running unweighted: out-ganak-7608600.pbs101-5
            # "out-ganak-7608600.pbs101-5", # sstd unweighted
            # "out-ganak-7623575.pbs101-", # sbva, tdstuff, etc.
            # "out-ganak-7623575.pbs101-0", # sbva, tdstuff, etc.
            # BEST:
            # dir: out-ganak-7623575.pbs101-0 # 3.5GB left over
            # call: --maxcache=16000 --sbva 1 --tdminw 18 --arjunoraclefindbins 4 --rdbclstarget 10000
            # "out-ganak-7626133.pbs101-", # new sbva configs
            # "out-ganak-7626133.pbs101-6", # new sbva configs
            # BEST:
            # dir: out-ganak-7626133.pbs101-6 # 3.3GB left over
            # call: --maxcache=16000 --arjunverb 2 --sbva 1 --sbvalitcut 6 --tdminw 5 --tdmaxw 50 --arjunoraclefindbins 4 --rdbclstarget 10000
            # "out-ganak-7643075.pbs101-", # double-arjun in case it's small, no bve expand on extra BVE -- BAD
            # "out-ganak-7644361.pbs101-", # FINAL
            # "out-ganak-7657686.pbs101-", # FINAL, fixed counting, basically the same
            # "out-ganak-7664158.pbs101-",  # FINAL, different bve size limits
            # "out-ganak-7665603.pbs101-", # FINAL, bvemaxsz 10-12-14
              # nice: --bveresolvmaxsz 12
              # --bveresolvmaxsz 12 --maxcache=16000 --sbva 1 --sbvalitcut 6 --tdminw 5 --tdmaxw 50 --arjunoraclefindbins 4 --rdbclstarget 10000
            # "out-ganak-7669419.pbs101-", # FINAL, doing double CNF minim -- NICE -- 3.2 GB remain



            ######################
            # 2024 track 3 (i.e. pmc) public instances
            ######################
            # "out-ganak-7559160.pbs101-",
            # # BEST: --maxcache=26000 -- YES, no minw/maxw!
            # had 2GB left over.
            # "out-ganak-7601744.pbs101-", # arjunoraclefindbins
            # BEST: --maxcache=26000 --arjunoraclefindbins 1
            # "out-ganak-7608684.pbs101-",# more arjunoraclefindbins
            # BEST: --maxcache=26000 --tdminw 8 --tdmaxw 70 --arjunoraclefindbins 2
            # "out-ganak-7623571.pbs101-", # sbva, tdstuff, etc
            # BEST seems:
            # "out-ganak-7623571.pbs101-5", # 2GB left over
            # call: --maxcache=24000 --arjunverb 2 --sbva 1 --tdminw 5 --tdmaxw 50 --arjunoraclefindbins 4 --rdbclstarget 10000
            # "out-ganak-7626141.pbs101-", # new sbva configs -- not too good
            # "out-ganak-7629018.pbs101-", # getting bins from cadiback, more bin finding
            # "out-ganak-7631789.pbs101-", # arjun total timeout checks
            # "out-ganak-7637633.pbs101-", # restarts -- bad (but may be because of mem/cpu)
            # "out-ganak-7639688.pbs101-", # different arjun setup
            # "out-ganak-7643068.pbs101-", # double-arjun in case it's small, no bve expand on extra BVE -- BAD
            # "out-ganak-7644550.pbs101-0", # FINAL
            # "out-ganak-7657684.pbs101-", # FINAL, fixed counting, not better.
            # "out-ganak-7665625.pbs101-", # FINAL, bvemaxsz 10-12-14
              # --bveresolvmaxsz 12 is NICE
              # --bveresolvmaxsz 12 --maxcache=24000 --sbva 1 --tdminw 5 --tdmaxw 50 --arjunoraclefindbins 4 --rdbclstarget 10000
            # "out-ganak-7670361.pbs101-", # FINAL, doing double CNF minim -- NICE. 2GB remain


            ######################
            # 2024 track 4 (i.e. pwmc) public instances
            ######################
            # "out-ganak-7601777.pbs101-",
            # BEST: --maxcache=18000 --arjunverb 2
            # had 7GB left over (with --precise 0, 8.5GB left over)
            # "out-ganak-7608937.pbs101", # more arjunoraclefindbins
            # BEST: --maxcache=18000 --arjunoraclefindbins 2
            # "out-ganak-7623679.pbs101-", # sbva, new cldel, etc.
            # "out-ganak-7623679.pbs101-8", # BEST
            # BEST: --maxcache=16000 --arjunverb 2 --sbva 0 --arjunoraclefindbins 4 --rdbclstarget 14000 # 12 GB left over
            # "out-ganak-7623679.pbs101-2", #good
            # GOOD: --maxcache=16000 --arjunverb 2 --sbva 1 --tdminw 18 --arjunoraclefindbins 4 --rdbclstarget 14000 # had 7 GB left over
            # "out-ganak-7626015.pbs101-", # new sbva configs
            # "out-ganak-7626015.pbs101-6", # new sbva configs
            # BEST:
            # dir: out-ganak-7626015.pbs101-6 # 4GB left over
            # call: --maxcache=20000 --arjunverb 2 --sbva 1 --tdminw 18 --sbvalitcut 6 --arjunoraclefindbins 4 --rdbclstarget 14000
            # "out-ganak-7643072.pbs101-", # double-arjun in case it's small, no bve expand on extra BVE
            # "out-ganak-7644569.pbs101-", # FINAL
            # "out-ganak-7657681.pbs101-", # FINAL, fixed counting, basically the same
            # "out-ganak-7664168.pbs101-", # FINAL different bve size limits
            # "out-ganak-7665622.pbs101-", #FINAL, bvemaxsz 10-12-14
               # nice: --bveresolvmaxsz 12
               # --bveresolvmaxsz 12 --maxcache=20000 --sbva 1 --tdminw 18 --sbvalitcut 6 --arjunoraclefindbins 4 --rdbclstarget 14000
            # "out-ganak-7670422.pbs101-", #-- FINAL, doing double CNF minim -- OK, 7GB remain

            ######################
            # paper -- Fixed job
            # "out-ganak-mc2324-13889246-0",  # all, 2023+24, major combos
            # "out-ganak-mc2324-13889246-1/", # all, 2023+24, major combos
            # "out-ganak-mc2324-13889246-2",  # all, 2023+24, major combos
            # "out-ganak-mc2324-13889246-3",  # all, 2023+24, major combos
            # "out-ganak-mc2324-13889246-4",  # all, 2023+24, major combos

            # update for chronoBT, different clause deletion
            # "out-ganak-mc2324-13890807-"

            # all, SAT combos
            # "out-ganak-mc2324-13889246-4",  # all, 2023+24, major combos
            # "out-ganak-mc2324-13889246-6",  # all, 2023+24, SAT combos
            # "out-ganak-mc2324-13889246-7",  # all, 2023+24, SAT combos
            # "out-ganak-mc2324-13889246-8",  # all, 2023+24, SAT combos
            # "out-ganak-mc2324-13889246-9",  # all, 2023+24, SAT combos

            # allindep
            # "out-ganak-mc2324-13889246-4",  # all, 2023+24, major combos
            # "out-ganak-mc2324-13889246-5",  # all, 2023+24, allindep

            # tuning
            # "out-ganak-mc2324-13889246-4",  # all, 2023+24, major combos
            # "out-ganak-mc2324-13889246-11",  # all, 2023+24
            # "out-ganak-mc2324-13889246-12",  # all, 2023+24

            # bdd and mini-ablation study
            # "out-ganak-mc2324-13911802-",  # bdd and mini-ablation study

            # holder+dirty new idea
            # "out-ganak-mc2324-13929441-", # dirty
            # "out-ganak-mc2324-13930324-", # nodirty
            # "out-ganak-mc2324-13931851-", # better nodirty
            # "out-ganak-mc2324-13935870-", # better nodirty, double bump
            # "out-ganak-mc2324-13943683-", # as above, but bin cl early abort, and mem layout for Comp -- BUGGY (used to be "GOOD")
            # "out-ganak-mc2324-13944875-", # as above, but bin cl (and long cl partially) tighter, cleanup
            # "out-ganak-mc2324-13945685-", # something went wrong above, so I went back 2x to "GOOD" but with a fix for bin that was maybe buggy
            # "out-ganak-mc2324-13946763-", # finally fixed binary counting
            # "out-ganak-mc2324-13972052-", # tighter datastructs via vec<>, float instead of double, double bump bin, fixed counting of removed literals for further clause minimization
            # "out-ganak-mc2324-13973015-", # double instead of float for activities
            # "out-ganak-mc2324-13977963-", # all experiments
            # "out-ganak-mc2324-13982113-", # new stamping
            # "out-ganak-mc2324-13982787-", # fixed score on satisfied clauses, blocked literal update
            # "out-ganak-mc2324-13994954-", # checking freq score divider (also move-to-front Cldata + 7 elems, no 8)

            # everything above is new2 now.
            # "out-ganak-mc2324-14063135-",
            # "out-ganak-mc2324-14065222-4",
            # "out-ganak-mc2324-14065222-0",

            # SAT solver combos
            # "out-ganak-mc2324-14065222-1",
            # "out-ganak-mc2324-14065222-2",
            # "out-ganak-mc2324-14065222-3",
            # "out-ganak-mc2324-14065222-4",

            # other solvers
            #"out-others-14064944", # sharpsat

            # "out-ganak-mc2324-14063135-1/", # baseline
            # "out-ganak-mc2324-14072353-1/", # basic SAT + ChronoBT
            # "out-ganak-mc2324-14063135-7/", # also enchanced SAT
            # "out-ganak-mc2324-14063135-4/", # also dual indep
            "out-ganak-mc2324-14063135-0/", # also extend decision (all in)
            # "out-others-gpmc",
            # "out-others-d4",
            # "out-others-14064944-0", # sharpsat-td

            # unused
            # "out-ganak-mc2324-14072353-0/", # all in, but TD contract off
            # "out-ganak-mc2324-14086287-" # sbva setup
            # "out-ganak-mc2324-14098333-", # sbva checks

            # checking revert, and cache del idea, and TD max limit
            # "out-ganak-mc2324-14145161-0", # revert, cache del, TD max limit

            # new ccnr ideas, oracle-sparsify timeout different, 0 indep set no puura
            # "out-ganak-mc2324-14164796-",

            # fixed bug in oracle-spars, more options
            # "out-ganak-mc2324-14165187-",

            # more vivification during minimization (and with it more binary XORs), apparently this matters a lot
            # "out-ganak-mc2324-14167470-",
            # c o GANAK SHA revision e246dcbba97bc937806f0863738478eec6c5c760
            # c o CMS revision: c15cd211779400ba8ed3ec45c24666487f0acec1
            # c o Arjun SHA revision: f24048a43411cb52410526b80eff2720d32bc942
            # c o Arjun SBVA SHA revision: 5912435affe8c77ecf364093cea29e0fc5c1b5cb

            # even more minimization during vivification, etc
            # "out-ganak-mc2324-14171543-",

            # fixed bugs related to ccnr backbone, and using purely CMS for very small-indep size (and no puura for <10)
            # this WAS the best
            # "out-ganak-mc2324-14195954-0",

            # complex numbers from Mei
            # "out-ganak-gpmc-mei-14294765-", # gpmc
            # "out-ganak-gpmc-mei-14295250-0", # ganak, fixed complex
            # "out-ganak-gpmc-mei-14295899-0", # mpfr ganak

            # NOTE! likely WRONG because it can overflow the timestamp
            # "out-ganak-mccomp2324-14366334-0", # update archetype only when sup comp has at least 20 vars that are still decidable

            # blocking
            # "out-ganak-mccomp2324-14299825",
            # "out-ganak-mccomp2324-14309534-0", # this IS the current best
            # "out-ganak-mccomp2324-14309534-1", # messed up
            # "out-ganak-mccomp2324-14321462-" # let's try blocking without messing it up
            # "out-ganak-mccomp2324-14362931-0", # new hash, timestamp archetype
            # "out-ganak-mccomp2324-14391524-", # changeable recompute archetype cutoff, better cache usage
            # "out-ganak-mccomp2324-14513393-0", # update arjun's order and add gate-based-eqlit, also fix gate-based-rem by not doing strengthening
            # "out-ganak-mccomp2324-14513746-0", # back to previous, except the strengthening fix for gate-based-rem
            # "out-ganak-mccomp2324-14519572-", # as previous 2, but fixed parsing of stupid weights
            # "out-ganak-mccomp2324-14528572-", # binaries in holder are no longer updated
            # "out-ganak-mccomp2324-14538562-1", # full-probe before everything in remove_eq_lits
            # "out-ganak-mccomp2324-14549970-1", # LTO
            # "out-ganak-mccomp2324-14552736-", # no LTO, but some change (e.g. sandybridge etc)
            # "out-ganak-mccomp2324-14558670-1", # some stuff enabled previously for no flto
            # "out-ganak-mccomp2324-14575065-", # occ-bve-resolv after bve
            # "out-ganak-mccomp2324-14576195-", # as 2 above, but now playing with TD iters
            # "out-ganak-mccomp2324-14621165", # float, also try different TD exp
            # "out-ganak-mccomp2324-14621994-0", # float fixed mem free, try different tdmaxw
            # "out-ganak-mccomp2324-14624577", # exact, no cache check, exp check
            # "out-ganak-mccomp2324-14624580", # exact, more exp check
            # "out-ganak-mccomp2324-14637456", # try vivif options
            # "out-ganak-mccomp2324-14648650-0", # much more precise cache
            # "out-ganak-mccomp2324-14655977-", # get bins from cadiback
            # "out-ganak-mccomp2324-14656347-0", # get bins, scc, enable weakening, more weakening -- GOOD
            # "out-ganak-mccomp2324-14661402-", # try different weakenings
            # "out-ganak-mccomp2324-14675861-0", # submitted to mccomp 2025 -- GOOD
            # "out-ganak-mccomp2324-14675861-", # default setup along WITHOUT appmc. Trying tditers -- ganak_7d97636055e_9104724fa_26d64aac --tdlookiters 20
            # "out-ganak-mccomp2324-14675861-1", # default setup along WITH appmc. Trying tditers
            # "out-ganak-mccomp2324-15010600-0", # TD start from 0
            # "out-ganak-mccomp2324-3382-0", # TRILLIUM -- TD start from 0
            # "out-ganak-mccomp2324-21238", # TRILLIUM --non-eq, and eq, and no probabilistic
            # "out-ganak-mccomp2324-21349-0", # TRILLIUM, ganak_7d97636055e_9104724fa_26d64aac (i.e. old run that was the fastest)
            # "out-ganak-mccomp2324-21481-0", # TRILLIUM, fixing memory usage for --prob 0
            # "out-ganak-mccomp2324-33393-0", # TRILLIUM, fixing td starting from 0, and maybe other TD issues too
            # "out-ganak-mccomp2324-35541-0", # TRILLIUM, fixing td starting from 0, now cutting disjoint components at toplevel for correct centroid
            # "out-ganak-mccomp2324-37704-0", # TRILLIUM, trying less scc upfront, weird non-committed code
            # "out-ganak-mccomp2324-41122-0", # TRILLIUM, new TD score, but BAD, not updating depth when all vars have been decided
            # "out-ganak-mccomp2324-41923-0", # TRILLIUM, new TD score, fixed depth update finally
            # "out-ganak-mccomp2324-43391-0", # TRILLIUM, new TD score, fixed depth update
            # "out-ganak-mccomp2324-44137-0", # TRILLIUM, finally good fix for TD
            ### "out-ganak-mccomp2324-209249", # TRILLIUM, with 64b fixes
            "out-ganak-mccomp2324-343419-0", # hypercut tries
            "out-ganak-mccomp2324-357608-0", # hypercut only when TD is not used
            ]
# only_dirs = ["out-ganak-6828273"]
# only_dirs = ["6606250"]
# not_calls = ["--nvarscutoffcache 20", "--nvarscutoffcache 3"]
# not_calls = ["--satsolver 0"]
not_versions = []
# only_dirs = []
# only_calls = ["vsadsadjust 128"]
only_calls = []
# not_calls = ["restart"]
not_calls = []
todo = versions

# all
fname_like = ""
# unproj only
# fname_like = " and (fname like '%track1%' or fname like '%track2%') "
# proj only
# fname_like = " and (fname like '%track3%' or fname like '%track4%') "

table_todo = []
for ver in todo :
    dirs_call = get_dirs(ver)
    for dir,call in dirs_call:
        bad = False
        for not_call in not_calls:
          if not_call in call:
            bad = True
        for not_version in not_versions:
          if not_version in ver:
            bad = True

        if len(only_calls) != 0:
          inside = False
          for only_call in only_calls:
            if only_call in call:
              inside = True
          if not inside: bad = True

        if len(only_dirs) != 0:
          inside = False
          for only_dir in only_dirs:
            if only_dir in (dir+"/"):
              inside = True
          if not inside: bad = True

        if bad:
          continue
        # print("----")
        # print("dir:", dir)
        # print("call:", call)
        # print("ver:", ver)

        # if "actexp 1.0" in call:
        #     continue
        # if "tdwithredbins 0" in call:
        #     continue
        # if "restart 1" in call:
        #     continue
        # if "vivif" not in call:
        #     continue
        # if "probe 1" in call:
        #     continue
        # if "polar" not in call:
        #     continue
        fname = "run-"+dir+".csv"
        with open("gencsv.sqlite", "w") as f:
            f.write(".headers off\n")
            f.write(".mode csv\n");
            f.write(".output "+fname+"\n")
            extra = ""
            f.write("select ganak_time from data where dirname='"+dir+"' and ganak_ver='"+ver+"'\n and ganak_time is not NULL "+fname_like)
        os.system("sqlite3 mydb.sql < gencsv.sqlite")
        os.unlink("gencsv.sqlite")

        fname2 = fname + ".gnuplotdata"
        num_solved = convert_to_cactus(fname, fname2)
        fname2_s.append([fname2, call, ver[:10], num_solved, dir])
        table_todo.append([dir, ver])

for only_counted in [False, True]:
  counted_req = ""
  if only_counted:
    print("::: --------- Data based on ONLY benchmarks that are COUNTED ------- :::")
    counted_req = " and ganak_time is not NULL "
  else:
    print("::: --------- Data based on ALSO UNCOUNTED benchmarks ------- :::")
  with open("gen_table.sqlite", "w") as f:
    f.write(".mode table\n")
    # f.write(".mode colum\n")
    # f.write(".headers off\n")
    dirs = ""
    vers = ""
    for dir,ver in table_todo:
      dirs += "'" + dir + "',"
      vers += "'" + ver + "',"
    dirs = dirs[:-1]
    vers = vers[:-1]
    extra = ""
    f.write("select \
        replace(dirname,'out-ganak-mc','') as dirname,\
        replace(ganak_call,'././ganak_','') as call,\
        sum(mem_out) as 'mem out', \
        sum(signal == 11) as 'sigSEGV', \
        sum(signal == 6) as 'sigABRT', \
        sum(signal == 14) as 'sigALRM', \
        sum(signal == 8) as 'sigFPE', \
        CAST(ROUND(avg(ganak_mem_MB), 0) AS INTEGER) as 'av memMB',\
        ROUND(avg(conflicts)/(1000.0*1000.0), 2) as 'av confM', \
        ROUND(avg(decisionsK)/(1000.0), 2) as 'av decM', \
        CAST(ROUND(avg(sat_called/1000.0),0) AS INTEGER) as 'av satcK',\
        sum(ganak_time is not null) as 'solved',\
        CAST(ROUND(sum(coalesce(ganak_time, 3600))/COUNT(*),0) AS INTEGER) as 'PAR2',\
        CAST(avg(opt_indep_sz-indep_sz) AS INTEGER) as 'avg-diff-opt-sz',\
        ROUND(avg(gates_extend_t), 3) as 'gates-ext-t',\
        ROUND(avg(padoa_extend_t), 3) as 'padoa-ext-t',\
        ROUND(avg(gates_extended), 3) as 'gates-ext',\
        ROUND(avg(padoa_extended), 3) as 'padoa-ext',\
        CAST(ROUND(max(cache_del_time), 0) AS INTEGER) as 'max cachdT',\
        CAST(ROUND(avg(backbone_time),0) AS INTEGER) as 'av backT',\
        CAST(ROUND(avg(arjun_time),0) AS INTEGER) as 'av arjT',\
        CAST(ROUND(avg(td_time),0) AS INTEGER) as 'av tdT',\
        ROUND(avg(td_width),0) as 'av tdw',\
        ROUND(avg(cache_miss_rate),2) as 'av cmiss',\
        ROUND(avg(cache_avg_hit_vars),2) as 'av chitvs',\
        ROUND(avg(ganak_mem_mb),2) as 'av memMB',\
        ROUND(avg(compsK/1000.0),2) as 'av compsM',\
        sum(fname is not null) as 'nfiles'\
        from data where dirname IN ("+dirs+") and ganak_ver IN ("+vers+") "+fname_like+" "+counted_req+"group by dirname order by solved asc")
  os.system("sqlite3 mydb.sql < gen_table.sqlite")
  os.unlink("gen_table.sqlite")

# out-ganak-mc2324-14063135-0|mc2023_track4_050.cnf|././ganak_d622244b9c9 --arjunverb 2 --maxcache 5000|4294967295
if True:
  for dir,ver in table_todo:
    with open("gen_table.sqlite", "w") as f:
      f.write(".mode table\n")
      f.write("select '"+dir+"', '"+ver+"'");
      for col in "indep_sz", "opt_indep_sz", "orig_proj_sz", "new_nvars", "ganak_mem_mb":
        f.write(", (SELECT "+col+" as 'median_"+col+"'\
        FROM data\
        where dirname IN ('"+dir+"') and ganak_ver IN ('"+ver+"') and "+col+" is not null"+fname_like+"\
        ORDER BY "+col+"\
        LIMIT 1\
        OFFSET (SELECT COUNT("+col+") FROM data\
          where dirname IN ('"+dir+"') and ganak_ver IN ('"+ver+"') \
          and "+col+" is not null) / 2) as median_"+col+" \
      ")
      for col in "gates_extended", "padoa_extended":
        f.write(", (SELECT "+col+" as 'median_"+col+"_NOZERO'\
        FROM data\
        where dirname IN ('"+dir+"') and ganak_ver IN ('"+ver+"') and "+col+" is not null "+fname_like+"\
                    and "+col+">0\
        ORDER BY "+col+"\
        LIMIT 1\
        OFFSET (SELECT COUNT("+col+") FROM data\
          where dirname IN ('"+dir+"') and ganak_ver IN ('"+ver+"') \
          and "+col+" is not null "+fname_like+" and "+col+">0) / 2) as median_"+col+" \
      ")
    os.system("sqlite3 mydb.sql < gen_table.sqlite")
    os.unlink("gen_table.sqlite")

gnuplotfn = "run-all.gnuplot"
with open(gnuplotfn, "w") as f:
    f.write("set term postscript eps color lw 1 \"Helvetica\" 12 size 9,5\n")
    f.write("set output \"run.eps\"\n")
    f.write("set title \"Counter ganak\"\n")
    f.write("set notitle\n")
    f.write("set key bottom right\n")
    # f.write("set xtics 200\n")
    f.write("unset logscale x\n")
    f.write("unset logscale y\n")
    f.write("set ylabel  \"Instances counted\"\n")
    f.write("set xlabel \"Time (s)\"\n")
    # f.write("plot [:][10:]\\\n")
    f.write("plot [:][800:]\\\n")
    i = 0
    # f.write(" \"runkcbox-prearjun.csv.gnuplotdata\" u 2:1 with linespoints  title \"KCBox\",\\\n")
    # f.write(" \"runsharptd-prearjun.csv.gnuplotdata\" u 2:1 with linespoints  title \"SharptTD\",\\\n")
    towrite = ""
    for fn,call,ver,num_solved,dir in fname2_s:
        # if "restart" not in call and num_solved > 142:
        if True:
            call = gnuplot_name_cleanup(call)
            dir  = gnuplot_name_cleanup(dir)
            ver  = gnuplot_name_cleanup(ver)
            oneline = "\""+fn+"\" u 2:1 with linespoints  title \""+ver+"-"+dir+"-"+call+"\""
            towrite += oneline
            towrite +=",\\\n"
    towrite = towrite[:(len(towrite)-4)]
    f.write(towrite)


if os.path.exists("run.eps"):
  os.unlink("run.eps")
if os.path.exists("run.pdf"):
  os.unlink("run.pdf")
if os.path.exists("run.png"):
  os.unlink("run.png")

def create_notebook():
  # Create a new notebook
  nb = nbf.v4.new_notebook()
  texts = []

  # Define the text content for the markdown cells
  text = """
# Step 1: Import necessary libraries
import pandas as pd
import sqlite3
import matplotlib.pyplot as plt
from functools import reduce
import numpy as np

dirs = ["""+dirs+"""]
"""
  texts.append(text)

  colnames = ["ganak_time", "ganak_mem_MB", "conflicts", "decisionsK", \
      "compsK", "td_width", "arjun_time", "backbone_time", "indep_sz", "td_time", "sat_called", "cache_miss_rate"]
  for colname in colnames:
    text= """
colname='"""+colname+"""'
names=[]
dfs = []
conn = sqlite3.connect('mydb.sql')
for d in dirs:
  # Step 3: Run the SQL query and load the results into a DataFrame
  query = "SELECT fname, "+colname+", ganak_call FROM data where "+colname+" is not NULL and dirname='"+d+"' order by "+colname
  df1 = pd.read_sql_query(query, conn)
  df1['num'] = range(len(df1))
  dfs.append(df1)
  names.append(d+" " +df1['ganak_call'][0])

for i in range(len(dfs)):
    for c in dfs[i].columns:
      if c == 'num': continue
      dfs[i].rename(columns={c: c+'_'+str(i)}, inplace=True)
result = dfs[0]

for i, df in enumerate(dfs[1:], start=1):
    result = pd.merge(result, df, on='num', how='outer', suffixes=(f'_{i-1}', f'_{i}'))

# Merge all DataFrames
conn.close()

# Step 5: Plot the data
plt.figure(figsize=(10, 6))
for i in range(len(dirs)):
  values = result[(colname + '_' + str(i))]
  min_non_zero = values[values > 0].min()
  adjusted_values = np.where(values == 0, min_non_zero, values)
  log_values = np.log10(adjusted_values)
  plt.plot(result['num'],log_values,marker='o')
plt.title('Plot of CNFs counted vs. log10 '+colname)
plt.xlabel('Number of CNFs counted')
plt.ylabel('log10 '+colname)
plt.legend(names,loc='center left', bbox_to_anchor=(0, -0.3))
plt.grid(True)
plt.show()
  """
    texts.append(text)

  for col1, col2 in [("td_width", "ganak_time"), ("td_width", "td_time"), ("td_width", "arjun_time"), ("td_width", "backbone_time"), ("compsK", "ganak_time")]:
    text= """
col1='"""+col1+"""'
col2='"""+col2+"""'
dirs = ["""+dirs+"""]

names = []
dfs = []
conn = sqlite3.connect('mydb.sql')

# Assign a color for each dirname
colors = plt.cm.tab10(np.linspace(0, 1, len(dirs)))

for d in dirs:
    # Run the SQL query and load the results into a DataFrame
    query = f"SELECT fname, {col1}, {col2}, ganak_call FROM data WHERE {col1} IS NOT NULL AND {col2} is not NULL and dirname='{d}' ORDER BY {colname}"
    df1 = pd.read_sql_query(query, conn)
    dfs.append(df1)
    names.append(d + " " + df1['ganak_call'][0])

conn.close()

# Plot the data
plt.figure(figsize=(10, 10))
for i, d in enumerate(dirs):
    plt.scatter(np.log10(dfs[i][col1]), np.log10(dfs[i][col2]), color=colors[i], label=names[i])

plt.title('Scatterplot')
plt.xlabel(f"{col1} (log10 scale)")
plt.ylabel(f"{col2} (log10 scale)")
plt.legend(loc='center left', bbox_to_anchor=(0, -0.2))
plt.show()
"""
    texts.append(text)

  # Create markdown cells
  markdown_cells = []
  for t in texts:
      markdown_cells.append(nbf.v4.new_code_cell(t))

  # Assign the cells to the notebook
  nb['cells'] = markdown_cells

  # Define the filename for the notebook
  filename = 'overview.ipynb'

  # Write the notebook to a file
  with open(filename, 'w') as f:
      nbf.write(nb, f)

  print(f"Notebook '{filename}' created successfully.")

create_notebook()
os.system("gnuplot "+gnuplotfn)
os.system("epstopdf run.eps run.pdf")
os.system("pdftoppm -png run.pdf run")
print("okular run.eps")
os.system("okular run.eps")


#### examples
# .mode table
# old = out-ganak-mccomp2324-14675861-0 -- mccomp2025
# new= out-ganak-mccomp2324-15010600-0 -- td start from 0

# old = "out-ganak-mccomp2324-21349-0", # TRILLIUM, ganak_7d97636055e_9104724fa_26d64aac (i.e. old run that was the fastest)
# new = "out-ganak-mccomp2324-35541-0", # TRILLIUM, fixing td starting from 0, now cutting disjoint components at toplevel for correct centroid

# select a1.dirname, a1.fname, a1.ganak_time as "old time", a2.dirname, a2.ganak_time as "new time" from data as a1, data as a2 where  a1.fname=a2.fname and a1.ganak_time is not null and a1.ganak_time is not null and a1.dirname like 'out-ganak-mccomp2324-21349-0' and a2.dirname like 'out-ganak-mccomp2324-35541-0' and a1.ganak_time < a2.ganak_time-100 and a1.ganak_time > 10 order by a1.ganak_time desc limit 50;

# cache hits
# select dirname, fname, round(ganak_time-arjun_time-td_time) as ganakT, cache_miss_rate from data where dirname like 'out-ganak-mccomp2324-14558670-1' and ganak_time is not null and cache_miss_rate is not null  and ganak_time-arjun_time > 100 order by cache_miss_rate asc limit 50;

# select a1.fname, a1.ganak_time as "gpmc time", a2.ganak_time from data as a1, data as a2 where  a1.fname=a2.fname and a1.ganak_time is not null and a1.ganak_time is not null and a1.fname not like 'mc%' and a1.solver like 'gpmc' and a2.solver like 'ganak' and a2.ganak_call like '%mode 6%' order by a2.ganak_time desc limit 50;

# select ganak_time, fname, dirname from data where ganak_time > 100 and dirname like '%out-ganak-gpmc-mei-14295899-0%' order by ganak_time desc limit 50;

# complex numbers
# select a1.dirname, a2.dirname, a1.fname, a1.ganak_time as "gpmc time", a2.ganak_time from data as a1, data as a2 where  a1.fname=a2.fname and a1.fname not like 'mc%' and a1.solver like 'gpmc' and a2.dirname like '%out-ganak-gpmc-mei-14295899-0%' and a2.ganak_call like '%mode 6%' order by a2.ganak_time desc limit 50;

# example fast problems
# select dirname, fname, ganak_time from data where dirname = 'out-ganak-mccomp2324-14309534-0' and ganak_time is not null and ganak_time < 200 and arjun_time < 20 and td_time < 5 order by ganak_time desc limit 50;

# small problems
# select dirname, fname, ganak_time,arjun_time,td_time from data where dirname = 'out-ganak-mccomp2324-14309534-0' and ganak_time is not null and ganak_time < 200 and arjun_time < 10 and td_time < 10 and new_nvars > 400 order by ganak_time desc limit 1000;

# new analyzer compare
# select a1.dirname, a1.fname, a1.ganak_time as "old time", a2.dirname, a2.ganak_time as "new time" from data as a1, data as a2 where  a1.fname=a2.fname and a1.ganak_time is not null and a1.ganak_time is not null and a1.dirname like 'out-ganak-mccomp2324-14362931-0' and a2.dirname like 'out-ganak-mccomp2324-14366334-0' and a1.ganak_time < a2.ganak_time-20 and a1.ganak_time < 500 order by a1.ganak_time desc limit 50;
