# bdnnf vs ddnnife: performance

From `scripts/bench_bdnnf.py` (ganak 23868bba, ddnnife 7ae78cd). Most queries
run on five full MC2023 track-1 circuits (t015, t066, t076, t084, t102;
0.15-3.2 MB as BDNNF); samples and the batch queries are 10000 each; 2-wise
runs on six shrunk instances. ddnnife reads `bdnnf to-d4`'s output. bdnnf
keeps the circuit in memory for multi-pass queries by default (`--cache`);
`--cache 0` streams the file on every pass. All answers agree, and every
t-wise sample passes `twise-check`.

| query | bdnnf total | bdnnf --cache 0 total | ddnnife total | bdnnf faster | timeouts (bdnnf / --cache 0 / ddnnife) |
|---|---|---|---|---|---|
| count | 2.41s | 2.41s | 15.9s | 5/5 | 0 / 0 / 0 |
| features | 5.37s | 5.36s | 28.6s | 5/5 | 0 / 0 / 0 |
| core | 5.41s | 5.34s | 15.7s | 5/5 | 0 / 0 / 0 |
| samples | 11.3s | 34.5s | 34.6s | 5/5 | 0 / 0 / 0 |
| count queries | 138.5s | 153.4s | 141.6s | 2/4 | 0 / 0 / 1 |
| SAT queries | 10.5s | 28.9s | 31.3s | 5/5 | 0 / 0 / 0 |
| atomic sets | 9.29s | 20.1s | 61.6s | 5/5 | 0 / 0 / 0 |
| 2-wise | 19.1s | 52.6s | 218.8s | 6/6 | 0 / 0 / 0 |

Peak memory: bdnnf 7-296 MB, with `--cache 0` 6-250 MB, ddnnife 38-3605 MB (one timed-out run not counted).

2-wise sample sizes (smaller is better):

| instance | bdnnf | ddnnife |
|---|---|---|
| auto1-s | 3109 | 809 |
| t015-s | 18 | 14 |
| t027-s | 451 | 94 |
| t028-s | 2672 | 2168 |
| t066-s | 381 | 184 |
| t003-s | 208 | 162 |
