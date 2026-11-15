# bdnnf vs ddnnife: performance

From `scripts/bench_bdnnf.py` (ganak c1be5e32, ddnnife 7ae78cd). Most queries
run on five full MC2023 track-1 circuits (t015, t066, t076, t084, t102;
0.15–3.2 MB as BDNNF), the batch queries on 10000 samples/queries; 2-wise on
six shrunk instances. ddnnife reads `bdnnf to-d4`'s output. All answers agree.

| Query | bdnnf total | ddnnife total | bdnnf faster | Notes |
|---|---|---|---|---|
| count | 2.9 s | 16.3 s | 5/5 | |
| feature counts | 6.2 s | 28.9 s | 5/5 | |
| core | 6.1 s | 16.2 s | 5/5 | |
| atomic sets | 34.7 s | 57.5 s | 4/5 | |
| 2-wise | 34.7 s | 223 s | 6/6 | bdnnf's samples 1.2–4.6× larger |
| 10000 samples | 35.2 s | 35.2 s | 2/5 | |
| 10000 count queries | 177 s | 137 s | 0/4 | ddnnife timed out on t084 |
| 10000 SAT queries | 74.8 s | 30.9 s | 1/5 | |

Peak memory: bdnnf 20–250 MB, ddnnife 150 MB–3.2 GB.
