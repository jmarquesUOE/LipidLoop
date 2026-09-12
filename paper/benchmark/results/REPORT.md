# NIST SRM 1950 benchmark

Consensus = Bowden 2017 (NIST interlaboratory exercise); Ref = >= 5 labs and COD <= 40 % (the LipidQC benchmark set), COD>40 = >= 5 labs with COD > 40 %, Ref + COD>40 = the paper's 339; Inf = 3-4 labs. Out-of-scope entries (eicosanoids, bile acids, cholesterol, S1P, free fatty acids) are excluded from the denominators.

| run | species (MS2 / all) | Ref recall MS2 | Ref recall MS2+RT | 339 recall MS2 | 339 recall MS2+RT | Inf recall MS2+RT | listed (MS2) | listed (RT model) | pooled rho (n) | median rep CV raw / column-normalised / interlab COD |
|---|---|---|---|---|---|---|---|---|---|---|
| v1 | 326 / 362 | 187/231 (81%) | 194/231 (84%) | 223/310 (72%) | 231/310 (75%) | 7/36 (19%) | 269/295 (91%) | 14/34 (41%) | 0.63 (178) | 6.6 % / 7.9 % / 24.0 % (n=178) |
| v2 | 350 / 391 | 191/231 (83%) | 196/231 (85%) | 225/310 (73%) | 233/310 (75%) | 6/36 (17%) | 271/310 (87%) | 17/37 (46%) | 0.59 (185) | 10.4 % / 7.1 % / 24.0 % (n=185) |
| v3 | 439 / 473 | 206/231 (89%) | 208/231 (90%) | 253/310 (82%) | 255/310 (82%) | 8/36 (22%) | 315/375 (84%) | 11/26 (42%) | 0.57 (200) | 6.5 % / 9.5 % / 24.0 % (n=200) |
| v4 | 326 / 411 | 188/231 (81%) | 199/231 (86%) | 222/310 (72%) | 235/310 (76%) | 8/36 (22%) | 271/291 (93%) | 27/75 (36%) | 0.54 (179) | 8.7 % / 6.3 % / 24.0 % (n=179) |
| v5 | 361 / 427 | 191/231 (83%) | 198/231 (86%) | 227/310 (73%) | 237/310 (76%) | 7/36 (19%) | 278/314 (89%) | 19/62 (31%) | 0.6 (183) | 18.7 % / 13.1 % / 24.0 % (n=183) |

## Consensus recall by family (Ref + COD>40, MS2 / MS2+RT of n)

| family | n | v1 | v2 | v3 | v4 | v5 |
|---|---|---|---|---|---|---|
| PC | 63 | 53 / 56 | 53 / 55 | 59 / 60 | 51 / 54 | 52 / 54 |
| TG | 59 | 46 / 48 | 47 / 49 | 55 / 55 | 43 / 47 | 44 / 46 |
| SM | 38 | 35 / 35 | 35 / 36 | 37 / 37 | 35 / 36 | 36 / 36 |
| PE | 35 | 22 / 24 | 24 / 24 | 25 / 25 | 22 / 22 | 24 / 24 |
| LPC | 25 | 20 / 20 | 20 / 21 | 20 / 20 | 20 / 22 | 19 / 19 |
| DG | 24 | 7 / 7 | 5 / 5 | 10 / 10 | 7 / 7 | 8 / 10 |
| CE | 19 | 15 / 16 | 13 / 14 | 15 / 16 | 15 / 17 | 16 / 17 |
| PI | 15 | 8 / 8 | 10 / 10 | 12 / 12 | 11 / 11 | 9 / 11 |
| Cer | 15 | 7 / 7 | 8 / 9 | 10 / 10 | 8 / 9 | 9 / 10 |
| LPE | 8 | 7 / 7 | 6 / 6 | 6 / 6 | 6 / 6 | 6 / 6 |
| HexCer | 5 | 3 / 3 | 4 / 4 | 4 / 4 | 4 / 4 | 4 / 4 |
| PG | 3 | 0 / 0 | 0 / 0 | 0 / 0 | 0 / 0 | 0 / 0 |
| PS | 1 | 0 / 0 | 0 / 0 | 0 / 0 | 0 / 0 | 0 / 0 |

## Precision proxy by family (our species found in any published list / our species), MS2 only

| family | v1 | v2 | v3 | v4 | v5 |
|---|---|---|---|---|---|
| AC | 5/7 | 4/8 | 6/11 | 6/7 | 5/7 |
| CE | 16/16 | 13/15 | 16/17 | 16/16 | 18/18 |
| Cer | 10/10 | 11/11 | 14/14 | 12/12 | 15/15 |
| DG | 7/7 | 5/5 | 10/10 | 7/7 | 8/8 |
| FA | 2/2 | 2/2 | 4/4 | 2/2 | 2/2 |
| HexCer | 3/3 | 5/5 | 6/7 | 5/6 | 5/5 |
| LPC | 21/22 | 20/22 | 21/23 | 21/22 | 19/20 |
| LPE | 7/7 | 6/8 | 8/12 | 6/8 | 7/9 |
| LPI | 1/1 | 4/4 | 5/5 | 4/4 | 5/5 |
| PA | - | 0/1 | 0/1 | 0/1 | 0/1 |
| PC | 72/80 | 72/81 | 81/94 | 69/74 | 70/82 |
| PE | 22/24 | 24/24 | 26/30 | 23/23 | 24/27 |
| PI | 8/8 | 10/10 | 12/13 | 11/11 | 9/9 |
| SM | 39/41 | 38/44 | 41/50 | 38/40 | 39/45 |
| TG | 56/67 | 57/70 | 65/84 | 51/58 | 52/61 |

## Within-family Spearman rho, log area vs log consensus (Ref tier, MS2, one polarity per family, families with >= 5 matches)

| family | v1 | v2 | v3 | v4 | v5 |
|---|---|---|---|---|---|
| CE | 0.75 (n=11, Pos) | 0.74 (n=11, Pos) | 0.73 (n=11, Pos) | 0.77 (n=11, Pos) | 0.72 (n=12, Pos) |
| Cer | 0.9 (n=5, Pos) | 0.94 (n=6, Neg) | 0.93 (n=8, Pos) | 0.94 (n=6, Pos) | 0.94 (n=6, Neg) |
| LPC | 0.82 (n=20, Pos) | 0.79 (n=20, Pos) | 0.8 (n=20, Pos) | 0.81 (n=20, Pos) | 0.79 (n=19, Pos) |
| LPE | 0.43 (n=6, Pos) | 0.77 (n=6, Pos) | 0.43 (n=6, Pos) | 0.26 (n=6, Pos) | 0.54 (n=6, Pos) |
| PC | 0.8 (n=43, Pos) | 0.87 (n=44, Pos) | 0.82 (n=49, Pos) | 0.87 (n=42, Pos) | 0.88 (n=46, Pos) |
| PE | 0.14 (n=19, Pos) | 0.71 (n=20, Neg) | 0.22 (n=21, Pos) | 0.05 (n=18, Pos) | 0.72 (n=20, Neg) |
| PI | 0.71 (n=6, Pos) | 0.73 (n=10, Neg) | 0.77 (n=10, Neg) | 0.77 (n=10, Neg) | 0.76 (n=8, Neg) |
| SM | 0.94 (n=27, Pos) | 0.95 (n=28, Neg) | 0.95 (n=30, Neg) | 0.95 (n=28, Neg) | 0.96 (n=28, Neg) |
| TG | 0.7 (n=35, Pos) | 0.76 (n=34, Pos) | 0.7 (n=38, Pos) | 0.53 (n=32, Pos) | 0.75 (n=32, Pos) |

## Out of scope (Bowden entries no spectral library can name)

11,12-DiHETrE, 11-HDoHE, 12,13-DiHOME, 12,13-EpOME, 12-HETE, 12-HHTrE, 13-HOTrE, 14-HDoHE, 15-HETE, 17-HDoHE, 18-HEPE, 20-HETE, 5,6-EET, 5-HEPE, 5-HETE, 8-HETE, 8-HETrE, 9,10-DiHOME, 9-HEPE, 9-HETE, 9-HODE, 9-OxoODE, CA, CDCA, Cholesterol, DCA, FFA 16:0, FFA 16:1, FFA 17:0, FFA 17:1, FFA 18:0, FFA 18:1, FFA 18:2, FFA 18:3, FFA 20:1, FFA 20:3, FFA 20:4, FFA 20:5, FFA 22:5, FFA 22:6, GCA, GCDCA, GDCA, GLCA, GUDCA, LCA, MCA, PGE2, S1P, TCA, TCDCA, TDCA, TLCA, TLCA-S, UDCA, dhS1P

## v1: MS2 species in covered families that no list carries (46)

AC 18:0, AC 18:2, Alkanyl-TG O-50:1, Alkenyl-TG P-52:1, Alkenyl-TG P-52:2, DG 46:8, LysoPC 19:0, PC 37:1, PC 39:3, PC 40:9, PC 42:2, PC 42:7, PC 42:9, PC O-40:3, PC O-42:2, PE 39:4, PE O-38:3, PE O-40:4, Plasmanyl-PC O-33:0, Plasmanyl-PC O-42:4, SM d36:4, SM d38:4, SM d40:4, SM d42:7, SM d42:8, SM d44:4, TG 18:1_18:0_19:0, TG 37:6, TG 45:0, TG 55:2, TG 55:4, TG 56:1, TG 57:1, TG 57:2, TG 57:3, TG 57:4, TG 57:5, TG 58:1, TG 59:3, TG 59:5, TG 59:6, TG 60:2, TG 60:3, TG 60:4, TG 68:9, TG 70:12

## v2: MS2 species in covered families that no list carries (59)

AC 18:2, AC 18:3, AC 24:0, AC 26:0, Alkanyl-TG O-52:1, Alkenyl-TG P-50:0, Alkenyl-TG P-50:1, Alkenyl-TG P-52:1, Alkenyl-TG P-52:2, Alkenyl-TG P-54:2, CE 17:2, CE 21:0, DG 35:7, DG 42:5, LysoPC 19:0, LysoPC 22:3, LysoPE 18:3, LysoPE 20:5, PA 20:6_26:1, PC 28:2, PC 37:1, PC 40:1, PC 42:2, PC 42:3, PC 42:7, PC O-39:0, PC O-40:3, PC O-42:2, PE 39:3, Plasmanyl-PC O-37:1, Plasmanyl-PC O-42:4, SM d36:4, SM d38:4, SM d38:5, SM d39:0, SM d40:4, SM d41:0, SM d42:0, SM d42:7, SM d42:8, SM d44:4, SM d44:5, TG 18:1_18:1_24:0, TG 45:0, TG 45:8, TG 47:10, TG 53:0, TG 55:2, TG 55:4, TG 56:0, TG 56:1, TG 57:3, TG 57:4, TG 57:5, TG 57:7, TG 58:1, TG 59:3, TG 59:5, TG 60:3

## v3: MS2 species in covered families that no list carries (75)

AC 18:0, AC 18:2, AC 24:0, AC 26:0, AC 26:1, Alkanyl-TG O-52:1, Alkenyl-TG P-50:0, Alkenyl-TG P-50:1, Alkenyl-TG P-52:1, Alkenyl-TG P-52:2, Alkenyl-TG P-54:1, Alkenyl-TG P-54:2, CE 17:2, CE 24:4, DG 44:6, GlcCer[NDS] d42:3, LysoPC 19:0, LysoPC 22:3, LysoPE 17:0, LysoPE 18:3, LysoPE 20:5, LysoPE 22:4, PA 44:7, PC 36:7, PC 37:1, PC 37:7, PC 39:3, PC 40:1, PC 41:3, PC 42:2, PC 42:3, PC 42:7, PC O-40:3, PE 38:7, PE O-36:0, PE O-36:1, PI 35:2, Plasmanyl-PC O-33:0, Plasmanyl-PC O-42:2, Plasmanyl-PC O-42:4, Plasmenyl-PE P-20:0_18:2, SM d36:4, SM d38:4, SM d39:0, SM d40:4, SM d40:7, SM d41:0, SM d42:0, SM d42:7, SM d42:8, SM d44:4, SM d44:5, TG 18:1_18:1_23:0, TG 45:0, TG 47:10, TG 47:5, TG 53:0, TG 55:1, TG 55:2, TG 55:4, TG 55:5, TG 56:0, TG 56:1, TG 57:2, TG 57:3, TG 57:4, TG 57:5, TG 57:6, TG 58:1, TG 59:3, TG 59:5, TG 59:6, TG 60:2, TG 60:3, TG 60:6

## v4: MS2 species in covered families that no list carries (68)

AC 17:4, AC 26:0, Alkanyl-TG O-52:1, Alkenyl-TG P-50:0, CE 19:0, CE 21:0, CE 24:4, DG 44:7, GlcCer[NDS] d42:3, LysoPC 19:0, LysoPE 18:3, LysoPE 20:5, PA 20:6_26:1, PC 33:5, PC 35:7, PC 37:1, PC 42:7, PC O-40:3, PC O-42:2, PC O-42:4, Plasmanyl-PC O-33:1, Plasmanyl-PC O-38:0, Plasmanyl-PC O-46:1, SM d31:2, SM d34:3, SM d36:4, SM d39:0, SM d40:4, SM d41:0, SM d42:0, SM d44:4, TG 18:1_18:0_20:0, TG 42:8, TG 43:5, TG 44:9, TG 45:0, TG 45:8, TG 47:10, TG 53:0, TG 55:0, TG 55:1, TG 55:2, TG 55:4, TG 55:5, TG 56:0, TG 57:1, TG 57:2, TG 57:3, TG 57:4, TG 57:5, TG 57:6, TG 58:0, TG 58:1, TG 59:4, TG 59:5, TG 59:7, TG 60:1, TG 60:2, TG 60:3, TG 60:4, TG 60:6, TG 60:7, TG 61:3, TG 62:2, TG 62:3, TG 62:4, TG 62:8, TG 66:6

## v5: MS2 species in covered families that no list carries (79)

AC 17:3, AC 18:0, AC 18:2, Alkanyl-TG O-50:0, Alkanyl-TG O-50:2, Alkenyl-TG P-50:0, Alkenyl-TG P-52:1, CE 19:0, CE 21:0, CE 24:4, LysoPC 19:0, LysoPC 22:3, LysoPE 18:3, LysoPE 20:5, PA 20:6_26:1, PC 37:1, PC 37:7, PC 40:1, PC 40:9, PC 42:2, PC 42:7, PC 42:9, PC 44:11, PC 44:5, PC 46:0, PC O-42:2, PC O-44:3, PE 38:7, PE 50:6, PE O-36:1, Plasmanyl-PC O-33:0, Plasmanyl-PC O-38:0, Plasmanyl-PC O-40:1, Plasmanyl-PC O-40:3, Plasmanyl-PC O-42:4, Plasmanyl-PC O-46:0, Plasmenyl-PE P-20:0_18:2, SM d34:3, SM d36:4, SM d38:4, SM d39:0, SM d40:4, SM d41:0, SM d42:0, SM d42:8, SM d44:4, SM d44:5, TG 18:1_18:1_23:0, TG 42:7, TG 42:8, TG 43:5, TG 44:9, TG 45:0, TG 47:10, TG 53:0, TG 55:0, TG 55:1, TG 55:2, TG 55:4, TG 56:0, TG 56:1, TG 57:1, TG 57:2, TG 57:3, TG 57:4, TG 58:1, TG 58:16, TG 59:3, TG 59:5, TG 59:6, TG 60:2, TG 60:3, TG 60:4, TG 60:5, TG 60:6, TG 60:7, TG 61:4, TG 62:4, TG 70:12
