# Protected Track B: technical legacy-to-v3 schema crosswalk

The archived B ledger is intact at SHA-256 `f0eb8b9e8a5b3b4a0ea6151ed75473b0be3d311e463824889a050d83beaaac2e`. A read-only scan checked all **5,514,399** rows against its 17-field historical schema and the published v3 10-field template. Every row was `TESTED`, had the declared within-pair denominator, finite statistics and valid P/q values; the ledger was sorted with no adjacent duplicate key. This diagnostic did not create a candidate v3 output or import B.

| v3 field | Proposed historical source or deterministic expression |
|---|---|
| `SNP` | `SNP` |
| `CHR` | `CHR` |
| `BP` | `BP` |
| `Z1` | `Z1` |
| `Z2` | `Z2` |
| `P_PLACO` | `P_PLACO_PLUS` |
| `Q_WITHIN_PAIR` | `PLACO_BH_Q` |
| `Q_BONFERRONI_ACROSS_PAIRS` | `min(1, 5 * PLACO_BH_Q)` |
| `status` | `analysis_status (TESTED only)` |
| `headline` | `P_PLACO_PLUS < 1e-8` |

The historical B headline threshold was P≤2.5e-8 across the original two primary Track B pairs (441 rows). The frozen Brain6 v3 five-track headline is P<1e-8, which would select **389** B rows under the proposed adapter. The historical within-pair BH q<0.05 count is 3,755; multiplying that q by five and capping at one yields 1,472 rows with proposed five-track q<0.05. These are **diagnostic counts, not admitted v3 family results**.

The conversion would copy SNP, chromosome, position and Z scores unchanged; rename `P_PLACO_PLUS` and `PLACO_BH_Q`; derive the five-track q and headline using the frozen v3 rules; and retain `TESTED`. The original A1/A2, P1/P2, statistic, error and family-size fields must remain preserved in the source ledger/provenance, even though the v3 display schema omits them. A read-only admission validator must additionally verify the exact historical gate/cleanup lineage, input/source identity, result hashes, and scientific compatibility before any import. No current frozen rule grants that admission; the protected slot remains blocked.
