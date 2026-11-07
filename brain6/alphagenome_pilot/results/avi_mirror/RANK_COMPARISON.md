# AVI and GWAS ranking comparison

This descriptive comparison asks whether overall AVI and the existing trait-specific ABF PIPs prioritize the same selected variants. The comparison is restricted to each region’s union of reported 95% credible sets; it is not genome-wide and cannot validate shared causality. No machine-learning model was trained.

Within each union, variants are ranked separately by raw AVI and conditional ABF PIP. The comparison uses a rank boundary of five (or the full set if smaller); ties share the same minimum rank and can increase the number selected. Trait rankings are **not estimated** when the original per-trait association-posterior gate did not support a credible set. Empty PIPs remain missing.

At the chromosome 5 case-study locus, rs2431108 ranks 1 by AVI, 1 by insomnia ABF PIP, and 3 by ADHD ABF PIP among the 71 selected variants. A shared high ranking is a hypothesis-prioritization cue, not proof of a shared causal variant.

| Region, GRCh37 | Trait | Candidate variants | AVI at boundary | Trait at boundary | Overlap | Variants shared |
|---|---|---:|---:|---:|---:|---|
| chr5:103447968–104447968 | insomnia | 71 | 5 | 5 | 1 | rs2431108 |
| chr5:103447968–104447968 | adhd | 71 | 5 | 5 | 1 | rs2431108 |
| chr11:28076261–29076261 | insomnia | 88 | 5 | not estimated | not estimated | not estimated |
| chr11:28076261–29076261 | adhd | 88 | 5 | 6 | 1 | rs2585817 |
| chr14:98137479–99137479 | insomnia | 69 | 5 | not estimated | not estimated | not estimated |
| chr14:98137479–99137479 | adhd | 69 | 5 | 5 | 0 | none |
| chr3:117142005–118142005 | insomnia | 128 | 5 | 5 | 0 | none |
| chr3:117142005–118142005 | adhd | 128 | 5 | 5 | 1 | rs62264764 |
| chr5:87185500–88185500 | insomnia | 3 | 3 | not estimated | not estimated | not estimated |
| chr5:87185500–88185500 | adhd | 3 | 3 | 3 | 3 | rs324886;rs4916723;rs62369151 |
| chr7:113571035–114571035 | insomnia | 66 | 5 | 5 | 1 | rs1456031 |
| chr7:113571035–114571035 | adhd | 66 | 5 | 5 | 0 | none |

These selected sets vary in size and composition. PIP ranks are conditioned on a single-causal-variant ABF model, and the AVI mirror is independently spot-checked but not verified record by record. Apparent overlap can reflect the initial selection, linkage disequilibrium or chance. The numbers have no multiple-testing-adjusted significance, and the original Brain6 evidence classes remain unchanged.
