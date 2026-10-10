# Pinned LDSC source attribution

The source files under `runtime/ldsc_munge_historical_text_v1/` derive from CBIIT LDSC commit `6c673952cee74bd5c57aef1555a03b1c015399a0`, retaining the original authorship and copyright notices. The exact upstream GNU GPL version3 license is supplied beside this notice. [The official repository](https://github.com/CBIIT/ldsc) identifies its GPL3 license; the local license bytes equal the LICENSE object at that pinned commit.

The separate historical text adapter changes the gzip input mode to `rt` for Python3 compatibility. Its unchanged numerical behavior and exact source differences are documented in the historical munge preparation and independent review receipts. No other scientific modification to that copied runtime is represented by this notice.

This notice and license were added outside the frozen runtime directory so that admitted input files and inventories remain unchanged. Original source, modifications and audit code remain available in the research package. This attribution does not grant rights to any GWAS source, genotype panel or participant data. Those assets remain excluded from Git under their separately assessed conditions.
