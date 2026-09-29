# Historical-to-current Track B gate delta

The two input_gate.lock.json objects have identical fields and values except readiness_gate_sha256. Both readiness hashes verify against their respective TSV bytes. The TSVs have the same component rows and columns; exactly two blockers cells differ. In both, current state appends ACTIVE_LAVA_DOWNLOAD_BLOCKS_LD_MATERIALIZATION to the historical string. All other fields, including the five dense input identities and contract lock, are identical.

This narrows the v1 gate mismatch to readiness blocker text, but it does not admit historical B to frozen Brain6 v3. Any exception for this historical gate lineage still requires the proposed, explicitly authorized admission amendment and independent validator. See gate_delta_audit.json for exact before/after values and hashes.
