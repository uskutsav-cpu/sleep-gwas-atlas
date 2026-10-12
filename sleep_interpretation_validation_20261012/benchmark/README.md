# Provisional AI interpretation benchmark

The authored inventory contains 60 new final tasks, six separate development tasks, ten balanced strata and a condition-specific private key. B/C carry identical canonical information. A uses the authentic frozen Human GWAS ATLAS resource and receives correct abstention gold when exact atlas source identity is not authenticated. The protected numeric evidence and gold are on the external SSD. Human scientific certification remains pending.

`PROTOCOL_AND_SAP.md` and `execution_registration.json` define an initial 90-cell, ten-final-task feasibility experiment using three OpenAI Codex aliases. It is not completion of the 60-task benchmark, a comparison of independent providers, or evidence of general interface superiority. No responses are produced by the builder, preparation tool, validator or scorer.

Run order is construction, per-cell preparation, mechanical validation, independent review, freeze, actual model execution, then provisional scoring. Do not supply gold or other-condition responses to any solver. Every fresh agent reads exactly its assigned complete JSON packet. Never substitute automatic parsing, scripted extraction, gold copying, or a code-generated answer for an actual model response.

For each registered cell, give a fresh `fork_turns=none` agent this instruction, filling only the assigned paths and registered values:

> Read only the complete assigned packet at PACKET_PATH. Do not read any gold/key, other condition packet, other response, or external source. Use only the packet's system_prompt, user_prompt, evidence and response_schema to produce your own answer. Tool access is restricted by instruction to reading this packet and writing RESPONSE_PATH; do not use code to derive or copy an answer. Return one response JSON envelope and save it at RESPONSE_PATH. Do not browse, delegate or message another agent. Record the actual alias requested and any observed backend/version; if unavailable use UNKNOWN. Do not invent temperature, sampling seed, timing, tokens or proof of isolation. Report incomplete packet reads, truncation, compaction, or unexpected access honestly.

The envelope shape is:

```json
{
  "task_id": "REGISTERED_TASK_ID",
  "condition": "REGISTERED_A_B_OR_C",
  "model_alias": "ACTUAL_REGISTERED_ALIAS",
  "provider": "OpenAI",
  "requested_reasoning_effort": "medium",
  "observed_backend_version": "UNKNOWN",
  "execution_status": "COMPLETED_RESPONSE",
  "response": {"status": "ANSWER_STATUS", "other_required_fields": "MODEL_ANSWER"},
  "attempt_index": 1,
  "temperature": "UNKNOWN_UNCONTROLLABLE",
  "sampling_seed": "UNKNOWN_UNCONTROLLABLE",
  "timing_seconds": null,
  "tokens": null,
  "packet_sha256": "REGISTERED_PACKET_SHA256",
  "access_declaration": {"read_paths": ["PACKET_PATH"], "write_paths": ["RESPONSE_PATH"], "other_access": [], "full_packet_read": true, "truncated": false, "compactions": 0},
  "observable_access_trace_status": "UNKNOWN_UNLESS_INDEPENDENTLY_CHECKED"
}
```

`response` must contain the actual closed structured model answer, without the illustrative `other_required_fields` placeholder. Use `TRANSPORT_FAILED`, `UNEVALUABLE` or `INTERRUPTED` execution status for an attempt that never yields a complete answer; preserve the attempt rather than creating an answer. A complete invalid JSON answer remains a completed response and is scored invalid. The orchestrator records timestamps and actual observable execution/tool evidence separately. Self-report does not certify isolation.

Aggregate actual first definitive response envelopes into a JSON list or JSONL file. Preserve all original attempts and any explicitly allowed transport retry. Pass the frozen registration, public final task inventory and private key to `score_benchmark.py`, writing results only into the private benchmark directory. Its strict correctness rate is mechanical frozen-key agreement, provisional pending human review. It reports A answerability/abstention separately and omits inferential intervals for the ten-task feasibility.

After `FREEZE_RECEIPT.json` exists, do not regenerate or edit frozen files. Corrections require an explicit versioned amendment preserving original evidence, gold and responses. The immutable receipt records full inventory, code, protocol, private key and all packet hashes before the first final response.
