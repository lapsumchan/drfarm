# DrFARM: instructions for coding agents

This is the DrFARM package lane for Student W. Read existing ancestor/nested instructions and retain applicable project rules. Current user instructions and authorization take precedence over this draft.

Before implementation, read `docs/agent/PROJECT_INSTRUCTIONS.md`, `docs/agent/CHECKPOINT.md` and `docs/agent/PACKAGE_CHECKS.md`. For substantive mathematical/numerical changes, consult `docs/agent/PLUMBING_STANDARD.md` and select the applicable checks. Record target → rewrite → justification → finite computation → output meaning. Use `docs/agent/EVIDENCE_TEMPLATE.md` for a concise evidence record, with PASS / FAIL / NOT RUN / N/A and exact evidence. Do not equate written tests with executed results or tests with a statistical theorem.

Retrieve compatible prior work through `docs/agent/SOURCE_MAP.md` and `docs/agent/SOURCES.json` before recreating a derivation, algorithm or simulation. Read only relevant sections. Preserve source identities, mathematical assumptions, shape/order/units and failure history. If private references are unavailable in this checkout, continue unblocked work and request the exact missing reference only if needed; never fabricate access.

Preserve historical behavior as a named reference. Label corrections and experimental extensions. Profile measured bottlenecks on matched inputs and accuracy; retain convergence/failure diagnostics, RNG semantics and cache bindings. Avoid unrelated tests and repeated whole-archive audits for low-impact changes.

Keep one current package checkpoint at `docs/agent/CHECKPOINT.md`. Local work is authorized within the user task; external publication and native computation use the session's actual authorization. Do not change other project repositories or publish private notes by copying this kit into Git wholesale.

These instructions guide agents. They are not an executable mathematical checker; CI must separately run meaningful package tests when implemented.

The `docs/agent/` support files are private project context supplied with the maintenance review bundle and are intentionally excluded from the public package/source patch. Restore that context from the bundle before later agent work; missing private references do not block independent maintenance.
