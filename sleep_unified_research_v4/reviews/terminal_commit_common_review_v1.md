The frozen v1 helper is rejected for new controller use. Two fresh metadata-only witnesses reproduced both false admissions, and the two original witnesses remain preserved and accepted by its read-only consumer. The inspected helper SHA256 is `4ed063eb62f5768455cf38957f38ff6b1f32f2d9649757d3c8544f63d67e7e90`.

A dangling `receipt.json.failure.json` symlink leaves `.exists()` false. V1 checks only `.exists()` at producer line 74 and consumer lines 101 and 103. The `broken_addendum` witness returns `commit=True` and consumer acceptance despite the failure-marker directory entry. This violates the intended any-failure-addendum veto.

V1 copies an empty stage binding at line 40 and accepts it at consumer line 98. The `empty_stage_binding` witness returns `commit=True` and consumer acceptance with `{}` as the stage identity. Exact nonempty receipt hashes alone do not bind a particular reviewed stage.

`terminal_commit_review_controls_v1.py` creates fresh fixtures exclusively in `reviews/terminal_commit_controls_v1_recheck`; it imports the frozen helper, never patches it, and separately rechecks the two untouched fixtures in `reviews/terminal_commit_controls_v1`. The receipt records all eight regular-file hashes and two literal symlink targets across both fixture sets. Initial and final helper hashes match. These controls establish the two rejection findings; they are not a complete audit of v1.

The additive v2 review is separate. No original evidence, helper, controller, research source or outcome was altered. GWAS/reference body reads, decompressions, workers and fits were all zero. This report grants no operational execution admission.
