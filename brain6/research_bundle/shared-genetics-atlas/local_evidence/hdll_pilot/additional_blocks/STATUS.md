# Additional blocks: preparation stopped

No additional block was downloaded and no additional HDL-L result was produced. After the protocol was written, the preparation script could not create its working directory on `/Volumes/Brain6Work`: the operating system returned `No space left on device`. That working volume is backed by the full external SSD. The first `python3` attempt also failed because that interpreter lacks `requests`; `/opt/miniconda3/bin/python3` has the dependency and was used for the second attempt. No pre-existing data were deleted.

The existing chr3.44 reference remains available, but the saved R script requires `dplyr` and `data.table`, which were not available in the currently detected R library paths. No further trait-comparison model was run. Availability of input summary statistics or a reference is not itself a model result.

When space is available, run `prepare.py` with a Python interpreter containing `requests`; confirm storage before proceeding, then implement and validate the planned HDL-L runs. The script stops if the selected lead does not map uniquely to a BIM block and records the selected ZIP-member checksums. The preparation script has not yet been successfully executed end-to-end. The protocol remains prospective for these additional blocks; the three lead regions were previously selected and are not independent of earlier analyses.
