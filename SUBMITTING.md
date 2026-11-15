# Submitting (team notes; not in the zip)

This file and `.gitattributes` are marked `export-ignore`, so `git archive` leaves them out.

## Building the zip

The zip is this branch (`submission`), archived:

```bash
git archive --format=zip --prefix=CS690_Project6/ -o CS690_Project6.zip submission
```

Only committed files go in, so commit everything first.

## Adding your files

Commit into `slides/`, `report/` or `code/` on the `submission` branch, then build the zip again. Keep one final
version of each file (no `_v2` or `final_final` copies); put slide sources next to the PDF in `slides/`, and the
report as `report/report.pdf`.
