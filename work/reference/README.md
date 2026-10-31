# work/reference

Text substrate (LaTeX source and bibliographies) for papers cited by
`prior_insights` in `astra.yaml`, one directory per arXiv ID. Kept so that a
`scope:` field can point at a file and line number and a reader can check it
without re-fetching from arXiv.

Figures and compiled PDFs are deliberately not kept here. The PDF each quote
is verified against lives in the ASTRA paper cache
(`astra paper path <doi>`), populated by `astra paper add <doi>`; quotes are
checked with `astra validate astra.yaml --verify-evidence`.
