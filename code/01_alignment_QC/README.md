# Genome_QC

This directory contains the main alignment quality-control workflow for comparative coding-sequence alignments.

## Main tasks

- remove low-coverage alignment columns
- detect and remove short suspicious fragments
- summarize fragment structure
- remove highly incomplete sequences
- translate, codon-align, concatenate, and sample alignments

## Main workflow

The clearest bundled workflow is recorded in `option.sh`:

1. `RemoveGap_cds.pl`
2. `RemoveFragmentSeq_cds.pl`
3. `RemoveFragmentSeq_cds.pl` again with stricter thresholds
4. `RemoveGap_cds.pl` again after fragment cleanup
5. `FragmentCount.pl` and optionally `FragmentScreen.pl`
6. `RemoveSmallFrag.pl`
7. `RemoveEmptySequences.pl`

See the repository-level workflow notes in `../../docs/workflow_overview.md` for a condensed command example.

## Key scripts

- `RemoveGap_cds.pl`: remove codon columns below a species coverage threshold
- `RemoveFragmentSeq_cds.pl`: mask short fragments flanked by large gaps
- `FragmentCount.pl`: summarize fragment structure for each sequence
- `FragmentScreen.pl`: identify specific small fragments for downstream removal
- `RemoveSmallFrag.pl`: remove fragments listed by the screening step
- `RemoveEmptySequences.pl`: drop sequences with too much missing data
- `Translate.pl` and `CodonAlign.pl`: translation and codon-alignment helpers
- `Concatenate.pl`: concatenate alignments

## Local files

- `Standard.txt`: codon translation table used by several scripts
- `option.sh`, `hg38.option.sh`: command examples / workflow notes
- `QC_compare/`: example or reference comparison outputs

## Cautions

- Many scripts assume they are run from inside this directory.
- Several scripts open `Standard.txt` with a relative path.
- Output files such as `Delete_info.txt`, `gap_info.txt`, and `empty_sequence_info.txt` may be written to the current working directory.
