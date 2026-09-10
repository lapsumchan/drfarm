# Cite DrFARM and preserve a reproducible release

The GitHub repository is the project's home. Its README is the main introduction;
the optional pkgdown site renders that README, R help and articles from the same
repository. “Function reference” means API help. Scholarly and software citations
belong here and in `citation("drfarm")`.

| Resource | Purpose | What to identify |
|---|---|---|
| GitHub repository | Evolving source, README, issues and review | Full commit for a development checkout |
| Package documentation site | Readable view of the repository's documentation | Package version and matching source download |
| Zenodo | Citable snapshots of selected software releases | Version DOI for the code actually used |

## Paper and historical archive

Chan, L. S., Li, G., Fauman, E. B., Yin, X., Laakso, M., Boehnke, M. &
Song, P. X. K. (2025). **DrFARM: identification of pleiotropic genetic variants
in genome-wide association studies.** *Nature Communications*, **16**, 5789.
[doi:10.1038/s41467-025-60439-4](https://doi.org/10.1038/s41467-025-60439-4).

Reference 66 in the paper identifies **lapsumchan/drfarm: 0.1.0**, published
20 April 2025 on Zenodo: [doi:10.5281/zenodo.15252156](https://doi.org/10.5281/zenodo.15252156).
That record links to the repository's [0.1.0 source tag](https://github.com/lapsumchan/drfarm/tree/0.1.0).
Its recorded creator is `lapsumchan`; the package's authors are Lap Sum Chan,
Gen Li and Peter X.K. Song. The seven-author article and the software have
distinct authorship records.

The verified [concept DOI, 10.5281/zenodo.15252155](https://doi.org/10.5281/zenodo.15252155),
identifies the software's version series. Use it for discovery; use the
**version DOI** to identify an analysis's exact release.
[Zenodo explains the distinction](https://support.zenodo.org/help/en-gb/1-upload-deposit/97-what-is-doi-versioning).

## Cite the code actually used

Cite the method paper and the installed software version. For historical 0.1.0,
use its version DOI above. For this unpublished **0.1.0.9002** development
candidate, record the repository URL and full source commit alongside the package
version; it has no newly assigned release DOI. The 0.1.0 archive does **not**
contain the weighted coefficient option or Gaussian ECM reference.

```sh
git rev-parse HEAD
```

```r
packageVersion("drfarm")
citation("drfarm")
sessionInfo()
```

For a source archive without Git history, record its SHA256 and the source commit
from its release or review manifest. Package version alone can be ambiguous
between development commits. A paper citation or an archive DOI does not validate
new algorithms or confer the paper's inferential claims on the ECM reference.

## Selected releases, not every commit

Git manages day-to-day history. Our proposed archive policy is to preserve
reviewed, tested releases in the existing Zenodo version series, keeping the
paper's 0.1.0 snapshot identifiable. Code changes after a released analysis
belong in a new release, with a new version DOI. Routine README edits do not
need a DOI each.

For each scientific release:

1. Bind the package version, Git tag, full commit and source archive SHA256.
2. Run the original minimal data example and applicable checks; retain output,
   warnings and actual stopping reasons, dependency versions and seeds.
3. Include runnable reproduction commands and identify the data and result
   artifacts. Keep private or restricted research data outside the public archive.
4. Review citation metadata, licensing and release notes; publish only after the
   maintainer approves that concrete release.
5. Record the new version DOI in the release and analysis manifests, verify its
   connection to the existing series, and keep the historical DOI unchanged.

The [GitHub–Zenodo integration](https://docs.github.com/en/repositories/archiving-a-github-repository/referencing-and-citing-content)
can archive GitHub releases; it is not continuous mirroring of every commit.
Before enabling or changing it, verify the repository's existing connection and
the maintainer's control of this Zenodo series. No integration or publication
has been changed in this documentation sprint.

`CITATION.cff` supplies machine-readable GitHub citation metadata; `inst/CITATION`
supplies R's citations. Keep their authors, version and article DOI consistent.
We use CFF without a second `.zenodo.json`: when both exist, Zenodo
[uses `.zenodo.json` and ignores CFF metadata](https://help.zenodo.org/docs/github/describe-software/citation-file/).
Do not assign the historical DOI to newer software metadata.

Before the next archive, reconcile one existing metadata difference: the 0.1.0
Zenodo record says GPL-3.0-only, while package `DESCRIPTION` says GPL (>= 3).
This documentation preserves the package license and the historical record;
any archive metadata correction requires checking the original license notices.
