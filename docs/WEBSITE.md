# Build and review the DrFARM documentation

The existing GitHub repository and its README are the primary user-facing home.
This optional site is a generated view of the same source, maintained in the same
repository. Editing README changes the GitHub landing page after publication;
rebuilding pkgdown updates the website. No second implementation or separate
website repository is required. The “Function reference” menu is API help;
“Cite DrFARM” contains scholarly/software citations.

The site uses **pkgdown 2.0.7**, the package's README, Rd help, and two articles.
Numerical R/Rcpp source is unchanged by this documentation sprint. The site is a
review candidate; the configured destination is **https://lapsumchan.github.io/drfarm/**.
It has not been published by this local build.

## Local build

From the candidate's complete source checkout, with R 4.3.3, Pandoc 3.1.3 and a
C++ toolchain available:

```sh
Rscript --vanilla tools/install-site-dependencies.R
Rscript --vanilla tools/build-site.R
```

The second line is the site-build command. The first reuses the recorded
60-package dependency manifest and adds 17 pinned documentation packages. These
pins are a reproducibility target, not a claim that newer versions are unsupported.
The installer accepts `--library=/absolute/path` and `--verify-only`; activate a
custom library through `R_LIBS_USER` before building. Ubuntu build dependencies
include libcurl4-openssl-dev, libssl-dev, libxml2-dev, libfontconfig1-dev,
libfreetype6-dev, libharfbuzz-dev, libfribidi-dev, libpng-dev, libtiff-dev,
libjpeg-dev, and gfortran. The workflow installs these explicitly.

Rendering disables optional package-timeline and cross-package metadata downloads;
installing dependencies is a separate, explicit network step.
The builder stages a public-file allowlist, builds and installs the source
package in `site-build/library`, and renders `_site/`. It never reads private
project handoffs. Build/install logs, installed exports and session information
are written to `site-build/`. The site includes its matching source archive and
standalone example scripts in `downloads/`; this avoids linking new APIs to the
older public-main installation. Generated output is ignored by Git.

Open `_site/index.html` to browse the static preview, or serve `_site/` with a
local static-file server. An optional development preview is provided by
`npm ci` followed by `npm run dev` (Node 20.19+ or 22.12+; Vite is pinned in
package-lock.json). This serves pkgdown output; Node is not needed to install
DrFARM or build its documentation. MathJax is loaded from a CDN for typesetting; an internet
connection is needed for equations in this first preview. GitHub renders the
README's native dollar-delimited math independently. Website source links are
rewritten only in the staging copy of README.

## Check the user journey

Run the full original example against the freshly installed candidate:

```sh
Rscript --vanilla tools/reproduce.R --mode full \
  --output site-build/original-example --label documentation-candidate \
  --library site-build/library
R_LIBS=site-build/library Rscript --vanilla inst/examples/quickstart.R
R CMD check --no-manual site-build/drfarm_0.1.0.9002.tar.gz
```

The original example uses 100 remMap and 25 DrFARM tuning cells and both
historical inference functions. Completion and numerical compatibility are
separate from convergence: preserve `loss_increase` warnings and returned
status. The article collects repeated warnings into a count, without altering
the computation. The four expensive 500-cycle ECM comparisons are retained,
source-bound evidence and are not a site build step.

Run `python3 tools/check-site.py _site` to check local links and anchors.
Check desktop and narrow layouts, menu navigation, function reference links,
citation, changelog, source archive and example downloads. Read the methods
article before interpreting any fitted or inferential output. Windows and macOS
have no local validation receipt in this sprint.

## Publication path

`.github/workflows/pkgdown.yaml` builds on pull requests, pushes to main and
manual dispatch, and uploads a review artifact. **Publishing is manual only**:
a dispatch with `publish=true` on main runs the separate Pages deployment job.
It requires a prior decision to publish this candidate, an approved merge, and
repository Pages configured to use GitHub Actions. The default input is false;
a local build or an ordinary merge does not enable public Pages.

The proposed site is the package documentation at
`https://lapsumchan.github.io/drfarm/`, not Student W's personal website.
Before publication, review the local implementation commits as well as this
sprint's documentation diff: public main lacks the weighted coefficient path,
Gaussian ECM reference, finite-budget diagnostics, tests and maintenance docs.
Update the candidate/public-main status text when the publication state changes.
No passing Actions run or deployment is claimed by the local receipt.

References: [pkgdown introduction](https://pkgdown.r-lib.org/articles/pkgdown.html),
[pkgdown build_site](https://pkgdown.r-lib.org/reference/build_site.html), and
[GitHub's custom Pages workflow](https://docs.github.com/en/pages/getting-started-with-github-pages/using-custom-workflows-with-github-pages).
