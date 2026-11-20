# AstroAmase Lite

AstroAmase Lite fits astronomical spectra using a selected subset of the Astro AMASE molecular catalogs. This branch contains a staged research draft. Open the [guide](docs/index.html) or download the [Colab notebook](docs/downloads/AstroAmase_Lite_Integrated.ipynb).

The standard CPU workflow runs in the user's own Google Colab session. It gets pinned source files from the public data deposit, builds the active database in the runtime, and records the catalog selection with each fit. The general Draft 02 database contains 333 catalogs for 231 identities. The PILS candidate contains 344 catalogs for 242 identities and applies only within 329150.112845–362899.138030 MHz. The notebook can select individually reviewed additions, exclude native entries, and accept new text catalogs.

## First fit

1. Open Google Colab and select **File → Upload notebook**. Choose the downloaded notebook and a standard CPU runtime.
2. Leave **Draft 02 general**, **Public archive to runtime**, and **Synthetic example** selected. Run steps 1 through 7. The observed units are MHz and kelvin.
3. For the small reference observation, expect ten assigned peaks and three fitted species. Save the results ZIP from step 7. It contains the input, catalog manifest, fit settings, tables and an offline spectrum plot.

The [practice upload](docs/assets/rotor-upload-example.zip) contains a fictional line list and a matching spectrum. The [guide](docs/index.html#catalogs) gives the fit settings and expected result. Custom uploads use JSON metadata with either an explicit Einstein A CSV or a supported fixed-width SPCAT file. The partition function and line weights must use the same convention.

## Source and reproducibility

The upstream code base for this branch is `zfried/astro_amase` revision `651a2e47ffe1acaca2c139a9bc963779bf9fbbcc`. The source bundle embedded in the notebook is checked against the branch's analysis source. Saved results identify the source bundle, base data, active catalogs and any uploaded source files by checksum. Re-running from the saved settings requires those same versions.

The catalog membership and the acceptance of an astronomical detection are separate scientific decisions. Review alternative carriers, blending, missing predicted lines, and catalog coverage when interpreting assignments. The Lite changes and the upstream data sources are documented in the [notices](LITE-NOTICES.md). Keep the applicable license texts and per-catalog scientific credits with any shared work.

## Checking a change

Run `python3 -m unittest discover -s tests -p 'test_lite*.py' -v` from the repository root. This checks that the notebook cells compile and have cleared outputs, its embedded analysis source matches the branch, the catalog recipe checksums agree, and the offline instructions link to the included files. The same checks are configured for GitHub Actions. Changes to catalog membership, intensities or fitting also need numerical comparisons on representative spectra.
