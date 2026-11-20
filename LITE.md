# AstroAmase Lite

Open the [Colab quick start](README.md) or the [instructions website](https://bmcguir2.github.io/astro_amase/) to choose a notebook and run a fit.

## Maintainer checks

Run `python3 -m unittest discover -s tests -p 'test_lite*.py' -v` from the repository root. The checks verify the four downloadable notebooks, the embedded analysis source, the catalog recipe checksums, and the instructions site. GitHub Actions runs the same checks. Changes to catalog selection, line intensities or fitting also require comparisons on representative spectra.

For published fits, record the notebook edition and catalog provenance from the result ZIP. See the [software license and spectroscopy credits](LITE-NOTICES.md).
