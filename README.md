# AstroAmase Lite

Fit a radio astronomy spectrum in your own Google Colab session. Colab installs the analysis software and retrieves the selected molecular catalogs when you run the notebook.

**[Open AstroAmase Lite in Google Colab](https://colab.research.google.com/github/bmcguir2/astro_amase/blob/astroamase-lite/docs/downloads/AstroAmase_Lite_Integrated.ipynb)** · **[Read the step-by-step instructions](https://bmcguir2.github.io/astro_amase/)** · [Download the notebook](https://bmcguir2.github.io/astro_amase/downloads/AstroAmase_Lite_Integrated.ipynb)

## Run your first fit

1. Open the notebook in Colab and sign in to Google. Select a standard CPU runtime.
2. Keep **Draft 02 general**, **Public archive to runtime**, and **Synthetic example** selected. Click the run button beside each numbered step, from 1 through 7. On the first run, catalog preparation can take 10–15 minutes.
3. The practice spectrum should yield **10 assigned peaks and 3 fitted species**. In step 7, download the result ZIP to keep the fit, plots, settings, and catalog list.
4. For your own fit, upload a plain-text file with two columns, frequency in MHz and intensity in K. Set your fit parameters in the notebook. [Input format and controls](https://bmcguir2.github.io/astro_amase/#spectrum)

The opening notebook contains the 333-catalog general reference. The guide also lets you open the general ground-state and PILS-band editions directly in Colab. [Choose a notebook](https://bmcguir2.github.io/astro_amase/#notebook-editions)

The source in this repository supports the Colab workflow. [Software license and spectroscopy credits](LITE-NOTICES.md)
