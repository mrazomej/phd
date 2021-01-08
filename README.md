# `phd`

<a rel="license" href="http://creativecommons.org/licenses/by/4.0/"><img alt="Creative Commons License" style="border-width:0" src="https://i.creativecommons.org/l/by/4.0/88x31.png" /></a><br />This work is licensed under a <a rel="license" href="http://creativecommons.org/licenses/by/4.0/">Creative Commons Attribution 4.0 International License</a>.

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.3815142.svg)](https://doi.org/10.5281/zenodo.3815142)


## Contact
If you have any questions about the scientific material or the actual structure
of this thesis repository please [feel free to open an issue](https://github.com/gchure/phd/issues) and I will
respond as soon as I can.

Each chapter of the thesis has a [Disqus](https://disqus.com/) discussion board at the bottom of
[each Chapter](https://gchure.github.io/phd/chapter_01/#references). If you'd like to have a more public conversation, please head there!

## About
This repository contains all files used to generate my PhD thesis "The Molecular
Biophysics of Evolutionary and Physiological Adaptation" defended at the
California Institute of Technology (hopefully) on June 1, 2020. This includes
all Python and Stan code used to generate figures, all markdown files used to
actually write the thesis, and all HTML files and shell scripts used to generate
the website and PDF versions of the thesis. 

## Repository Structure

### `docs`
This folder contains all live files displayed on the associated thesis website
`http://gchure.github.io/phd`. This folder is recognized by GitHub Pages which
serves the web content. Everything in this folder is generated by executing
`make html` in the repository root. There is nothing of *direct* use in this
folder of the repository. 

### `doks-theme`
This folder contains all materials for the [Doks website theme](https://doks.themejack.com/blue/) which is
parsed and built using the [Jekyll](https://jekyllrb.com/) blogging framework. This theme was
purchased and subsequently *heavily* modified to match my aesthetic
requirements. If you are interested in using a similar theme for your project,
please support the developers of the theme and [purchase a license](https://themeforest.net/item/doks-jekyll-theme-for-project-documentation/21102900).

### `src`
This folder houses and meat-and-potatoes of the actual PhD thesis. There are
many subdirectories, which I outline below. 

* **`chapter_0X`** \| There is a single folder for each chapter of the thesis.
  Within each folder is a single markdown document for each portion of the
  chapter (e.g. `section_01_abstract.md`, `section_02_model.md` etc.). These
  files contain the text of the chapter.
  + `figs` \| This folder contains PDF and PNG versions of all images used in that chapter. 
  + `code` \| This folder contains all Python and/or Stan files used in the
    generation of all quantitative figures of the chapter. Each script relies on
    the use of the `phd` custom software module which is described in the next
    section of this README file.
* **`frontmatter`** \| This folder contains all of the thesis frontmatter
  markdown files including the abstract, acknowledgments, list of published
  contributions, and a `lua` filter that allows the use of short-captions for
  the PDF version of the thesis. 
* **`styles`** \| This folder contains the bibliography format used in this
  thesis (`cell.csl`) as well as a [Pandoc markdown
  template](https://pandoc.org/MANUAL.html#templates) which follows the
  guidelines provided by the [Caltech Thesis
  website](https://libguides.caltech.edu/theses).

Outside of these directories, there is a smattering of other files that are used
in the actual generation of the thesis. These are listed below. 

* **`caltech_thesis.cls`** \| This is the LaTeX template file provided from the
  Caltech library.
* **`caltech.png`** \| A high-resolution PNG file of the Caltech logo which is
  present on the PDF version of the thesis.
* **`metadata.yaml`** \| A yaml file containing metadata about the thesis. This
  includes author name, thesis title, rights statement, and a variety of other
  configuration details for the HTML and PDF generation of the thesis. 
* **`ref_formatting.yaml`** \| This is a yaml file which contains configurations
  for how figure, equation, and template references should be displayed. 
* **`ref_heading.md`** \| A markdown file with a single line that reads `#
  References`. This is a hack I used to allow chapter-level references for the
  website, but not for the PDF.
* **`references.bib`** \| The bibliography file containing information of all
  references used in the thesis. 
* **`simple_letters.ist`** \| A index style file provided by the Caltech
  Library. 

To execute much of the code in the `chapter_0X/code` directories, you must have
the data sets downloaded locally to a `src/data` folder. Due to size
constraints, these data have not been directly stored on GitHub and you must
download them externally. These data are available via the CaltechDATA research
data repository under the [DOI: 10.22002/D1.1426](https://doi.org/10.22002/D1.1426).

### `phd`
This is a custom Python package used in the analysis and generation of all
figures presented in the thesis. If you want to regenerate any of the figures in
this thesis, you will need to install this package locally. You can do so by
executing the following command of the root directory of the repository. 

```
pip install -e ./
```

Within the `phd` folder is a slew of Python modules. These modules serve the
following functions. 

* **`__init__.py`** \| An init file used to collect the package for import. 
* **`_fit_bivariate_normal_AstroML.py`** \| This contains several functions
  adapted from [AstroPy](https://www.astropy.org/) to facilitate fitting of a bivariate gaussian
  distribution to flow cytometry data for reproducible gating. 
* **`bayes.py`** \| Contains myriad functions for performing Bayesian
  statistical inference including interaction with the
  [Stan](http://mc-stan.org) probabalistic programming language. 
* **`flow.py`** \| A variety of functions for interacting with flow cytometry
  data. 
* **`image.py`** \| A variety of functions for processing image files including
  segmentation, filtering, and measurement of object properties. 
* **`io.py`** \| A slew of functions for file IO.
* **`mscl.py`** \| A set of highly-custom functions used in the analysis of
  single-cell microscopy data discussed in Chapters 5 and 9 of the thesis. 
* **`stats.py`** \| Helper functions for easy calculation of statistical
  properties including Highest Probability Density intervals which are reported
  frequently throughout the thesis. 
* **`thermo.py`** \| Contains two classes `MWC` and `SimpleRepression` which
  define the thermodynamic models discussed in the bulk of the thesis. 
* **`viz.py`** \| A collection of functions which format [matplotlib](https://matplotlib.org/) and
  [Altair](https://altair-viz.github.io/) plots to my personal style.

### `dst`
This directory contains the PDF version of the thesis and is populated by
executing `make pdf` from the root directory. This folder exists solely so I
didn't confuse myself on where the most recent version of the PDF lived. 

## `./`
Outside of the directories described above, there are a variety of other files.
Below, I describe the purpose of each one. 

* **`_config.yml`** \| A configuration file for website generation via
  [Jekyll](https://jekyllrb.com/).
* **`.gitignore`** \| A text file with information of files that Git should ignore.
* **`abstract.md`** \| A markdown file that specifies the [abstract
  page](http://gchure.github.io/phd/abstract) of the thesis website. 
* **`acknowledgements.md`** \| Just like `abstract.md`, this file generates the
  webpage for the
  [acknowledgments](http://gchure.github.io/phd/acknowledgements)
* **`chapter_0X.md`** \| Markdown files that contain all information necessary
  to build the chapter pages of the thesis website. 
* **`coppyfigs.sh`** \| A shell script to be executed by the `makefile` which
  copies all `.png` versions of the figures from the `src/chapter_0X/figs/`
  folder and copies them to the appropriate place in the `docs/` folder. 
* **`favicon.ico, Gemfile, Gemfile.lock`** \| Files which specify details about
  the Ruby requirements necessary to build the website. 
* **`index.md`** \| The [homepage](http://gchure.github.io/phd) of the thesis
  website. 
* **`Makefile`** \| A makefile used to assemble the thesis in to a
  human-friendly form. This file can be executed to make either the website
  components (`make html`) or the PDF components (`make pdf`) of the thesis. 
* **`setup.py`** \| Python file necessary for installation of the `phd` Python
  module. 


# License
All creative work in this repository is licensed under a [Creative Commons CC-BY
4.0 Attribution license](https://creativecommons.org/licenses/by/4.0/). All
software is issued under a standard MIT license as follows:

```
Copyright 2021 Manuel Razo-Mejia

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
```

