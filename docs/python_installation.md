# Python installation notes

- Install Anaconda
    - https://www.continuum.io/downloads
      - Consider `miniconda` if you don't need all the dev tools in the full `anaconda` distribution.
      - Recommended to download python v3.6 as your default distribution.

      - *DO NOT INSTALL AS ROOT OR ADMIN!*
- Modify path as necessary to point to Anaconda distribution (if not using pyenv)

- Installation of `pyenv` is recommended using `homebrew`:
```
  brew update
  brew install pyenv
```

- If you use `miniconda`, then install critical parts of python:
```
    conda install ipython
    conda install jupyter    
```

- Updates/installation of fundamental packages:
```
    conda install scipy
    conda install numpy
    conda install astropy
```

- pyND (this distribution; be sure to put it somewhere in your `PYTHONPATH`)
```
    git clone https://github.com/jchowk/pyND.git
```

- linetools
```
    git clone https://github.com/linetools/linetools.git
    cd linetools
    python setup.py develop
```

- pyigm
```
      git clone https://github.com/pyigm/pyigm.git
      cd pyigm
      python setup.py develop
```

- specdb
```    
    git clone https://github.com/specdb/specdb.git
    cd specdb
    python setup.py develop
```
    - grab the database (see installation notes in `igmspec`)


- astroquery  [required for finder charts]

     `pip install astroquery`

- emcee and corner:      
```
  pip install emcee  
  pip install corner```

- Install non-standard packages
    - aplpy  [recommended; not essential]

      `pip install aplpy`

    - pyregion [maybe?]

      `pip install pyregion`

    - pymc3 [only required for some f(N) analysis in pyigm]
```  
        git clone https://github.com/pymc-devs/pymc3.git
        cd pymc3
        python setup.py install```


- For even more specific cases:
  - fitsio [only for DESI codes]
    - Grab from github (my fork)
    - git clone https://github.com/profxj/fitsio.git
    - cd fitsio
    - python setup.py install
  - matplotlib Basemap  [For projection plots only]
    - conda install basemap
  - seaborn [Statistical plots; optional, but recommended]
      - conda install seaborn
  - bokeh [java figures; optional, but recommended]
      - conda install bokeh
  - ginga [Excellent imaging package; recommended]
      - git clone https://github.com/ejeschke/ginga.git
      - cd ginga
      - python setup.py install
  - PYPIT
    - See installation notes
    - https://github.com/PYPIT/PYPIT/blob/master/doc/installing.rst
  - pandoc — For notebook conversions to rst [optional]
    - https://github.com/jgm/pandoc/releases/tag/1.15.2
  - scikit

   - conda install scikit-learn
   - conda install scikit-image   #  For images


  - Progressbar (nice progress bars for long loops)
    - conda install progressbar

  - alis
    - git clone https://github.com/rcooke-ast/ALIS
    - ...

On rebuilds of astropy,
`git clean -fxd` or `sudo git clean -fxd`
is your friend.


**WARNING:**   `sudo python` may give you a version of Python that  you aren’t expecting!
