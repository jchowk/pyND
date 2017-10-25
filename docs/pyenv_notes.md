# Using `pyenv`

The routine `pyenv` gives you some options for "easily" using multiple python versions.
## I. Installation
### A. Install with homebrew
```
brew install pyenv
brew install pyenv-virtualenv
```
### B. Ensure your BASH `.profile` is set up to use `pyenv`
Insert the following into your `.profile` near the bottom.
```
if command -v pyenv 1>/dev/null 2>&1; then
    eval "$(pyenv init -)"
    eval "$(pyenv virtualenv-init -)"
fi
```
## II. Initialization
In order to define your root or system version of python.

```
pyenv init
```



## III. Installing python versions
I used the following to install python v2.7 using the *miniconda* distribution.

```
pyenv install miniconda2-latest
```

## IV. Switching to the new version

There are several ways to do this. These include:

1. `pyenv activate miniconda2-latest` —  This activates my python 2.7 distribution, which I downloaded using pyenv as `miniconda2-latest`. To get rid of it do `pyenv deactivate`.

2. `pyenv shell miniconda2-latest` — This changes python versions for your current shell. As far as I know you have to run it again with your root distribution name to go back.

3. `pyenv local miniconda2-latest` — This is slick. If you have a directory for which you always want to run python 2.7, you can do this in that directory. It then forces it to change python versions to the noted one everytime you run something from this directory.
