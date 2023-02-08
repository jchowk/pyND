# pyND: Installing / using `pyND`

There are several ways to make use of `pyND`. Currently the principal approaches are:

## **Editable pip installation:** 

This installation approach allows edits to the code / `git pull` updates to be directly accessible. To enable this installation, invoke the following in the terminal from within the `pyND` directory:

```
pip install -e .
```

## **Install from GitHub via `pip`:**

```
pip install git+https://github.com/jchowk/pyND.git
```

## **Include `pyND` in your `$PYTHONPATH`:**

This approach also makes changes to the base code immediately accessible. Add the full path to the `pyND` code to your `$PYTHONPATH` variable by invoking something like (from the `bash` terminal or within your shell configuration file):

```
export PYTHONPATH="$PYTHONPATH:/path/to/pyND/pyND/"
```

Note the path has to point to the subdirectory `pyND/pyND/`. 
