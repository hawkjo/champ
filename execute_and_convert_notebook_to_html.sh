out_fpath="~/Downloads/`basename $1 .ipynb`.html"
ipython nbconvert --to=html --ExecutePreprocessor.enabled=True --ClearOutputsPreprocessor.enabled=True --output=$out_fpath $1
