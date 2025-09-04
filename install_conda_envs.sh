# in your project folder
for yml in conda_setup/envs/*.yml; do
    echo "Creating conda env from $yml"
    conda env create -f $yml || conda env update -f $yml
done