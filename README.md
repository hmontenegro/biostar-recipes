# Biostar Recipes

Life sciences are experiencing a reproducibility crisis. 
It has become unexpectedly challenging to even understand how analyses 
are performed and even more difficult to repeat the same process again.

The Biostar Recipes were created to help scientists share, document, 
distribute and rerun data analysis scripts. 

We call these scripts recipes. 

Recipes are designed to facilitate 
the distribution, the sharing and the reuse of bioinformatics pipelines
Visit the [Biostar Engine][engine] for installation of each recipe into into the engine.

## Installation

#### 1\. Create a virtual environment

[conda]: https://conda.io/docs/

    conda create -y --name engine python=3.6
    source activate engine

#### 2\. Clone the recipe code:

There are different repositories for the engine and the recipes.

    # This repository stores the various data analysis recipes.
    git clone https://github.com/biostars/biostar-recipes.git

### 3\. Install the requirements:

    # Install the conda requirements.
    conda install --file conf/conda_requirements.txt

    # Add the recipes to the python path.
    python setup.py develop
    
## What is a recipe made of?

A recipe consists of two files: 

1. The interface specification. This is a text file in the JSON format that describes the parameters
and how these are rendered.
2. The script template. This specify the actions that the script performs.

## How to write recipes?

See the page 

* [How to write recipes?](docs/how-to-write-recipes.md)

## How do I learn how to write recipes?

Investigate the tutorial recipes at:

* https://github.com/biostars/biostar-recipes/tree/master/recipes/tutorial
  
Alternatively look at these same recipes deployed on the main site at:

* https://www.bioinformatics.recipes/project/view/tutorial/
    
## What are data collections?
 
The Biostar Engine operates with the concept of a *data collection*.
 
A collection may though of as a directory that contains one or more (any number) of files.

The `value` attribute of the data provides the first file of the collection.
This is handy when the collection contains a single file. The `toc` (table of contents) attribute 
is a filelist of all files in the collection. It is a list of all files in the collection.

Usage:

    echo parameter.value
    
to get the first file of the collection. Use
    
    cat parameter.toc
      
to obtain a list of all files of the collection.
    
## What is a data type?

Types can be thought of as tags that allow you to filter the input data.
For example a parameter of `FASTA` type will only list data that is tagged
with the word `FASTA`. Data may have more than one tag listed in the type.
A case insensitive regular expression match is performed to 
match the parameter type to the data type.

## Recipe requirements

We recommend the following best practices:

* Every recipe must be documented and fully operational.
* Every recipe must have test data and results associated with it to
demonstrate the input requirements as well as the results.

The test data should be small so that it can be readily 
executed to allow users to investigate the outputs of it.
