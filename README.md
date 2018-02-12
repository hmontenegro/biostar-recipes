# Biostar Recipes


Life sciences are experiencing a reproducibility crisis. 
It has become unexpectedly challenging to even understand how analyses 
are performed and even more difficult to repeat the same process again.

This [Biostar Engine][engine] was designed to help scientists share, document, 
distribute and rerun data analysis scripts. We call these scripts recipes. 

Recipes are designed to facilitate 
the distribution, the sharing and the reuse of bioinformatics pipelines

[engine]: https://github.com/biostars/biostar-engine

## Installation

A typical installation links the current package into the current python environment

    python setup.py develop
    

## What is a recipe made of?

A recipe consists of two files: 

1. The interface specification. This is a text file in the JSON format that describes the parameters
and how these are rendered.
2. The script template. This specify the actions that the script performs.

# What are collections?
 
The Biostar Engine operates with the concept of a *data collection*.
 
A collection may though of as a directory that contains one or more (any number) of files.

The `value` attribute of the data provides the first file of the collection 
(this is handy when the collection contains a single file) the `toc` (table of contents) attribute 
provides access to all files in the collection.

## Recipe requirements

* Every recipe must be documented and fully operational.
* Every recipe must have test data and results associated with it to
demonstrate the input requirements as well as the results.
* Every recipe may be modified, changed, reused and 
shared across different projects.
            
Notably the software is able to generate a graphical user interface
from each recipe when these contain a JSON specification syntax.
