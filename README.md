# Biostar Recipes


Life sciences are experiencing a reproducibility crisis. 
It has become unexpectedly challenging to even understand how analyses 
are performed and even more difficult to repeat the same process again.

The Biostar Recipes were created to help scientists share, document, 
distribute and rerun data analysis scripts. 

We call these scripts recipes. 

Recipes are designed to facilitate 
the distribution, the sharing and the reuse of bioinformatics pipelines

Visit the [Biostar Engine][engine] for installation of each recipe into  
into the engine.

    https://github.com/biostars/biostar-engine

[engine]: https://github.com/biostars/biostar-engine

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
