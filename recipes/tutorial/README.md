## Tutorial

A recipe consists of two files: 

1. The interface specification. This is a text file in the JSON format that describes the parameters
and how these are rendered.
2. The script template. This specify the actions that the script performs.


Example:

* `hello1.hjson` is the inteface file
* `hello1.sh` is the script template for the first recipe

The interface file `hello1.hjson` contains the bare minimum information
needed to provide a simple context for the tool. This recipe does not take
any external parameters. 

    
    {
      settings: {
        name: Hello World
        summary: This recipe prints: Hello World!
        help:
          '''
          # Help
        
          This recipe prints **Hello World!** to the console.
        
          Note how the console output is captured in the standard output of the results.
          '''
      }
    }


Whereas the `hello1.sh` contains:
    
    # This is a regular bash script.
    
    echo 'Hello World!'

