# Reference Data

A list of data  that need to be installed 

### Directory setup

    # This will store reference datasets.
    mkdir -p /export/refs
    
    # This will store source code not installed through conda.
    mkdir -p /export/src
    
    # Link executables into this directory.
    mkdir -p /export/bin
    
### Centrifuge

    DIR=/export/refs/centrifuge/
    mkdir -p $DIR
    (cd $DIR && curl -O ftp://ftp.ccb.jhu.edu/pub/infphilo/centrifuge/data/p_compressed+h+v.tar.gz)
    (cd $DIR && tar xzvf p_compressed+h+v.tar.gz)
    
### Augustus

    DIR=/export/src/
    mkdir -p $DIR
    (cd $DIR && curl -O http://bioinf.uni-greifswald.de/augustus/binaries/augustus.current.tar.gz)
    (cd $DIR && tar xzvf augustus.current.tar.gz)
    sudo apt-get install bamtools libbamtools-dev
    (cd $DIR/augustus && make)
    
then one of:
    
    #
    #(cd $DIR/augustus && sudo make install)
    #
    # Alternatively link each required executable.
    #ln -fs $DIR/augustus/bin/augustus /export/bin
    # ln -fs $DIR/augustus/bin/augustus /export/bin
