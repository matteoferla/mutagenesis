This folder has the original C++ code written by Andrew Firth that can
be found at [the Pedel/Driver/Glue webpage](http://guinevere.otago.ac.nz/aef/STATS/)

They all compile fine, except for two `driver.cxx` and `glue_mc.cxx`

    $ g++ driver.cxx -o driver
    $ ./driver
    Segmentation fault: 11
    
    $ g++ driver.cxx -o driver0x -std=c++0x
    $ ./driver0x
    Segmentation fault: 11
    $ g++ driver.cxx -o driver98 -std=c++98
    $ ./driver98
    Segmentation fault: 11

Which means there is a stackOverflow, but I am unfamiliar with C++.