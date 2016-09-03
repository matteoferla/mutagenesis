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

Which means there is a stackOverflow, but I am unfamiliar with C++ and I am 
confused why a decade old script fails due to it.
The offending line is:

    define maxndaugh 524288
    
    double          PBk[maxndaugh],
                    Xk[maxndaugh],
                    Xk2[maxndaugh];

Namely, the code tries to make an array too large.
SO has this to say about this SO error: [link](StackOverflow has this to say about the stackOverflow error: [link]()http://stackoverflow.com/questions/571945/getting-a-stack-overflow-exception-when-declaring-a-large-array .) .

For now the following has been changed:

    #define maxpos 15
    #define maxndaugh 32768

In glue I changed:

    #define maxvar 1000000