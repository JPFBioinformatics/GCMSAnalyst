To convert from .ms to .mzML use proteowizard
    there is also mzXML, but this is a dated file format and only to be used if you are using
    legacy software
    proteowizard was downloaded from github and is written in C++ so you have to use VScode
    to write it to a library so you can use it in the project
    I do not have the C++ knowledge to do this, so instead we will download the end user app
    version of proteowizard and call it as a subprocess of our application
        this will slightly reduce efficiency but due to small data load this will likely not
        be significant

