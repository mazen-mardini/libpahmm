setup.py uses a library called pcpp to parse the preprosessor 
directives in include/cpahmm.h. The results are then fed into
cffi to compile the foreign interface layer of
the library. Sometimes external dependencies are
included using "#include", these must be ignored during the
parsing process, and the way we do this is by providing fake 
files to silence the parser, effectively making it omit any 
"#include" directive.

If you include any external header-files in include/cpahmm.h, 
it is very important that you also create an empty dummy-file 
with the same name in this directory.