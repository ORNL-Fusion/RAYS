# A collection of examples running the RAYS code, the post_process_RAYS code,
and ray plotting codes

The RAYS code requires one input namelist file called rays.in.  The 
post\_process\_RAYS code also requires an input namelist file called 
post\_process\_rays.in.  Typically I give the input files a more specific name
and do a soft link to the generic filenames.  The rays.in list contains a variable
\'run\_label\' to help identify the specific run.  This run label is attached to the
names of output files again to help keep things straight and to make it convenient to
group related runs in a single directory.

I make myself a directory for RAYS runs outside github tree and run from subdirectories
in there.  I make an aliases for the codes to make execution simple. e.g.


alias rays=\'/Users/dbb/Lab\_work/git\_dbb/RAYS/RAYS\_code/RAYS\'

alias pp\_rays=\'/Users/dbb/Lab\_work/git\_dbb/RAYS/post\_process\_RAYS/post\_process\_RAYS\'

alias plt\_rays=\'python /Users/dbb/Lab\_work/git\_dbb/RAYS/graphics\_RAYS/plot\_RAYS\_slab.py\'

