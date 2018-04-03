# IDL Scripts for IRTF Time Series Analysis

## Dependencies:
### BJD UTC Tools. 
You need to create an environmental variable `$ASTRO_DATA` and where you will put the leap year info:

	wget http://www.physics.wisc.edu/~craigm/idl/down/JPLEPH.405
	
(as in the BJD tools here: <a href=http://astroutils.astronomy.ohio-state.edu/time/pro/README>http://astroutils.astronomy.ohio-state.edu/time/pro/README</a>)


## Troubleshooting:

 - If you have a problem running `get_profile_widths`, it's probably the type of files in the file lists. You may want to run `get_profile_widths,/esXtract` to use the `_es_ms.fits` files.
