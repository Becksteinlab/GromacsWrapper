#!/usr/bin/env tclsh
# $Id: update_package.tcl 3905 2009-08-24 14:25:52Z oliver $
#set VMDSCRIPTDIR $env(HOME)/Biop/Library/vmd/scripts
set VMDSCRIPTDIR $env(PWD)
cd $VMDSCRIPTDIR
pkg_mkIndex -verbose $VMDSCRIPTDIR 
 
