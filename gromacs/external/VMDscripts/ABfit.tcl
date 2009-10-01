# $Id$
# Copyright (c) 2009 Oliver Beckstein
# Released under GNU Public License 3 (or higher, your choice)
#
# Align a trajectory A to a reference structure B.
# Based on suggestions by Philip Fowler and Ben Hall.
#
# Cope with big trajectories by writing dcds in chunks and
# joint them later with catdcd.

package provide ABfit 0.1

package require bigdcd_GW

proc _rmsfit {frame} {
    global ref sel all
    $all move [measure fit $sel $ref]
    puts "$frame: [measure rmsd $sel $ref]"
}

proc ABfit {args} {
    global ref sel all
    set usage "usage: ABsetup -ref <atomselection>  -sel <atomselection> -all <atomselection>
                              -traj <list of input trajectories>  -type <filetype>
                              -step <int> -prefix <string>

Align *sel* on *ref* and write a dcd of *all*. If *all* is not
supplied then it defaults to everything in the top molecule.

*traj* is a list of trajectory files that are loaded sequentially; the
type---dcd, xtc, ...---can be set with the *type* keyword. The
trajectories are loaded into molecule top.

*prefix* is the prefix for the output dcd files; they are numbered
sequentially as <prefix>_0001.dcd, <prefix>_0002.dcd, ..."

    set ref "NONE"
    set sel "NONE"
    set all [atomselect top {all}]

    set prefix "ABfit_"
    set type "auto"
    set step 1
    set trajectories {}

    foreach {argname argval} $args {
	switch $argname {
	    "-ref"    {set ref $argval}
	    "-sel"    {set sel $argval}
	    "-all"    {set all $argval}
	    "-prefix" {set prefix $argval}
	    "-type"   {set type $argval}
	    "-step"   {set step $argval}
	    "-traj"   {set trajectories $argval}   
	    default  {puts "ABfit_setup: Unknown argument $argname!\n$usage"
		return}
	}
    }

    if {$ref == "NONE"} {
	puts "Set ref to the reference atomselection"
	return
    }
    if {$sel == "NONE"} {
	puts "Set sel to the selection atomselection"
	return
    }

    bigdcd2dcd _rmsfit $type $step $prefix $trajectories
    
}
