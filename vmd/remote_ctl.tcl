# $Id$
# vmd remote control --- client/server scripts to run VMD from python
# Copyright (c) 2007-2009 Oliver Beckstein <orbeckst@gmail.com>
# Released under the GNU Lesser Public License, version 3 or later.
#
# based on bounce.tcl and vmdcollab.tcl
# from http://www.ks.uiuc.edu/Research/vmd/script_library/scripts/vmdcollab/
# with help from Justin Gullingsrud <justin@ks.uiuc.edu>
#
# Start this script in VMD and send commands to the listening port to have VMD 
# execute them remotely (eg 'telnet localhost 5555'). However, this is really 
# meant to be used with vmd_ctl.py. That is why we have a special 
# RESPONSE_TERMINATOR to ease the communication (see Protocol)
#
# Usage: vmd -e remote_ctl.tcl
#    or  vmd> source remote_ctl.tcl
#
# Security: we only allow connections from localhost (see acpt)
#  
# Protocol: port 5555
#           RESPONSE_TERMINATOR to signify end of VMDs response
#
# Commands: some transmitted commands have special meaning and are interpreted
#           by the server and not passed on to Tcl:
#
#           close             close the current connection (socket)
#           exit              exit the server (no more connections possible,
#                             but current connections are still open)
#           loglevel N        set LOGLEVEL to value N (0<=N<=2)
  
namespace eval remote_ctl {
    variable main
    variable clients
    variable default_vmd_port 5555

    variable RESPONSE_TERMINATOR "__END_OF_VMD_RESPONSE__\n"
    # important: puts translates '\n' -> '\r\n' (which MUST be used in vmd_ctl.py)

    variable LOGLEVEL 1
    # loglevel: 0 quiet
    #           1 startup and quit message, unauthorised connections
    #           2 socket closures

    proc start "{port $default_vmd_port}" {
	variable main
	set main [socket -server remote_ctl::acpt $port]  
	putlog "Listening on port $port" 1
    }

    proc acpt { sock addr port } {
	variable clients
	if {[string compare $addr "127.0.0.1"] != 0} {
	    putlog "Unauthorized connection attempt from $addr port $port" 1
	    close $sock
	    return
	}
	putlog "Accept $sock from $addr port $port" 2
	set clients($sock) 1
	fconfigure $sock -buffering line
	fileevent $sock readable [list remote_ctl::recv $sock]
    }

    proc recv { sock } {
	variable main
	variable clients
	variable RESPONSE_TERMINATOR
	variable LOGLEVEL
	if { [eof $sock] || [catch {gets $sock line}] || \
		 [string compare $line "close"] == 0} {
	    # end of file, abnormal connection drop, or requested:
	    # shut down this connection
	    close $sock
	    putlog "Closing $sock" 2
	    unset clients($sock)
	} else {
	    if {[string compare $line "exit"] == 0} {
		# prevent new connections
		# existing connections stay open
		putlog "Disallowing further incoming connections by request of $sock" 1
		close $main
		return
	    } elseif {[regexp {loglevel *([0-9])} $line dummy LOGLEVEL]} {
		putlog "Changed loglevel to $LOGLEVEL" 1
		return
	    }
		
	    # execute the received commands and send back result
	    set rc [catch $line result]
	    if { $rc } {
		puts $sock "Error executing command '$line': \n$result"
		puts       "Error executing command '$line': \n$result"
	    } else {
		puts $sock $result
		puts       $result
	    }
	    # terminator -- must correspond to RESPONSE_TERMINATOR in vmd_ctl.py
	    puts -nonewline $sock $RESPONSE_TERMINATOR
	}
    }
    
    proc putlog { text {threshold 0} } {
	variable LOGLEVEL
	if {[expr $LOGLEVEL >= $threshold]}  {
	    puts $text
	}
	return
    }
}



# assuming that loading this file means that we actually want to run
# the server, we just  start it:
remote_ctl::putlog "Starting remote_ctl server in VMD."
remote_ctl::start
